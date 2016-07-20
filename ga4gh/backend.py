"""
Module responsible for handling protocol requests and returning
responses.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol
import ga4gh.registry as registry

# We must import these modules before queries are run, or we won't be
# able to load the polymorphic types at run time. These should
# probably be imported in the frontend so that we are more
# configurable and open the possibility of a plugin architecture.
import ga4gh.datasource.simulator  # noqa
import ga4gh.datasource.htslib  # noqa
import ga4gh.datasource.sql  # noqa
import ga4gh.datasource.obo  # noqa


def _parseIntegerArgument(args, key, defaultValue):
    """
    Attempts to parse the specified key in the specified argument
    dictionary into an integer. If the argument cannot be parsed,
    raises a BadRequestIntegerException. If the key is not present,
    return the specified default value.
    """
    ret = defaultValue
    try:
        if key in args:
            try:
                ret = int(args[key])
            except ValueError:
                raise exceptions.BadRequestIntegerException(key, args[key])
    except TypeError:
        raise Exception((key, args))
    return ret


def _parsePageToken(pageToken, numValues):
    """
    Parses the specified pageToken and returns a list of the specified
    number of values. Page tokens are assumed to consist of a fixed
    number of integers seperated by colons. If the page token does
    not conform to this specification, raise a InvalidPageToken
    exception.
    """
    tokens = pageToken.split(":")
    if len(tokens) != numValues:
        msg = "Invalid number of values in page token"
        raise exceptions.BadPageTokenException(msg)
    try:
        values = map(int, tokens)
    except ValueError:
        msg = "Malformed integers in page token"
        raise exceptions.BadPageTokenException(msg)
    return values


class IntervalIterator(object):
    """
    Implements generator logic for types which accept a start/end
    range to search for the object. Returns an iterator over
    (object, pageToken) pairs. The pageToken is a string which allows
    us to pick up the iteration at any point, and is None for the last
    value in the iterator.
    """
    def __init__(self, request, parentContainer):
        self._request = request
        self._parentContainer = parentContainer
        self._searchIterator = None
        self._currentObject = None
        self._nextObject = None
        self._searchAnchor = None
        self._distanceFromAnchor = None
        if not request.page_token:
            self._initialiseIteration()
        else:
            # Set the search start point and the number of records to skip from
            # the page token.
            searchAnchor, objectsToSkip = _parsePageToken(
                request.page_token, 2)
            self._pickUpIteration(searchAnchor, objectsToSkip)

    def _extractProtocolObject(self, obj):
        """
        Returns the protocol object from the object passed back by iteration.
        """
        return obj

    def _initialiseIteration(self):
        """
        Starts a new iteration.
        """
        self._searchIterator = self._search(
            self._request.start,
            self._request.end if self._request.end != 0 else None)
        self._currentObject = next(self._searchIterator, None)
        if self._currentObject is not None:
            self._nextObject = next(self._searchIterator, None)
            self._searchAnchor = self._request.start
            self._distanceFromAnchor = 0
            firstObjectStart = self._getStart(self._currentObject)
            if firstObjectStart > self._request.start:
                self._searchAnchor = firstObjectStart

    def _pickUpIteration(self, searchAnchor, objectsToSkip):
        """
        Picks up iteration from a previously provided page token. There are two
        different phases here:
        1) We are iterating over the initial set of intervals in which start
        is < the search start coorindate.
        2) We are iterating over the remaining intervals in which start >= to
        the search start coordinate.
        """
        self._searchAnchor = searchAnchor
        self._distanceFromAnchor = objectsToSkip
        self._searchIterator = self._search(
            searchAnchor,
            self._request.end if self._request.end != 0 else None)
        obj = next(self._searchIterator)
        if searchAnchor == self._request.start:
            # This is the initial set of intervals, we just skip forward
            # objectsToSkip positions
            for _ in range(objectsToSkip):
                obj = next(self._searchIterator)
        else:
            # Now, we are past this initial set of intervals.
            # First, we need to skip forward over the intervals where
            # start < searchAnchor, as we've seen these already.
            while self._getStart(obj) < searchAnchor:
                obj = next(self._searchIterator)
            # Now, we skip over objectsToSkip objects such that
            # start == searchAnchor
            for _ in range(objectsToSkip):
                if self._getStart(obj) != searchAnchor:
                    raise exceptions.BadPageTokenException
                obj = next(self._searchIterator)
        self._currentObject = obj
        self._nextObject = next(self._searchIterator, None)

    def next(self):
        """
        Returns the next (object, nextPageToken) pair.
        """
        if self._currentObject is None:
            raise StopIteration()
        nextPageToken = None
        if self._nextObject is not None:
            start = self._getStart(self._nextObject)
            # If start > the search anchor, move the search anchor. Otherwise,
            # increment the distance from the anchor.
            if start > self._searchAnchor:
                self._searchAnchor = start
                self._distanceFromAnchor = 0
            else:
                self._distanceFromAnchor += 1
            nextPageToken = "{}:{}".format(
                self._searchAnchor, self._distanceFromAnchor)
        ret = self._extractProtocolObject(self._currentObject), nextPageToken
        self._currentObject = self._nextObject
        self._nextObject = next(self._searchIterator, None)
        return ret

    def __iter__(self):
        return self


class ReadsIntervalIterator(IntervalIterator):
    """
    An interval iterator for reads
    """
    def __init__(self, request, parentContainer, reference):
        self._reference = reference
        super(ReadsIntervalIterator, self).__init__(request, parentContainer)

    def _search(self, start, end):
        return self._parentContainer.getReadAlignments(
            self._reference, start, end)

    @classmethod
    def _getStart(cls, readAlignment):
        if readAlignment.alignment.position.position == 0:
            # unmapped read with mapped mate; see SAM standard 2.4.1
            return readAlignment.next_mate_position.position
        else:
            # usual case
            return readAlignment.alignment.position.position

    @classmethod
    def _getEnd(cls, readAlignment):
        return (
            cls._getStart(readAlignment) +
            len(readAlignment.aligned_sequence))


class VariantsIntervalIterator(IntervalIterator):
    """
    An interval iterator for variants
    """

    def _search(self, start, end):
        return self._parentContainer.getVariants(
            self._request.reference_name, start, end,
            self._request.call_set_ids)

    @classmethod
    def _getStart(cls, variant):
        return variant.start

    @classmethod
    def _getEnd(cls, variant):
        return variant.end


class VariantAnnotationsIntervalIterator(IntervalIterator):
    """
    An interval iterator for annotations
    """

    def __init__(self, request, parentContainer):
        super(VariantAnnotationsIntervalIterator, self).__init__(
            request, parentContainer)
        # TODO do input validation somewhere more sensible
        if self._request.effects is None:
            self._effects = []
        else:
            self._effects = self._request.effects

    def _search(self, start, end):
        return self._parentContainer.getVariantAnnotations(
            self._request.reference_name, start, end)

    def _extractProtocolObject(self, pair):
        variant, annotation = pair
        return annotation

    @classmethod
    def _getStart(cls, pair):
        variant, annotation = pair
        return variant.start

    @classmethod
    def _getEnd(cls, pair):
        variant, annotation = pair
        return variant.end

    def next(self):
        while True:
            ret = super(VariantAnnotationsIntervalIterator, self).next()
            vann = ret[0]
            if self.filterVariantAnnotation(vann):
                return self._removeNonMatchingTranscriptEffects(vann), ret[1]
        return None

    def filterVariantAnnotation(self, vann):
        """
        Returns true when an annotation should be included.
        """
        # TODO reintroduce feature ID search
        ret = False
        if len(self._effects) != 0 and not vann.transcript_effects:
            return False
        elif len(self._effects) == 0:
            return True
        for teff in vann.transcript_effects:
            if self.filterEffect(teff):
                ret = True
        return ret

    def filterEffect(self, teff):
        """
        Returns true when any of the transcript effects
        are present in the request.
        """
        ret = False
        for effect in teff.effects:
            ret = self._matchAnyEffects(effect) or ret
        return ret

    def _checkIdEquality(self, requestedEffect, effect):
        """
        Tests whether a requested effect and an effect
        present in an annotation are equal.
        """
        return self._idPresent(requestedEffect) and (
            effect.id == requestedEffect.id)

    def _idPresent(self, requestedEffect):
        return requestedEffect.id != ""

    def _matchAnyEffects(self, effect):
        ret = False
        for requestedEffect in self._effects:
            ret = self._checkIdEquality(requestedEffect, effect) or ret
        return ret

    def _removeNonMatchingTranscriptEffects(self, ann):
        newTxE = []
        if len(self._effects) == 0:
            return ann
        for txe in ann.transcript_effects:
            add = False
            for effect in txe.effects:
                if self._matchAnyEffects(effect):
                    add = True
            if add:
                newTxE.append(txe)
        ann.transcript_effects.extend(newTxE)
        return ann


class Backend(object):
    """
    Backend for handling the server requests.
    This class provides methods for all of the GA4GH protocol end points.
    """
    def __init__(self, registry_db):
        self._requestValidation = False
        self._responseValidation = False
        self._defaultPageSize = 100
        self._maxResponseLength = 2**20  # 1 MiB
        self._registry_db = registry_db

    def get_registry_db(self):
        """
        Get the registry DB used by this backend
        """
        return self._registry_db

    def setRequestValidation(self, requestValidation):
        """
        Set enabling request validation
        """
        self._requestValidation = requestValidation

    def setResponseValidation(self, responseValidation):
        """
        Set enabling response validation
        """
        self._responseValidation = responseValidation

    def setDefaultPageSize(self, defaultPageSize):
        """
        Sets the default page size for request to the specified value.
        """
        self._defaultPageSize = defaultPageSize

    def setMaxResponseLength(self, maxResponseLength):
        """
        Sets the approximate maximum response length to the specified
        value.
        """
        self._maxResponseLength = maxResponseLength

    def startProfile(self):
        """
        Profiling hook. Called at the start of the runSearchRequest method
        and allows for detailed profiling of search performance.
        """
        pass

    def endProfile(self):
        """
        Profiling hook. Called at the end of the runSearchRequest method.
        """
        pass

    def validateRequest(self, jsonDict, requestClass):
        """
        Ensures the jsonDict corresponds to a valid instance of requestClass
        Throws an error if the data is invalid
        """
        if self._requestValidation:
            if not protocol.validate(jsonDict, requestClass):
                raise exceptions.RequestValidationFailureException(
                    jsonDict, requestClass)

    ###########################################################
    #
    # Iterators over the data hierarchy. These methods help to
    # implement the search endpoints by providing iterators
    # over the objects to be returned to the client.
    #
    ###########################################################

    def _topLevelObjectGenerator(self, request, numObjects, getByIndexMethod):
        """
        Returns a generator over the results for the specified request, which
        is over a set of objects of the specified size. The objects are
        returned by call to the specified method, which must take a single
        integer as an argument. The returned generator yields a sequence of
        (object, nextPageToken) pairs, which allows this iteration to be picked
        up at any point.
        """
        currentIndex = 0
        if request.page_token:
            currentIndex, = _parsePageToken(request.page_token, 1)
        while currentIndex < numObjects:
            object_ = getByIndexMethod(currentIndex)
            currentIndex += 1
            nextPageToken = None
            if currentIndex < numObjects:
                nextPageToken = str(currentIndex)
            yield object_.get_protobuf(), nextPageToken

    def _protocolObjectGenerator(self, request, numObjects, getByIndexMethod):
        """
        Returns a generator over the results for the specified request, from
        a set of protocol objects of the specified size. The objects are
        returned by call to the specified method, which must take a single
        integer as an argument. The returned generator yields a sequence of
        (object, nextPageToken) pairs, which allows this iteration to be picked
        up at any point.
        """
        currentIndex = 0
        if request.page_token:
            currentIndex, = _parsePageToken(request.page_token, 1)
        while currentIndex < numObjects:
            object_ = getByIndexMethod(currentIndex)
            currentIndex += 1
            nextPageToken = None
            if currentIndex < numObjects:
                nextPageToken = str(currentIndex)
            yield object_, nextPageToken

    def _protocolListGenerator(self, request, objectList):
        """
        Returns a generator over the objects in the specified list using
        _protocolObjectGenerator to generate page tokens.
        """
        return self._protocolObjectGenerator(
            request, len(objectList), lambda index: objectList[index])

    def _objectListGenerator(self, request, objectList):
        """
        Returns a generator over the objects in the specified list using
        _topLevelObjectGenerator to generate page tokens.
        """
        return self._topLevelObjectGenerator(
            request, len(objectList), lambda index: objectList[index])

    def _singleObjectGenerator(self, datamodelObject):
        """
        Returns a generator suitable for a search method in which the
        result set is a single object.
        """
        yield (datamodelObject.toProtocolElement(), None)

    def _noObjectGenerator(self):
        """
        Returns a generator yielding no results
        """
        return iter([])

    def _run_sql_query(self, request, query, responseBuilder):
        start = 0
        if request.page_token:
            try:
                start = int(request.page_token)
            except ValueError:
                raise exceptions.BadPageTokenException(request.page_token)
        end = start + request.page_size
        offset = start
        for result in query[start:end]:
            offset += 1
            responseBuilder.addValue(result.get_protobuf())
            if responseBuilder.isFull():
                break
        total_rows = query.count()
        if offset < total_rows:
            responseBuilder.setNextPageToken(str(offset))

    def datasetsGenerator(self, request, responseBuilder):
        """
        Returns a generator over the (dataset, nextPageToken) pairs
        defined by the specified request
        """
        query = self._registry_db.get_datasets_search_query(request)
        self._run_sql_query(request, query, responseBuilder)

    def bioSamplesGenerator(self, request, responseBuilder):
        # Ensure that the dataset exists.
        self._registry_db.get_dataset(request.dataset_id)
        query = self._registry_db.get_bio_samples_search_query(request)
        self._run_sql_query(request, query, responseBuilder)

        # dataset = self.getDataRepository().getDataset(request.dataset_id)
        # results = []
        # for obj in dataset.getBioSamples():
        #     include = True
        #     if request.name:
        #         if request.name != obj.getLocalId():
        #             include = False
        #     if request.individual_id:
        #         if request.individual_id != obj.getIndividualId():
        #             include = False
        #     if include:
        #         results.append(obj)
        # return self._objectListGenerator(request, results)

    def individualsGenerator(self, request, responseBuilder):
        # Ensure that the dataset exists.
        self._registry_db.get_dataset(request.dataset_id)
        query = self._registry_db.get_individuals_search_query(request)
        self._run_sql_query(request, query, responseBuilder)
        # dataset = self.getDataRepository().getDataset(request.dataset_id)
        # results = []
        # for obj in dataset.getIndividuals():
        #     include = True
        #     if request.name:
        #         if request.name != obj.getLocalId():
        #             include = False
        #     if include:
        #         results.append(obj)
        # return self._objectListGenerator(request, results)

    def readGroupSetsGenerator(self, request, responseBuilder):
        # """
        # Returns a generator over the (readGroupSet, nextPageToken) pairs
        # defined by the specified request.
        # """
        # dataset = self.getDataRepository().getDataset(request.dataset_id)
        # results = []
        # for obj in dataset.getReadGroupSets():
        #     include = True
        #     rgsp = obj.toProtocolElement()
        #     if request.name:
        #         if request.name != obj.getLocalId():
        #             include = False
        #     if request.bio_sample_id:
        #         rgsp.ClearField("read_groups")
        #         for readGroup in obj.getReadGroups():
        #             if request.bio_sample_id == readGroup.getBioSampleId():
        #                 rgsp.read_groups.extend(
        #                     [readGroup.toProtocolElement()])
        #         # If none of the biosamples match and the readgroupset
        #         # contains reagroups, don't include in the response
        #         if len(rgsp.read_groups) == 0 and \
        #                 len(obj.getReadGroups()) != 0:
        #             include = False
        #         else:
        #             include = True and include
        #     if include:
        #         results.append(rgsp)
        # return self._protocolListGenerator(request, results)
        """
        Returns a generator over the (readGroupSet, nextPageToken) pairs
        defined by the specified request.
        """
        # Ensure that the dataset exists.
        self._registry_db.get_dataset(request.dataset_id)
        query = self._registry_db.get_read_group_sets_search_query(request)
        self._run_sql_query(request, query, responseBuilder)

    def referenceSetsGenerator(self, request, responseBuilder):
        """
        Returns a generator over the (referenceSet, nextPageToken) pairs
        defined by the specified request.
        """
        query = self._registry_db.get_reference_sets_search_query(request)
        self._run_sql_query(request, query, responseBuilder)

    def referencesGenerator(self, request, responseBuilder):
        """
        Returns a generator over the (reference, nextPageToken) pairs
        defined by the specified request.
        """
        # Ensure that the reference set exists.
        self._registry_db.get_reference_set(request.reference_set_id)
        query = self._registry_db.get_references_search_query(request)
        self._run_sql_query(request, query, responseBuilder)

    def variantSetsGenerator(self, request, responseBuilder):
        """
        Returns a generator over the (variantSet, nextPageToken) pairs defined
        by the specified request.
        """
        # Ensure that the dataset exists.
        self._registry_db.get_dataset(request.dataset_id)
        query = self._registry_db.get_variant_sets_search_query(request)
        self._run_sql_query(request, query, responseBuilder)

    def variantAnnotationSetsGenerator(self, request, responseBuilder):
        """
        Returns a generator over the (variantAnnotationSet, nextPageToken)
        pairs defined by the specified request.
        """
        # Ensure that the variantset exists.
        self._registry_db.get_variant_set(request.variant_set_id)
        query = self._registry_db.get_variant_annotation_sets_search_query(
            request)
        self._run_sql_query(request, query, responseBuilder)

    def readsGenerator(self, request, response_builder):
        """
        Returns a generator over the (read, nextPageToken) pairs defined
        by the specified request
        """
        if not request.reference_id:
            raise exceptions.UnmappedReadsNotSupported()
        if len(request.read_group_ids) < 1:
            raise exceptions.BadRequestException(
                "At least one readGroupId must be specified")
        db = self._registry_db
        read_group_set = None
        read_groups = []
        for read_group_id in request.read_group_ids:
            read_group = db.get_read_group(read_group_id)
            if read_group_set is None:
                read_group_set = read_group.read_group_set
            elif read_group_set != read_group.read_group_set:
                raise exceptions.BadRequestException(
                    "All ReadGroups must come from the same ReadGroupSet")
            read_groups.append(read_group)
        reference = db.get_reference(request.reference_id)
        read_group_set.run_search(
            request, reference, read_groups, response_builder)

    def variantsGenerator(self, request, response_builder):
        """
        Returns a generator over the (variant, nextPageToken) pairs defined
        by the specified request.
        """
        variant_set = self._registry_db.get_variant_set(
                request.variant_set_id)
        call_sets = []
        # Avoid the query to find call_sets unless we need to.
        if len(request.call_set_ids) > 0:
            call_sets = variant_set.call_sets.filter(
                registry.CallSet.id.in_(request.call_set_ids)).all()
        if len(call_sets) != len(request.call_set_ids):
            # TODO get the offending call set ID.
            raise exceptions.CallSetNotInVariantSetException(
                    "TODO", variant_set.id)
        variant_set.run_search(request, call_sets, response_builder)

    def variantAnnotationsGenerator(self, request, response_builder):
        """
        Returns a generator over the (variantAnnotaitons, nextPageToken) pairs
        defined by the specified request.
        """
        annotation_set = self._registry_db.get_variant_annotation_set(
            request.variant_annotation_set_id)
        annotation_set.run_search(request, response_builder)

    def featuresGenerator(self, request, response_builder):
        """
        Returns a generator over the (features, nextPageToken) pairs
        defined by the (JSON string) request.
        """
        feature_set = self._registry_db.get_feature_set(request.feature_set_id)
        feature_set.run_search(request, response_builder)

        # compoundId = None
        # parentId = None
        # if request.feature_set_id is not None:
        #     compoundId = datamodel.FeatureSetCompoundId.parse(  # noqa
        #         request.feature_set_id)
        # if request.parent_id and request.parent_id != "":
        #     compoundParentId = datamodel.FeatureCompoundId.parse(  # noqa
        #         request.parent_id)
        #     parentId = compoundParentId.featureId
        #     # A client can optionally specify JUST the (compound) parentID,
        #     # and the server needs to derive the dataset & featureSet
        #     # from this (compound) parentID.
        #     if compoundId is None:
        #         compoundId = compoundParentId
        #     else:
        #         # check that the dataset and featureSet of the parent
        #         # compound ID is the same as that of the featureSetId
        #         mismatchCheck = (
        #             compoundParentId.dataset_id != compoundId.dataset_id or
        #             compoundParentId.feature_set_id !=
        #             compoundId.feature_set_id)
        #         if mismatchCheck:
        #             raise exceptions.ParentIncompatibleWithFeatureSet()

        # if compoundId is None:
        #     raise exceptions.FeatureSetNotSpecifiedException()

        # dataset = self._registry_db.getDataset(
        #     compoundId.dataset_id)
        # featureSet = dataset.getFeatureSet(compoundId.feature_set_id)
        # return featureSet.getFeatures(
        #     request.reference_name, request.start, request.end,
        #     request.page_token, request.page_size,
        #     request.feature_types, parentId)

    def callSetsGenerator(self, request, responseBuilder):
        """
        Returns a generator over the (callSet, nextPageToken) pairs defined
        by the specified request.
        """
        # Ensure that the variant set exists.
        self._registry_db.get_variant_set(request.variant_set_id)
        query = self._registry_db.get_call_sets_search_query(request)
        self._run_sql_query(request, query, responseBuilder)
        # compoundId = datamodel.VariantSetCompoundId.parse(
        #     request.variant_set_id)
        # dataset = self._registry_db.getDataset(compoundId.dataset_id)
        # variantSet = dataset.getVariantSet(compoundId.variant_set_id)
        # results = []
        # for obj in variantSet.getCallSets():
        #     include = True
        #     if request.name:
        #         if request.name != obj.getLocalId():
        #             include = False
        #     if request.bio_sample_id:
        #         if request.bio_sample_id != obj.getBioSampleId():
        #             include = False
        #     if include:
        #         results.append(obj)
        # return self._objectListGenerator(request, results)

    def featureSetsGenerator(self, request, responseBuilder):
        """
        Returns a generator over the (featureSet, nextPageToken) pairs
        defined by the specified request.
        """
        # Ensure that the dataset exists.
        self._registry_db.get_dataset(request.dataset_id)
        query = self._registry_db.get_feature_sets_search_query(request)
        self._run_sql_query(request, query, responseBuilder)

    ###########################################################
    #
    # Public API methods. Each of these methods implements the
    # corresponding API end point, and return data ready to be
    # written to the wire.
    #
    ###########################################################

    def runGetRequest(self, obj):
        """
        Runs a get request by converting the specified datamodel
        object into its protocol representation.
        """
        protocolElement = obj.get_protobuf()
        jsonString = protocol.toJson(protocolElement)
        return jsonString

    def runSearchRequest(
            self, requestStr, requestClass, responseClass, objectGenerator):
        """
        Runs the specified request. The request is a string containing
        a JSON representation of an instance of the specified requestClass.
        We return a string representation of an instance of the specified
        responseClass in JSON format. Objects are filled into the page list
        using the specified object generator, which must return
        (object, nextPageToken) pairs, and be able to resume iteration from
        any point using the nextPageToken attribute of the request object.
        """
        self.startProfile()
        try:
            request = protocol.fromJson(requestStr, requestClass)
        except protocol.json_format.ParseError:
            raise exceptions.InvalidJsonException(requestStr)
        # TODO How do we detect when the page size is not set?
        if not request.page_size:
            request.page_size = self._defaultPageSize
        if request.page_size < 0:
            raise exceptions.BadPageSizeException(request.page_size)
        responseBuilder = protocol.SearchResponseBuilder(
            responseClass, request.page_size, self._maxResponseLength)
        # nextPageToken = None
        # for obj, nextPageToken in objectGenerator(request):
        #     responseBuilder.addValue(obj)
        #     if responseBuilder.isFull():
        #         break
        # responseBuilder.setNextPageToken(nextPageToken)
        objectGenerator(request, responseBuilder)
        responseString = responseBuilder.getSerializedResponse()
        self.endProfile()
        return responseString

    def runListReferenceBases(self, id_, requestArgs):
        """
        Runs a listReferenceBases request for the specified ID and
        request arguments.
        """
        reference = self._registry_db.get_reference(id_)
        start = _parseIntegerArgument(requestArgs, 'start', 0)
        end = _parseIntegerArgument(requestArgs, 'end', reference.length)
        if end == 0:  # assume meant "get all"
            end = reference.length
        if 'pageToken' in requestArgs:
            pageTokenStr = requestArgs['pageToken']
            if pageTokenStr != "":
                start = _parsePageToken(pageTokenStr, 1)[0]

        chunkSize = self._maxResponseLength
        nextPageToken = None
        if start + chunkSize < end:
            end = start + chunkSize
            nextPageToken = str(start + chunkSize)
        reference.check_query_range(start, end)
        sequence = reference.run_get_bases(start, end)

        # build response
        response = protocol.ListReferenceBasesResponse()
        response.offset = start
        response.sequence = sequence
        if nextPageToken is not None:
            response.next_page_token = nextPageToken
        return protocol.toJson(response)

    # Get requests.

    def runGetCallSet(self, id_):
        """
        Returns a callset with the given id
        """
        call_set = self._registry_db.get_call_set(id_)
        return self.runGetRequest(call_set)

    def runGetVariant(self, id_):
        """
        Returns a variant with the given id
        """
        compound_id = registry.VariantCompoundId.parse(id_)
        variant_set = self._registry_db.get_variant_set(
                int(compound_id.variant_set_id))
        variant = variant_set.run_get_variant(compound_id)
        # TODO variant is a special case here, as it's returning a
        # protocol element rather than a datamodel object. We should
        # fix this for consistency.
        json_string = protocol.toJson(variant)
        return json_string

    def runGetBioSample(self, id_):
        """
        Runs a getBioSample request for the specified ID.
        """
        bio_sample = self._registry_db.get_bio_sample(id_)
        return self.runGetRequest(bio_sample)

    def runGetIndividual(self, id_):
        """
        Runs a getIndividual request for the specified ID.
        """
        individual = self._registry_db.get_individual(id_)
        return self.runGetRequest(individual)

    def runGetFeature(self, id_):
        """
        Returns JSON string of the feature object corresponding to
        the feature compoundID passed in.
        """
        compoundId = datamodel.FeatureCompoundId.parse(id_)  # noqa
        dataset = self._registry_db.getDataset(compoundId.dataset_id)
        featureSet = dataset.getFeatureSet(compoundId.feature_set_id)
        gaFeature = featureSet.getFeature(compoundId)
        jsonString = protocol.toJson(gaFeature)
        return jsonString

    def runGetReadGroupSet(self, id_):
        """
        Returns a readGroupSet with the given id_
        """
        read_group_set = self._registry_db.get_read_group_set(id_)
        return self.runGetRequest(read_group_set)

    def runGetReadGroup(self, id_):
        """
        Returns a read group with the given id_
        """
        read_group = self._registry_db.get_read_group(id_)
        return self.runGetRequest(read_group)

    def runGetReference(self, id_):
        """
        Runs a getReference request for the specified ID.
        """
        reference = self._registry_db.get_reference(id_)
        return self.runGetRequest(reference)

    def runGetReferenceSet(self, id_):
        """
        Runs a getReferenceSet request for the specified ID.
        """
        reference_set = self._registry_db.get_reference_set(id_)
        return self.runGetRequest(reference_set)

    def runGetVariantSet(self, id_):
        """
        Runs a getVariantSet request for the specified ID.
        """
        variant_set = self._registry_db.get_variant_set(id_)
        return self.runGetRequest(variant_set)

    def runGetFeatureSet(self, id_):
        """
        Runs a getFeatureSet request for the specified ID.
        """
        feature_set = self._registry_db.get_feature_set(id_)
        return self.runGetRequest(feature_set)

    def runGetDataset(self, id_):
        """
        Runs a getDataset request for the specified ID.
        """
        dataset = self._registry_db.get_dataset(id_)
        return self.runGetRequest(dataset)

    def runGetVariantAnnotationSet(self, id_):
        """
        Runs a getVariantSet request for the specified ID.
        """
        variant_annotation_set = self._registry_db.get_variant_annotation_set(
            id_)
        return self.runGetRequest(variant_annotation_set)

    # Search requests.

    def runSearchReadGroupSets(self, request):
        """
        Runs the specified SearchReadGroupSetsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchReadGroupSetsRequest,
            protocol.SearchReadGroupSetsResponse,
            self.readGroupSetsGenerator)

    def runSearchIndividuals(self, request):
        """
        Runs the specified search SearchIndividualsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchIndividualsRequest,
            protocol.SearchIndividualsResponse,
            self.individualsGenerator)

    def runSearchBioSamples(self, request):
        """
        Runs the specified SearchBioSamplesRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchBioSamplesRequest,
            protocol.SearchBioSamplesResponse,
            self.bioSamplesGenerator)

    def runSearchReads(self, request):
        """
        Runs the specified SearchReadsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchReadsRequest,
            protocol.SearchReadsResponse,
            self.readsGenerator)

    def runSearchReferenceSets(self, request):
        """
        Runs the specified SearchReferenceSetsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchReferenceSetsRequest,
            protocol.SearchReferenceSetsResponse,
            self.referenceSetsGenerator)

    def runSearchReferences(self, request):
        """
        Runs the specified SearchReferenceRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchReferencesRequest,
            protocol.SearchReferencesResponse,
            self.referencesGenerator)

    def runSearchVariantSets(self, request):
        """
        Runs the specified SearchVariantSetsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantSetsRequest,
            protocol.SearchVariantSetsResponse,
            self.variantSetsGenerator)

    def runSearchVariantAnnotationSets(self, request):
        """
        Runs the specified SearchVariantAnnotationSetsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantAnnotationSetsRequest,
            protocol.SearchVariantAnnotationSetsResponse,
            self.variantAnnotationSetsGenerator)

    def runSearchVariants(self, request):
        """
        Runs the specified SearchVariantRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantsRequest,
            protocol.SearchVariantsResponse,
            self.variantsGenerator)

    def runSearchVariantAnnotations(self, request):
        """
        Runs the specified SearchVariantAnnotationsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantAnnotationsRequest,
            protocol.SearchVariantAnnotationsResponse,
            self.variantAnnotationsGenerator)

    def runSearchCallSets(self, request):
        """
        Runs the specified SearchCallSetsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchCallSetsRequest,
            protocol.SearchCallSetsResponse,
            self.callSetsGenerator)

    def runSearchDatasets(self, request):
        """
        Runs the specified SearchDatasetsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchDatasetsRequest,
            protocol.SearchDatasetsResponse,
            self.datasetsGenerator)

    def runSearchFeatureSets(self, request):
        """
        Returns a SearchFeatureSetsResponse for the specified
        SearchFeatureSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchFeatureSetsRequest,
            protocol.SearchFeatureSetsResponse,
            self.featureSetsGenerator)

    def runSearchFeatures(self, request):
        """
        Returns a SearchFeaturesResponse for the specified
        SearchFeaturesRequest object.

        :param request: JSON string representing searchFeaturesRequest
        :return: JSON string representing searchFeatureResponse
        """
        return self.runSearchRequest(
            request, protocol.SearchFeaturesRequest,
            protocol.SearchFeaturesResponse,
            self.featuresGenerator)
