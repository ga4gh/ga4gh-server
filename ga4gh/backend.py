"""
Module responsible for handling protocol requests and returning
responses.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.datamodel as datamodel
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol


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
    def __init__(self, dataRepository):
        self._requestValidation = False
        self._responseValidation = False
        self._defaultPageSize = 100
        self._maxResponseLength = 2**20  # 1 MiB
        self._dataRepository = dataRepository

    def getDataRepository(self):
        """
        Get the data repository used by this backend
        """
        return self._dataRepository

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
            yield object_.toProtocolElement(), nextPageToken

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

    def datasetsGenerator(self, request):
        """
        Returns a generator over the (dataset, nextPageToken) pairs
        defined by the specified request
        """
        return self._topLevelObjectGenerator(
            request, self.getDataRepository().getNumDatasets(),
            self.getDataRepository().getDatasetByIndex)

    def readGroupSetsGenerator(self, request):
        """
        Returns a generator over the (readGroupSet, nextPageToken) pairs
        defined by the specified request.
        """
        dataset = self.getDataRepository().getDataset(request.dataset_id)
        if request.name == "":
            return self._topLevelObjectGenerator(
                request, dataset.getNumReadGroupSets(),
                dataset.getReadGroupSetByIndex)
        else:
            try:
                readGroupSet = dataset.getReadGroupSetByName(request.name)
            except exceptions.ReadGroupSetNameNotFoundException:
                return self._noObjectGenerator()
            return self._singleObjectGenerator(readGroupSet)

    def referenceSetsGenerator(self, request):
        """
        Returns a generator over the (referenceSet, nextPageToken) pairs
        defined by the specified request.
        """
        results = []
        for obj in self.getDataRepository().getReferenceSets():
            include = True
            if request.md5checksum:
                if request.md5checksum != obj.getMd5Checksum():
                    include = False
            if request.accession:
                if request.accession not in obj.getSourceAccessions():
                    include = False
            if request.assembly_id:
                if request.assembly_id != obj.getAssemblyId():
                    include = False
            if include:
                results.append(obj)
        return self._objectListGenerator(request, results)

    def referencesGenerator(self, request):
        """
        Returns a generator over the (reference, nextPageToken) pairs
        defined by the specified request.
        """
        referenceSet = self.getDataRepository().getReferenceSet(
            request.reference_set_id)
        results = []
        for obj in referenceSet.getReferences():
            include = True
            if request.md5checksum:
                if request.md5checksum != obj.getMd5Checksum():
                    include = False
            if request.accession:
                if request.accession not in obj.getSourceAccessions():
                    include = False
            if include:
                results.append(obj)
        return self._objectListGenerator(request, results)

    def variantSetsGenerator(self, request):
        """
        Returns a generator over the (variantSet, nextPageToken) pairs defined
        by the specified request.
        """
        dataset = self.getDataRepository().getDataset(request.dataset_id)
        return self._topLevelObjectGenerator(
            request, dataset.getNumVariantSets(),
            dataset.getVariantSetByIndex)

    def variantAnnotationSetsGenerator(self, request):
        """
        Returns a generator over the (variantAnnotationSet, nextPageToken)
        pairs defined by the specified request.
        """
        compoundId = datamodel.VariantSetCompoundId.parse(
            request.variant_set_id)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        variantSet = dataset.getVariantSet(request.variant_set_id)
        return self._topLevelObjectGenerator(
            request, variantSet.getNumVariantAnnotationSets(),
            variantSet.getVariantAnnotationSetByIndex)

    def readsGenerator(self, request):
        """
        Returns a generator over the (read, nextPageToken) pairs defined
        by the specified request
        """
        if not request.reference_id:
            raise exceptions.UnmappedReadsNotSupported()
        if len(request.read_group_ids) < 1:
            raise exceptions.BadRequestException(
                "At least one readGroupId must be specified")
        elif len(request.read_group_ids) == 1:
            return self._readsGeneratorSingle(request)
        else:
            return self._readsGeneratorMultiple(request)

    def _readsGeneratorSingle(self, request):
        compoundId = datamodel.ReadGroupCompoundId.parse(
            request.read_group_ids[0])
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        readGroupSet = dataset.getReadGroupSet(compoundId.read_group_set_id)
        referenceSet = readGroupSet.getReferenceSet()
        if referenceSet is None:
            raise exceptions.ReadGroupSetNotMappedToReferenceSetException(
                    readGroupSet.getId())
        reference = referenceSet.getReference(request.reference_id)
        readGroup = readGroupSet.getReadGroup(compoundId.read_group_id)
        intervalIterator = ReadsIntervalIterator(
            request, readGroup, reference)
        return intervalIterator

    def _readsGeneratorMultiple(self, request):
        compoundId = datamodel.ReadGroupCompoundId.parse(
            request.read_group_ids[0])
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        readGroupSet = dataset.getReadGroupSet(compoundId.read_group_set_id)
        referenceSet = readGroupSet.getReferenceSet()
        if referenceSet is None:
            raise exceptions.ReadGroupSetNotMappedToReferenceSetException(
                    readGroupSet.getId())
        reference = referenceSet.getReference(request.reference_id)
        readGroupIds = readGroupSet.getReadGroupIds()
        if set(readGroupIds) != set(request.read_group_ids):
            raise exceptions.BadRequestException(
                "If multiple readGroupIds are specified, "
                "they must be all of the readGroupIds in a ReadGroup")
        intervalIterator = ReadsIntervalIterator(
            request, readGroupSet, reference)
        return intervalIterator

    def variantsGenerator(self, request):
        """
        Returns a generator over the (variant, nextPageToken) pairs defined
        by the specified request.
        """
        compoundId = datamodel.VariantSetCompoundId \
            .parse(request.variant_set_id)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        variantSet = dataset.getVariantSet(compoundId.variant_set_id)
        intervalIterator = VariantsIntervalIterator(request, variantSet)
        return intervalIterator

    def variantAnnotationsGenerator(self, request):
        """
        Returns a generator over the (variantAnnotaitons, nextPageToken) pairs
        defined by the specified request.
        """
        compoundId = datamodel.VariantAnnotationSetCompoundId.parse(
            request.variant_annotation_set_id)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        variantSet = dataset.getVariantSet(compoundId.variant_set_id)
        variantAnnotationSet = variantSet.getVariantAnnotationSet(
            request.variant_annotation_set_id)
        intervalIterator = VariantAnnotationsIntervalIterator(
            request, variantAnnotationSet)
        return intervalIterator

    def featuresGenerator(self, request):
        """
        Returns a generator over the (features, nextPageToken) pairs
        defined by the (JSON string) request.
        """
        compoundId = None
        parentId = None
        if request.feature_set_id is not None:
            compoundId = datamodel.FeatureSetCompoundId.parse(
                request.feature_set_id)
        if request.parent_id and request.parent_id != "":
            compoundParentId = datamodel.FeatureCompoundId.parse(
                request.parent_id)
            parentId = compoundParentId.featureId
            # A client can optionally specify JUST the (compound) parentID,
            # and the server needs to derive the dataset & featureSet
            # from this (compound) parentID.
            if compoundId is None:
                compoundId = compoundParentId
            else:
                # check that the dataset and featureSet of the parent
                # compound ID is the same as that of the featureSetId
                mismatchCheck = (
                    compoundParentId.dataset_id != compoundId.dataset_id or
                    compoundParentId.feature_set_id !=
                    compoundId.feature_set_id)
                if mismatchCheck:
                    raise exceptions.ParentIncompatibleWithFeatureSet()

        if compoundId is None:
            raise exceptions.FeatureSetNotSpecifiedException()

        dataset = self.getDataRepository().getDataset(
            compoundId.dataset_id)
        featureSet = dataset.getFeatureSet(compoundId.feature_set_id)
        return featureSet.getFeatures(
            request.reference_name, request.start, request.end,
            request.page_token, request.page_size,
            request.feature_types, parentId)

    def callSetsGenerator(self, request):
        """
        Returns a generator over the (callSet, nextPageToken) pairs defined
        by the specified request.
        """
        compoundId = datamodel.VariantSetCompoundId.parse(
            request.variant_set_id)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        variantSet = dataset.getVariantSet(compoundId.variant_set_id)
        if request.name == "":
            return self._topLevelObjectGenerator(
                request, variantSet.getNumCallSets(),
                variantSet.getCallSetByIndex)
        else:
            try:
                callSet = variantSet.getCallSetByName(request.name)
            except exceptions.CallSetNameNotFoundException:
                return self._noObjectGenerator()
            return self._singleObjectGenerator(callSet)

    def featureSetsGenerator(self, request):
        """
        Returns a generator over the (featureSet, nextPageToken) pairs
        defined by the specified request.
        """
        dataset = self.getDataRepository().getDataset(request.dataset_id)
        return self._topLevelObjectGenerator(
            request, dataset.getNumFeatureSets(),
            dataset.getFeatureSetByIndex)

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
        protocolElement = obj.toProtocolElement()
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
        nextPageToken = None
        for obj, nextPageToken in objectGenerator(request):
            responseBuilder.addValue(obj)
            if responseBuilder.isFull():
                break
        responseBuilder.setNextPageToken(nextPageToken)
        responseString = responseBuilder.getSerializedResponse()
        self.endProfile()
        return responseString

    def runListReferenceBases(self, id_, requestArgs):
        """
        Runs a listReferenceBases request for the specified ID and
        request arguments.
        """
        compoundId = datamodel.ReferenceCompoundId.parse(id_)
        referenceSet = self.getDataRepository().getReferenceSet(
            compoundId.reference_set_id)
        reference = referenceSet.getReference(id_)
        start = _parseIntegerArgument(requestArgs, 'start', 0)
        end = _parseIntegerArgument(requestArgs, 'end', reference.getLength())
        if end == 0:  # assume meant "get all"
            end = reference.getLength()
        if 'pageToken' in requestArgs:
            pageTokenStr = requestArgs['pageToken']
            if pageTokenStr != "":
                start = _parsePageToken(pageTokenStr, 1)[0]

        chunkSize = self._maxResponseLength
        nextPageToken = None
        if start + chunkSize < end:
            end = start + chunkSize
            nextPageToken = str(start + chunkSize)
        sequence = reference.getBases(start, end)

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
        compoundId = datamodel.CallSetCompoundId.parse(id_)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        variantSet = dataset.getVariantSet(compoundId.variant_set_id)
        callSet = variantSet.getCallSet(id_)
        return self.runGetRequest(callSet)

    def runGetVariant(self, id_):
        """
        Returns a variant with the given id
        """
        compoundId = datamodel.VariantCompoundId.parse(id_)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        variantSet = dataset.getVariantSet(compoundId.variant_set_id)
        gaVariant = variantSet.getVariant(compoundId)
        # TODO variant is a special case here, as it's returning a
        # protocol element rather than a datamodel object. We should
        # fix this for consistency.
        jsonString = protocol.toJson(gaVariant)
        return jsonString

    def runGetFeature(self, id_):
        """
        Returns JSON string of the feature object corresponding to
        the feature compoundID passed in.
        """
        compoundId = datamodel.FeatureCompoundId.parse(id_)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        featureSet = dataset.getFeatureSet(compoundId.feature_set_id)
        gaFeature = featureSet.getFeature(compoundId)
        jsonString = protocol.toJson(gaFeature)
        return jsonString

    def runGetReadGroupSet(self, id_):
        """
        Returns a readGroupSet with the given id_
        """
        compoundId = datamodel.ReadGroupSetCompoundId.parse(id_)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        readGroupSet = dataset.getReadGroupSet(id_)
        return self.runGetRequest(readGroupSet)

    def runGetReadGroup(self, id_):
        """
        Returns a read group with the given id_
        """
        compoundId = datamodel.ReadGroupCompoundId.parse(id_)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        readGroupSet = dataset.getReadGroupSet(compoundId.read_group_set_id)
        readGroup = readGroupSet.getReadGroup(id_)
        return self.runGetRequest(readGroup)

    def runGetReference(self, id_):
        """
        Runs a getReference request for the specified ID.
        """
        compoundId = datamodel.ReferenceCompoundId.parse(id_)
        referenceSet = self.getDataRepository().getReferenceSet(
            compoundId.reference_set_id)
        reference = referenceSet.getReference(id_)
        return self.runGetRequest(reference)

    def runGetReferenceSet(self, id_):
        """
        Runs a getReferenceSet request for the specified ID.
        """
        referenceSet = self.getDataRepository().getReferenceSet(id_)
        return self.runGetRequest(referenceSet)

    def runGetVariantSet(self, id_):
        """
        Runs a getVariantSet request for the specified ID.
        """
        compoundId = datamodel.VariantSetCompoundId.parse(id_)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        variantSet = dataset.getVariantSet(id_)
        return self.runGetRequest(variantSet)

    def runGetFeatureSet(self, id_):
        """
        Runs a getFeatureSet request for the specified ID.
        """
        compoundId = datamodel.FeatureSetCompoundId.parse(id_)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        featureSet = dataset.getFeatureSet(id_)
        return self.runGetRequest(featureSet)

    def runGetDataset(self, id_):
        """
        Runs a getDataset request for the specified ID.
        """
        dataset = self.getDataRepository().getDataset(id_)
        return self.runGetRequest(dataset)

    def runGetVariantAnnotationSet(self, id_):
        """
        Runs a getVariantSet request for the specified ID.
        """
        compoundId = datamodel.VariantAnnotationSetCompoundId.parse(id_)
        dataset = self.getDataRepository().getDataset(compoundId.dataset_id)
        variantSet = dataset.getVariantSet(compoundId.variant_set_id)
        variantAnnotationSet = variantSet.getVariantAnnotationSet(id_)
        return self.runGetRequest(variantAnnotationSet)

    # Search requests.

    def runSearchReadGroupSets(self, request):
        """
        Runs the specified SearchReadGroupSetsRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchReadGroupSetsRequest,
            protocol.SearchReadGroupSetsResponse,
            self.readGroupSetsGenerator)

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
