"""
Module responsible for handling protocol requests and returning
responses.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import json
import random

import ga4gh.protocol as protocol
import ga4gh.datamodel.references as references
import ga4gh.exceptions as exceptions
import ga4gh.datamodel as datamodel
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.variants as variants


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


def _getVariantSet(request, variantSetIdMap):
    variantSetId = request.variantSetId
    variantSet = _safeMapQuery(
        variantSetIdMap, variantSetId,
        exceptionClass=exceptions.VariantSetNotFoundException)
    return variantSet


def _safeMapQuery(idMap, id_, exceptionClass=None, idErrorString=None):
    """
    Attempt to retrieve a value from a map, throw an appropriate error
    if the key is not present
    """
    try:
        obj = idMap[id_]
    except KeyError:
        if idErrorString is None:
            idErrorString = id_
        if exceptionClass is None:
            exceptionClass = exceptions.ObjectWithIdNotFoundException
        raise exceptionClass(idErrorString)
    return obj


class IntervalIterator(object):
    """
    Implements generator logic for types which accept a start/end
    range to search for the object. Returns an iterator over
    (object, pageToken) pairs. The pageToken is a string which allows
    us to pick up the iteration at any point, and is None for the last
    value in the iterator.
    """
    def __init__(self, request, containerIdMap):
        self._request = request
        self._containerIdMap = containerIdMap
        self._container = self._getContainer()
        self._searchIterator = None
        self._currentObject = None
        self._nextObject = None
        self._searchAnchor = None
        self._distanceFromAnchor = None
        if request.pageToken is None:
            self._initialiseIteration()
        else:
            # Set the search start point and the number of records to skip from
            # the page token.
            searchAnchor, objectsToSkip = _parsePageToken(request.pageToken, 2)
            self._pickUpIteration(searchAnchor, objectsToSkip)

    def _initialiseIteration(self):
        """
        Starts a new iteration.
        """
        self._searchIterator = self._search(
            self._request.start, self._request.end)
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
        self._searchIterator = self._search(searchAnchor, self._request.end)
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
                assert self._getStart(obj) == searchAnchor
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
        ret = self._currentObject, nextPageToken
        self._currentObject = self._nextObject
        self._nextObject = next(self._searchIterator, None)
        return ret

    def __iter__(self):
        return self


class ReadsIntervalIterator(IntervalIterator):
    """
    An interval iterator for reads
    """
    def _getContainer(self):
        if len(self._request.readGroupIds) != 1:
            if len(self._request.readGroupIds) == 0:
                msg = "Read search requires a readGroup to be specified"
            else:
                msg = "Read search over multiple readGroups not supported"
            raise exceptions.NotImplementedException(msg)
        readGroupId = self._request.readGroupIds[0]
        readGroup = _safeMapQuery(
            self._containerIdMap, readGroupId,
            exceptions.ReadGroupNotFoundException)
        return readGroup

    def _search(self, start, end):
        return self._container.getReadAlignments(
            self._request.referenceId, start, end)

    @classmethod
    def _getStart(cls, readAlignment):
        return readAlignment.alignment.position.position

    @classmethod
    def _getEnd(cls, readAlignment):
        return cls._getStart(readAlignment) + \
            len(readAlignment.alignedSequence)


class VariantsIntervalIterator(IntervalIterator):
    """
    An interval iterator for variants
    """
    def _getContainer(self):
        return _getVariantSet(self._request, self._containerIdMap)

    def _search(self, start, end):
        return self._container.getVariants(
            self._request.referenceName, start, end, self._request.variantName,
            self._request.callSetIds)

    @classmethod
    def _getStart(cls, variant):
        return variant.start

    @classmethod
    def _getEnd(cls, variant):
        return variant.end


class AbstractBackend(object):
    """
    An abstract GA4GH backend.
    This class provides methods for all of the GA4GH protocol end points.
    """
    def __init__(self):
        self._referenceSetIdMap = {}
        self._referenceSetIds = []
        self._referenceIdMap = {}
        self._referenceIds = []
        self._requestValidation = False
        self._responseValidation = False
        self._defaultPageSize = 100
        self._maxResponseLength = 2**20  # 1 MiB
        self._datasetIdMap = {}
        self._datasetIds = []

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

    def getDatasetIds(self):
        """
        Returns a list of datasets in this backend
        """
        return self._datasetIds

    def getDataset(self, datasetId):
        """
        Returns a dataset with id datasetId
        """
        return _safeMapQuery(
            self._datasetIdMap, datasetId,
            exceptions.DatasetNotFoundException)

    def getReferenceSets(self):
        """
        Returns the list of ReferenceSets in this backend
        """
        return list(self._referenceSetIdMap.values())

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
            if not requestClass.validate(jsonDict):
                raise exceptions.RequestValidationFailureException(
                    jsonDict, requestClass)

    def validateResponse(self, jsonString, responseClass):
        """
        Ensures the jsonDict corresponds to a valid instance of responseClass
        Throws an error if the data is invalid
        """
        if self._responseValidation:
            jsonDict = json.loads(jsonString)
            if not responseClass.validate(jsonDict):
                raise exceptions.ResponseValidationFailureException(
                    jsonDict, responseClass)

    ###########################################################
    #
    # Iterators over the data hierarchy. These methods help to
    # implement the search endpoints by providing iterators
    # over the objects to be returned to the client.
    #
    ###########################################################

    def _topLevelObjectGenerator(self, request, idMap, idList):
        """
        Generalisation of the code to iterate over the objects at the top
        of the data hierarchy.
        """
        currentIndex = 0
        if request.pageToken is not None:
            currentIndex, = _parsePageToken(request.pageToken, 1)
        while currentIndex < len(idList):
            objectId = idList[currentIndex]
            object_ = idMap[objectId]
            currentIndex += 1
            nextPageToken = None
            if currentIndex < len(idList):
                nextPageToken = str(currentIndex)
            yield object_.toProtocolElement(), nextPageToken

    def datasetsGenerator(self, request):
        """
        Returns a generator over the (dataset, nextPageToken) pairs
        defined by the specified request
        """
        return self._topLevelObjectGenerator(
            request, self._datasetIdMap, self._datasetIds)

    def readGroupSetsGenerator(self, request):
        """
        Returns a generator over the (readGroupSet, nextPageToken) pairs
        defined by the specified request.
        """
        dataset = self.getDataset(request.datasetId)
        return self._topLevelObjectGenerator(
            request, dataset.getReadGroupSetIdMap(),
            dataset.getReadGroupSetIds())

    def referenceSetsGenerator(self, request):
        """
        Returns a generator over the (referenceSet, nextPageToken) pairs
        defined by the specified request.
        """
        return self._topLevelObjectGenerator(
            request, self._referenceSetIdMap, self._referenceSetIds)

    def referencesGenerator(self, request):
        """
        Returns a generator over the (reference, nextPageToken) pairs
        defined by the specified request.
        """
        return self._topLevelObjectGenerator(
            request, self._referenceIdMap, self._referenceIds)

    def variantSetsGenerator(self, request):
        """
        Returns a generator over the (variantSet, nextPageToken) pairs defined
        by the specified request.
        """
        dataset = self.getDataset(request.datasetId)
        return self._topLevelObjectGenerator(
            request, dataset.getVariantSetIdMap(),
            dataset.getVariantSetIds())

    def readsGenerator(self, request):
        """
        Returns a generator over the (read, nextPageToken) pairs defined
        by the specified request
        """
        if request.referenceId is None:
            raise exceptions.UnmappedReadsNotSupported()
        if len(request.readGroupIds) != 1:
            raise exceptions.NotImplementedException(
                "Exactly one read group id must be specified")
        compoundId = reads.CompoundReadGroupId(request.readGroupIds[0])
        dataset = self.getDataset(compoundId.datasetId)
        intervalIterator = ReadsIntervalIterator(
            request, dataset.getReadGroupIdMap())
        return intervalIterator

    def variantsGenerator(self, request):
        """
        Returns a generator over the (variant, nextPageToken) pairs defined
        by the specified request.
        """
        compoundId = variants.CompoundVariantSetId(request.variantSetId)
        dataset = self.getDataset(compoundId.datasetId)
        intervalIterator = VariantsIntervalIterator(
            request, dataset.getVariantSetIdMap())
        return intervalIterator

    def callSetsGenerator(self, request):
        """
        Returns a generator over the (callSet, nextPageToken) pairs defined
        by the specified request.
        """
        if request.name is not None:
            raise exceptions.NotImplementedException(
                "Searching over names is not supported")
        compoundId = variants.CompoundVariantSetId(request.variantSetId)
        dataset = self.getDataset(compoundId.datasetId)
        variantSet = _getVariantSet(
            request, dataset.getVariantSetIdMap())
        return self._topLevelObjectGenerator(
            request, variantSet.getCallSetIdMap(),
            variantSet.getCallSetIds())

    ###########################################################
    #
    # Public API methods. Each of these methods implements the
    # corresponding API end point, and return data ready to be
    # written to the wire.
    #
    ###########################################################

    def runGetRequest(self, idMap, id_):
        """
        Runs a get request by indexing into the provided idMap and
        returning a json string of that object
        """
        obj = _safeMapQuery(idMap, id_)
        protocolElement = obj.toProtocolElement()
        jsonString = protocolElement.toJsonString()
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
            requestDict = json.loads(requestStr)
        except ValueError:
            raise exceptions.InvalidJsonException(requestStr)
        self.validateRequest(requestDict, requestClass)
        request = requestClass.fromJsonDict(requestDict)
        if request.pageSize is None:
            request.pageSize = self._defaultPageSize
        if request.pageSize <= 0:
            raise exceptions.BadPageSizeException(request.pageSize)
        responseBuilder = protocol.SearchResponseBuilder(
            responseClass, request.pageSize, self._maxResponseLength)
        nextPageToken = None
        for obj, nextPageToken in objectGenerator(request):
            responseBuilder.addValue(obj)
            if responseBuilder.isFull():
                break
        responseBuilder.setNextPageToken(nextPageToken)
        responseString = responseBuilder.getJsonString()
        self.validateResponse(responseString, responseClass)
        self.endProfile()
        return responseString

    def runGetCallset(self, id_):
        """
        Returns a callset with the given id
        """
        compoundId = variants.CompoundCallsetId(id_)
        dataset = self.getDataset(compoundId.datasetId)
        variantSet = _getVariantSet(
            compoundId, dataset.getVariantSetIdMap())
        return self.runGetRequest(variantSet.getCallSetIdMap(), id_)

    def runListReferenceBases(self, id_, requestArgs):
        """
        Runs a listReferenceBases request for the specified ID and
        request arguments.
        """
        # parse arguments
        reference = _safeMapQuery(
            self._referenceIdMap, id_,
            exceptions.ObjectWithIdNotFoundException)
        start = 0
        end = datamodel.PysamDatamodelMixin.fastaMax
        if 'start' in requestArgs:
            startString = requestArgs['start']
            try:
                start = int(startString)
            except ValueError:
                raise exceptions.BadRequestIntegerException(
                    'start', startString)
        if 'end' in requestArgs:
            endString = requestArgs['end']
            try:
                end = int(endString)
            except ValueError:
                raise exceptions.BadRequestIntegerException(
                    'end', endString)
        if 'pageToken' in requestArgs:
            pageTokenStr = requestArgs['pageToken']
            start = _parsePageToken(pageTokenStr, 1)[0]
        chunkSize = self._maxResponseLength

        # get reference bases
        gbEnd = min(start + chunkSize, end)
        sequence = reference.getBases(start, gbEnd)

        # determine nextPageToken
        if len(sequence) == chunkSize:
            nextPageToken = start + chunkSize
        elif len(sequence) > chunkSize:
            raise exceptions.ServerError()  # should never happen
        else:
            nextPageToken = None

        # build response
        response = protocol.ListReferenceBasesResponse()
        response.offset = start
        response.sequence = sequence
        response.nextPageToken = nextPageToken
        return response.toJsonString()

    # Get requests.

    def runGetReadGroupSet(self, id_):
        """
        Returns a readGroupSet with the given id_
        """
        compoundId = reads.CompoundReadGroupSetId(id_)
        dataset = self.getDataset(compoundId.datasetId)
        return self.runGetRequest(dataset.getReadGroupSetIdMap(), id_)

    def runGetReadGroup(self, id_):
        """
        Returns a read group with the given id_
        """
        compoundId = reads.CompoundReadGroupId(id_)
        dataset = self.getDataset(compoundId.datasetId)
        return self.runGetRequest(dataset.getReadGroupIdMap(), id_)

    def runGetReference(self, id_):
        """
        Runs a getReference request for the specified ID.
        """
        return self.runGetRequest(self._referenceIdMap, id_)

    def runGetReferenceSet(self, id_):
        """
        Runs a getReferenceSet request for the specified ID.
        """
        return self.runGetRequest(self._referenceSetIdMap, id_)

    def runGetVariantSet(self, id_):
        """
        Runs a getVariantSet request for the specified ID.
        """
        compoundId = variants.CompoundVariantSetId(id_)
        dataset = self.getDataset(compoundId.datasetId)
        return self.runGetRequest(
            dataset.getVariantSetIdMap(), str(compoundId))

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

    def runSearchVariants(self, request):
        """
        Runs the specified SearchVariantRequest.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantsRequest,
            protocol.SearchVariantsResponse,
            self.variantsGenerator)

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


class EmptyBackend(AbstractBackend):
    """
    A GA4GH backend that contains no data.
    """


class SimulatedBackend(AbstractBackend):
    """
    A GA4GH backend backed by no data; used mostly for testing
    """
    def __init__(self, randomSeed=0, numCalls=1, variantDensity=0.5,
                 numVariantSets=1, numReferenceSets=1,
                 numReferencesPerReferenceSet=1, numAlignments=2):
        super(SimulatedBackend, self).__init__()

        # Datasets
        dataset1 = datasets.SimulatedDataset(
            "simulatedDataset1", randomSeed, numCalls,
            variantDensity, numVariantSets, numAlignments)
        dataset2 = datasets.SimulatedDataset(
            "simulatedDataset2", randomSeed, numCalls,
            variantDensity, numVariantSets, numAlignments)
        self._datasetIdMap[dataset1.getId()] = dataset1
        self._datasetIdMap[dataset2.getId()] = dataset2
        self._datasetIds = sorted(self._datasetIdMap.keys())

        # References
        randomGenerator = random.Random()
        randomGenerator.seed(randomSeed)
        for i in range(numReferenceSets):
            referenceSetId = "referenceSet{}".format(i)
            referenceSetSeed = randomGenerator.getrandbits(32)
            referenceSet = references.SimulatedReferenceSet(
                referenceSetId, referenceSetSeed,
                numReferencesPerReferenceSet)
            self._referenceSetIdMap[referenceSetId] = referenceSet
            for reference in referenceSet.getReferences():
                referenceId = reference.getId()
                self._referenceIdMap[referenceId] = reference
        self._referenceSetIds = sorted(self._referenceSetIdMap.keys())
        self._referenceIds = sorted(self._referenceIdMap.keys())


class FileSystemBackend(AbstractBackend):
    """
    A GA4GH backend backed by data on the file system
    """
    def __init__(self, dataDir):
        super(FileSystemBackend, self).__init__()
        self._dataDir = dataDir
        # TODO this code is very ugly and should be regarded as a temporary
        # stop-gap until we deal with iterating over the data tree properly.

        # References
        referencesDirName = "references"
        referenceSetDir = os.path.join(self._dataDir, referencesDirName)
        for referenceSetId in os.listdir(referenceSetDir):
            relativePath = os.path.join(referenceSetDir, referenceSetId)
            if os.path.isdir(relativePath):
                referenceSet = references.HtslibReferenceSet(
                    referenceSetId, relativePath)
                self._referenceSetIdMap[referenceSetId] = referenceSet
                for reference in referenceSet.getReferences():
                    referenceId = reference.getId()
                    self._referenceIdMap[referenceId] = reference
        self._referenceSetIds = sorted(self._referenceSetIdMap.keys())
        self._referenceIds = sorted(self._referenceIdMap.keys())

        # Datasets
        datasetDirs = [
            os.path.join(self._dataDir, directory)
            for directory in os.listdir(self._dataDir)
            if os.path.isdir(os.path.join(self._dataDir, directory)) and
            directory != referencesDirName]
        for datasetDir in datasetDirs:
            dataset = datasets.FileSystemDataset(datasetDir)
            self._datasetIdMap[dataset.getId()] = dataset
        self._datasetIds = sorted(self._datasetIdMap.keys())
