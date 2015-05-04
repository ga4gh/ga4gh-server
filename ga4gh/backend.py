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
import ga4gh.datamodel.reads as reads
import ga4gh.exceptions as exceptions
import ga4gh.datamodel.variants as variants
import ga4gh.datamodel as datamodel


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
    if len(request.variantSetIds) != 1:
        if len(request.variantSetIds) == 0:
            msg = "Variant search requires specifying a variantSet"
        else:
            msg = ("Variant search over multiple variantSets "
                   "not supported")
        raise exceptions.NotImplementedException(msg)
    variantSetId = request.variantSetIds[0]
    try:
        variantSet = variantSetIdMap[variantSetId]
    except KeyError:
        raise exceptions.VariantSetNotFoundException(variantSetId)
    return variantSet


class IntervalIterator(object):
    """
    Implements generator logic for types which accept a start/end
    range to search for the object
    """
    def __init__(self, request, containerIdMap):
        self._request = request
        self._containerIdMap = containerIdMap
        self._badPageTokenExceptionMessage = (
            "Inconsistent page token provided")
        self._container = self._getContainer()
        self._startPosition, self._equalPositionsToSkip = \
            self._getIntervalCounters()
        self._iterator = self._getIterator()
        self._generator = self._internalIterator()

    def __iter__(self):
        return self._generator

    def next(self):
        obj = next(self._generator)
        return obj

    def _raiseBadPageTokenException(self):
        raise exceptions.BadPageTokenException(
            self._badPageTokenExceptionMessage)

    def _getIntervalCounters(self):
        startPosition = self._request.start
        equalPositionsToSkip = 0
        if self._request.pageToken is not None:
            startPosition, equalPositionsToSkip = _parsePageToken(
                self._request.pageToken, 2)
        return startPosition, equalPositionsToSkip

    def _internalIterator(self):
        obj = next(self._iterator, None)
        if self._request.pageToken is not None:
            # First, skip any records with getStart < startPosition
            # or getEnd < request.start
            while (self._getStart(obj) < self._startPosition or
                   self._getEnd(obj) < self._request.start):
                obj = next(self._iterator, None)
                if obj is None:
                    self._raiseBadPageTokenException()
            # Now, skip equalPositionsToSkip records which have getStart
            # == startPosition
            equalPositionsSkipped = 0
            while equalPositionsSkipped < self._equalPositionsToSkip:
                if self._getStart(obj) != self._startPosition:
                    self._raiseBadPageTokenException()
                equalPositionsSkipped += 1
                obj = next(self._iterator, None)
                if obj is None:
                    self._raiseBadPageTokenException()
        # iterator is now positioned to start yielding valid records
        while obj is not None:
            nextObj = next(self._iterator, None)
            nextPageToken = None
            if nextObj is not None:
                if self._getStart(obj) == self._getStart(nextObj):
                    self._equalPositionsToSkip += 1
                else:
                    self._equalPositionsToSkip = 0
                nextPageToken = "{}:{}".format(
                    self._getStart(nextObj), self._equalPositionsToSkip)
            yield obj, nextPageToken
            obj = nextObj


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
        try:
            readGroup = self._containerIdMap[self._request.readGroupIds[0]]
        except KeyError:
            raise exceptions.ReadGroupNotFoundException(readGroupId)
        return readGroup

    def _getIterator(self):
        iterator = self._container.getReadAlignments(
            self._request.referenceId,
            self._startPosition, self._request.end)
        return iterator

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

    def _getIterator(self):
        iterator = self._container.getVariants(
            self._request.referenceName, self._startPosition,
            self._request.end, self._request.variantName,
            self._request.callSetIds)
        return iterator

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
        self._variantSetIdMap = {}
        self._variantSetIds = []
        self._referenceSetIdMap = {}
        self._referenceSetIds = []
        self._referenceIdMap = {}
        self._referenceIds = []
        self._readGroupSetIdMap = {}
        self._readGroupSetIds = []
        self._readGroupIds = []
        self._readGroupIdMap = {}
        self._requestValidation = False
        self._responseValidation = False
        self._defaultPageSize = 100
        self._maxResponseLength = 2**20  # 1 MiB

    def getVariantSets(self):
        """
        Returns the list of VariantSets in this backend.
        """
        return list(self._variantSetIdMap.values())

    def getReadGroupSets(self):
        """
        Returns the list of ReadGroupSets in this backend.
        """
        return list(self._readGroupSetIdMap.values())

    def getReferenceSets(self):
        """
        Returns the list of ReferenceSets in this backend
        """
        return list(self._referenceSetIdMap.values())

    def getReference(self, id_):
        """
        Returns a reference with the given id_
        """
        return self.runGetRequest(self._referenceIdMap, id_)

    def getReferenceSet(self, id_):
        """
        Returns a referenceSet with the given id_
        """
        return self.runGetRequest(self._referenceSetIdMap, id_)

    def listReferenceBases(self, id_, requestArgs):
        # parse arguments
        try:
            reference = self._referenceIdMap[id_]
        except KeyError:
            raise exceptions.ObjectWithIdNotFoundException(id_)
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

    def runGetRequest(self, idMap, id_):
        """
        Runs a get request by indexing into the provided idMap and
        returning a json string of that object
        """
        try:
            obj = idMap[id_]
        except KeyError:
            raise exceptions.ObjectWithIdNotFoundException(id_)
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

    def searchReadGroupSets(self, request):
        """
        Returns a GASearchReadGroupSetsResponse for the specified
        GASearchReadGroupSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReadGroupSetsRequest,
            protocol.SearchReadGroupSetsResponse,
            self.readGroupSetsGenerator)

    def searchReads(self, request):
        """
        Returns a GASearchReadsResponse for the specified
        GASearchReadsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReadsRequest,
            protocol.SearchReadsResponse,
            self.readsGenerator)

    def searchReferenceSets(self, request):
        """
        Returns a GASearchReferenceSetsResponse for the specified
        GASearchReferenceSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReferenceSetsRequest,
            protocol.SearchReferenceSetsResponse,
            self.referenceSetsGenerator)

    def searchReferences(self, request):
        """
        Returns a GASearchReferencesResponse for the specified
        GASearchReferencesRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReferencesRequest,
            protocol.SearchReferencesResponse,
            self.referencesGenerator)

    def searchVariantSets(self, request):
        """
        Returns a GASearchVariantSetsResponse for the specified
        GASearchVariantSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantSetsRequest,
            protocol.SearchVariantSetsResponse,
            self.variantSetsGenerator)

    def searchVariants(self, request):
        """
        Returns a GASearchVariantsResponse for the specified
        GASearchVariantsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantsRequest,
            protocol.SearchVariantsResponse,
            self.variantsGenerator)

    def searchCallSets(self, request):
        """
        Returns a GASearchCallSetsResponse for the specified
        GASearchCallSetsRequest Object.
        """
        return self.runSearchRequest(
            request, protocol.SearchCallSetsRequest,
            protocol.SearchCallSetsResponse,
            self.callSetsGenerator)

    # Iterators over the data hieararchy

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

    def readGroupSetsGenerator(self, request):
        """
        Returns a generator over the (readGroupSet, nextPageToken) pairs
        defined by the specified request.
        """
        return self._topLevelObjectGenerator(
            request, self._readGroupSetIdMap, self._readGroupSetIds)

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
        return self._topLevelObjectGenerator(
            request, self._variantSetIdMap, self._variantSetIds)

    def readsGenerator(self, request):
        """
        Returns a generator over the (read, nextPageToken) pairs defined
        by the specified request
        """
        intervalIterator = ReadsIntervalIterator(
            request, self._readGroupIdMap)
        return intervalIterator

    def variantsGenerator(self, request):
        """
        Returns a generator over the (variant, nextPageToken) pairs defined
        by the specified request.
        """
        intervalIterator = VariantsIntervalIterator(
            request, self._variantSetIdMap)
        return intervalIterator

    def callSetsGenerator(self, request):
        """
        Returns a generator over the (callSet, nextPageToken) pairs defined
        by the specified request.
        """
        if request.name is not None:
            raise exceptions.NotImplementedException(
                "Searching over names is not supported")
        variantSet = _getVariantSet(request, self._variantSetIdMap)
        return self._topLevelObjectGenerator(
            request, variantSet.getCallSetIdMap(),
            variantSet.getCallSetIds())

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


class EmptyBackend(AbstractBackend):
    """
    A GA4GH backend that contains no data.
    """


class SimulatedBackend(AbstractBackend):
    """
    A GA4GH backend backed by no data; used mostly for testing
    """
    def __init__(self, randomSeed=0, numCalls=1, variantDensity=0.5,
                 numVariantSets=1):
        super(SimulatedBackend, self).__init__()
        self._randomSeed = randomSeed
        self._randomGenerator = random.Random()
        self._randomGenerator.seed(self._randomSeed)
        for i in range(numVariantSets):
            variantSetId = "simVs{}".format(i)
            seed = self._randomGenerator.randint(0, 2**32 - 1)
            variantSet = variants.SimulatedVariantSet(
                seed, numCalls, variantDensity, variantSetId)
            self._variantSetIdMap[variantSetId] = variantSet
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

        # Reads
        readGroupSetId = "aReadGroupSet"
        readGroupSet = reads.SimulatedReadGroupSet(readGroupSetId)
        self._readGroupSetIdMap[readGroupSetId] = readGroupSet
        for readGroup in readGroupSet.getReadGroups():
            self._readGroupIdMap[readGroup.getId()] = readGroup
        self._readGroupSetIds = sorted(self._readGroupSetIdMap.keys())
        self._readGroupIds = sorted(self._readGroupIdMap.keys())

        # References
        referenceSetId = "aReferenceSet"
        referenceSet = references.SimulatedReferenceSet(referenceSetId)
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
        # Variants
        variantSetDir = os.path.join(self._dataDir, "variants")
        for variantSetId in os.listdir(variantSetDir):
            relativePath = os.path.join(variantSetDir, variantSetId)
            if os.path.isdir(relativePath):
                self._variantSetIdMap[variantSetId] = \
                    variants.HtslibVariantSet(variantSetId, relativePath)
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

        # References
        referenceSetDir = os.path.join(self._dataDir, "references")
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

        # Reads
        readGroupSetDir = os.path.join(self._dataDir, "reads")
        for readGroupSetId in os.listdir(readGroupSetDir):
            relativePath = os.path.join(readGroupSetDir, readGroupSetId)
            if os.path.isdir(relativePath):
                readGroupSet = reads.HtslibReadGroupSet(
                    readGroupSetId, relativePath)
                self._readGroupSetIdMap[readGroupSetId] = readGroupSet
                for readGroup in readGroupSet.getReadGroups():
                    self._readGroupIdMap[readGroup.getId()] = readGroup
        self._readGroupSetIds = sorted(self._readGroupSetIdMap.keys())
        self._readGroupIds = sorted(self._readGroupIdMap.keys())
