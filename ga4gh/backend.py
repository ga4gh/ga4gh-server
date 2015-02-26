"""
Module responsible for handling protocol requests and returning
responses.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.protocol as protocol
import ga4gh.datamodel.references as references
import ga4gh.datamodel.reads as reads


class Backend(object):
    """
    The GA4GH backend. This class provides methods for all of the GA4GH
    protocol end points.
    """
    def __init__(self, dataDir, variantSetClass):
        self._dataDir = dataDir
        # TODO this code is very ugly and should be regarded as a temporary
        # stop-gap until we deal with iterating over the data tree properly.
        # Variants
        self._variantSetIdMap = {}
        variantSetDir = os.path.join(self._dataDir, "variants")
        for variantSetId in os.listdir(variantSetDir):
            relativePath = os.path.join(variantSetDir, variantSetId)
            if os.path.isdir(relativePath):
                self._variantSetIdMap[variantSetId] = variantSetClass(
                    variantSetId, relativePath)
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

        # References
        self._referenceSetIdMap = {}
        referenceSetDir = os.path.join(self._dataDir, "references")
        for referenceSetId in os.listdir(referenceSetDir):
            relativePath = os.path.join(referenceSetDir, referenceSetId)
            if os.path.isdir(relativePath):
                referenceSet = references.ReferenceSet(
                    referenceSetId, relativePath)
                self._referenceSetIdMap[referenceSetId] = referenceSet
        self._referenceSetIds = sorted(self._referenceSetIdMap.keys())
        # Reads
        self._readGroupSetIdMap = {}
        readGroupSetDir = os.path.join(self._dataDir, "reads")
        for readGroupSetId in os.listdir(readGroupSetDir):
            relativePath = os.path.join(readGroupSetDir, readGroupSetId)
            if os.path.isdir(relativePath):
                readGroupSet = reads.ReadGroupSet(readGroupSetId, relativePath)
                self._readGroupSetIdMap[readGroupSetId] = readGroupSet
        self._readGroupSetIds = sorted(self._readGroupSetIdMap.keys())

        # callSets
        self._callSetIdMap = {}
        for variantSet in self._variantSetIdMap.values():
            sampleNames = variantSet.getSampleNames()
            variantSetId = variantSet.getVariantSetId()
            for name in sampleNames:
                if name not in self._callSetIdMap:
                    self._callSetIdMap[name] = []
                self._callSetIdMap[name].append(variantSetId)

    def getVariantSetIdMap(self):
        return self._variantSetIdMap

    def getCallSetIdMap(self):
        return self._callSetIdMap

    def parsePageToken(self, pageToken, numValues):
        """
        Parses the specified pageToken and returns a list of the specified
        number of values. Page tokens are assumed to consist of a fixed
        number of integers seperated by colons. If the page token does
        not conform to this specification, raise a InvalidPageToken
        exception.
        """
        tokens = pageToken.split(":")
        # TODO define exceptions.InvalidPageToken and raise here.
        if len(tokens) != numValues:
            raise Exception("Invalid number of values in page token")
        # TODO catch a ValueError here when bad integers are passed and
        # convert this into the appropriate InvalidPageToken exception.
        values = map(int, tokens)
        return values

    def runSearchRequest(
            self, requestStr, requestClass, responseClass, pageListName,
            objectGenerator):
        """
        Runs the specified request. The request is a string containing
        a JSON representation of an instance of the specified requestClass
        in which the page list variable has the specified pageListName.
        We return a string representation of an instance of the specified
        responseClass in JSON format. Objects are filled into the page list
        using the specified object generator, which must return
        (object, nextPageToken) pairs, and be able to resume iteration from
        any point using the nextPageToken attribute of the request object.
        """
        self.startProfile()
        # TODO change this to fromJsonDict and validate
        request = requestClass.fromJsonString(requestStr)
        pageList = []
        nextPageToken = None
        for obj, nextPageToken in objectGenerator(request):
            pageList.append(obj)
            if len(pageList) >= request.pageSize:
                break
        response = responseClass()
        response.nextPageToken = nextPageToken
        setattr(response, pageListName, pageList)
        self.endProfile()
        return response.toJsonString()

    def searchReadGroupSets(self, request):
        """
        Returns a GASearchReadGroupSetsResponse for the specified
        GASearchReadGroupSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.GASearchReadGroupSetsRequest,
            protocol.GASearchReadGroupSetsResponse, "readGroupSets",
            self.readGroupSetsGenerator)

    def searchReads(self, request):
        """
        Returns a GASearchReadsResponse for the specified
        GASearchReadsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.GASearchReadsRequest,
            protocol.GASearchReadsResponse, "reads",
            self.readsGenerator)

    def searchReferenceSets(self, request):
        """
        Returns a GASearchReferenceSetsResponse for the specified
        GASearchReferenceSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.GASearchReferenceSetsRequest,
            protocol.GASearchReferenceSetsResponse, "referenceSets",
            self.referenceSetsGenerator)

    def searchReferences(self, request):
        """
        Returns a GASearchReferencesResponse for the specified
        GASearchReferencesRequest object.
        """
        return self.runSearchRequest(
            request, protocol.GASearchReferencesRequest,
            protocol.GASearchReferencesResponse, "references",
            self.referencesGenerator)

    def searchVariantSets(self, request):
        """
        Returns a GASearchVariantSetsResponse for the specified
        GASearchVariantSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.GASearchVariantSetsRequest,
            protocol.GASearchVariantSetsResponse, "variantSets",
            self.variantSetsGenerator)

    def searchVariants(self, request):
        """
        Returns a GASearchVariantsResponse for the specified
        GASearchVariantsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.GASearchVariantsRequest,
            protocol.GASearchVariantsResponse, "variants",
            self.variantsGenerator)

    def searchCallSets(self, request):
        """
        Returns a GASearchCallSetsResponse for the specified
        GASearchCallSetsRequest Object.
        """
        return self.runSearchRequest(
            request, protocol.GASearchCallSetsRequest,
            protocol.GASearchCallSetsResponse, "callSets",
            self.callSetsGenerator)

    # Iterators over the data hieararchy

    def _topLevelObjectGenerator(self, request, idMap, idList):
        """
        Generalisation of the code to iterate over the objects at the top
        of the data hierarchy.
        """
        currentIndex = 0
        if request.pageToken is not None:
            currentIndex, = self.parsePageToken(request.pageToken, 1)
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

    def variantSetsGenerator(self, request):
        """
        Returns a generator over the (variantSet, nextPageToken) pairs defined
        by the specified request.
        """
        return self._topLevelObjectGenerator(
            request, self._variantSetIdMap, self._variantSetIds)

    def variantsGenerator(self, request):
        """
        Returns a generator over the (variant, nextPageToken) pairs defined by
        the specified request.
        """
        variantSetIds = request.variantSetIds
        startVariantSetIndex = 0
        startPosition = request.start
        if request.pageToken is not None:
            startVariantSetIndex, startPosition = self.parsePageToken(
                request.pageToken, 2)
        for variantSetIndex in range(startVariantSetIndex, len(variantSetIds)):
            variantSetId = variantSetIds[variantSetIndex]
            if variantSetId in self._variantSetIdMap:
                variantSet = self._variantSetIdMap[variantSetId]
                iterator = variantSet.getVariants(
                    request.referenceName, startPosition, request.end,
                    request.variantName, request.callSetIds)
                for variant in iterator:
                    nextPageToken = "{0}:{1}".format(
                        variantSetIndex, variant.start + 1)
                    yield variant, nextPageToken

    def callSetsGenerator(self, request):
        """
        Returns a generator over the (callSet, nextPageToken) pairs defined by
        the specified request.
        """
        # if no variantSetIds are input from client,
        # then set variantSetIds to all variantSetIds
        if request.variantSetIds == []:
            variantSetIds = self._variantSetIds
        else:
            variantSetIds = request.variantSetIds
        name = request.name
        if request.pageToken is not None:
            startVariantSetIndex, startCallSetPosition = self.parsePageToken(
                request.pageToken, 2)
        else:
            startVariantSetIndex = 0
            startCallSetPosition = 0

        for variantSetIndex in range(startVariantSetIndex, len(variantSetIds)):
            variantSetId = variantSetIds[variantSetIndex]
            if variantSetId in self._variantSetIdMap:
                variantSet = self._variantSetIdMap[variantSetId]
                callSetIterator = variantSet.getCallSets(
                    name, startCallSetPosition)
                for callSet, callSetPosition in callSetIterator:
                    callSet.variantSetIds = self._callSetIdMap[callSet.name]
                    nextPageToken = "{0}:{1}".format(
                        variantSetIndex, callSetPosition + 1)
                    yield callSet, nextPageToken

    def startProfile(self):
        pass

    def endProfile(self):
        pass


class MockBackend(Backend):
    """
    A mock Backend class for testing.
    """
    def __init__(self, dataDir=None):
        # TODO make a superclass of backend that does this
        # automatically without needing to know about the internal
        # details of the backend.
        self._dataDir = None
        self._variantSetIdMap = {}
        self._variantSetIds = []
