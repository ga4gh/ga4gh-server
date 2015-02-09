"""
Module responsible for handling protocol requests and returning
responses.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.protocol as protocol


class Backend(object):
    """
    The GA4GH backend. This class provides methods for all of the GA4GH
    protocol end points.
    """
    def __init__(self, dataDir, variantSetClass):
        self._dataDir = dataDir
        self._variantSetIdMap = {}
        # All directories in datadir are assumed to correspond to VariantSets.
        for variantSetId in os.listdir(self._dataDir):
            relativePath = os.path.join(self._dataDir, variantSetId)
            if os.path.isdir(relativePath):
                self._variantSetIdMap[variantSetId] = variantSetClass(
                    variantSetId, relativePath)
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

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

    def variantSetsGenerator(self, request):
        """
        Returns a generator over the (variantSet, nextPageToken) pairs defined
        by the speficied request.
        """
        currentIndex = 0
        if request.pageToken is not None:
            currentIndex, = self.parsePageToken(request.pageToken, 1)
        while currentIndex < len(self._variantSetIds):
            variantSet = protocol.GAVariantSet()
            variantSet.id = self._variantSetIds[currentIndex]
            variantSet.datasetId = "NotImplemented"
            variantSet.metadata = self._variantSetIdMap[
                variantSet.id].getMetadata()
            currentIndex += 1
            nextPageToken = None
            if currentIndex < len(self._variantSetIds):
                nextPageToken = str(currentIndex)
            yield variantSet, nextPageToken

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
