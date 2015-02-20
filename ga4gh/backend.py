"""
Module responsible for handling protocol requests and returning
responses.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import bintrees

import ga4gh.protocol as protocol


class DataObjectCollection(object):
    """
    A class representing a collection of data objects allowing for
    simple and efficient paging using pageTokens.
    """
    def __init__(self):
        self._idMap = bintrees.AVLTree()

    def add(self, dataObject):
        """
        Adds the specified object to this collection.
        """
        self._idMap[dataObject.getId()] = dataObject

    def parsePageToken(self, pageToken, numValues):
        """
        Parses the specified pageToken and returns a list of the specified
        number of values. Page tokens are assumed to consist of a fixed
        number of strings seperated by colons. If the page token does
        not conform to this specification, raise a InvalidPageToken
        exception.
        """
        tokens = pageToken.split(":")
        # TODO define exceptions.InvalidPageToken and raise here.
        if len(tokens) != numValues:
            raise Exception("Invalid number of values in page token")
        return tokens

    def search(self, request, depth):
        """
        Returns an iterator over the ProtocolElements defined by
        the specified request.
        """
        if depth == 1:
            maxKey = self._idMap.max_key()
            if request.pageToken is None:
                nextId = self._idMap.min_key()
            else:
                nextId, = self.parsePageToken(request.pageToken, 1)
            for dataObject in self._idMap.value_slice(nextId, None):
                protocolElement = dataObject.toProtocolElement()
                nextId = None
                if dataObject.getId() != maxKey:
                    nextId = self._idMap.succ_key(dataObject.getId())
                yield protocolElement, nextId
        else:
            assert depth == 2
            if request.pageToken is None:
                variantSetIds = request.variantSetIds
                nextId = variantSetIds[0]
            else:
                nextId, subtoken = self.parsePageToken(request.pageToken, 2)
            for dataObject in self._idMap.value_slice(nextId, None):
                # DataObject is also a data collection. So, we iterate over this.
                for v in dataObject.search(request, 1):
                    yield v

                # if dataObject.getId() != maxKey:
                # nextId = self._idMap.succ_key(dataObject.getId())

#         variantSetIds = request.variantSetIds
#         startVariantSetIndex = 0
#         startPosition = request.start
#         if request.pageToken is not None:
#             startVariantSetIndex, startPosition = self.parsePageToken(
#                 request.pageToken, 2)
#         for variantSetIndex in range(startVariantSetIndex, len(variantSetIds)):
#             variantSetId = variantSetIds[variantSetIndex]
#             if variantSetId in self._variantSetIdMap:
#                 variantSet = self._variantSetIdMap[variantSetId]
#                 iterator = variantSet.getVariants(
#                     request.referenceName, startPosition, request.end,
#                     request.variantName, request.callSetIds)
#                 for variant in iterator:
#                     nextPageToken = "{0}:{1}".format(
#                         variantSetIndex, variant.start + 1)
#                     yield variant, nextPageToken



class Backend(object):
    """
    The GA4GH backend. This class provides methods for all of the GA4GH
    protocol end points.
    """
    def __init__(self, dataDir, variantSetClass):
        self._dataDir = dataDir
        self._variantSets = DataObjectCollection()
        # All directories in datadir are assumed to correspond to VariantSets.
        for variantSetId in os.listdir(self._dataDir):
            relativePath = os.path.join(self._dataDir, variantSetId)
            if os.path.isdir(relativePath):
                variantSet = variantSetClass(variantSetId, relativePath)
                self._variantSets.add(variantSet)


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
        by the specified request.
        """
        for protocolElement in self._variantSets.search(request, 1):
            yield protocolElement

    def variantsGenerator(self, request):
        """
        Returns a generator over the (variant, nextPageToken) pairs defined by
        the specified request.
        """
        for protocolElement in self._variantSets.search(request, 2):
            yield protocolElement

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
