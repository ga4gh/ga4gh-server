"""
Module responsible for handling protocol requests and returning
responses.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import ga4gh.protocol as protocol
import ga4gh.datamodel.datasets as datasets


class Backend(object):
    """
    The GA4GH backend. This class provides methods for all of the GA4GH
    protocol end points.
    """
    def __init__(self):
        self._datasets = datasets.DatasetContainer()

    def addDataset(self, dataDir):
        """
        Adds a new dataset for data stored in the specifed directory.
        """
        dataset = datasets.Dataset(dataDir)
        self._datasets.addChild(dataset)


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
        # TODO change this to fromJSONDict and validate
        request = requestClass.fromJSONString(requestStr)
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
        return response.toJSONString()

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
            self._datasets.variantsGenerator)


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
