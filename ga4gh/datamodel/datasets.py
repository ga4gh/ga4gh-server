"""
Module responsible for translating variant data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.datamodel.variants as variants
import ga4gh.datamodel.datacontainers as datacontainers


class DatasetContainer(datacontainers.DataContainer):
    """
    Class representing a set of datasets.
    """
    def __init__(self):
        super(DatasetContainer, self).__init__()

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

    def variantSetsGenerator(self, request):
        """
        Returns a generator over the (variantSet, nextPageToken) pairs defined
        by the speficied request.
        """
        datasetId = ""
        variantSetId = ""
        if pageToken not None:
            datasetId, variantSetId = self.parsePageToken(request.pageToken, 2)



        while currentIndex < len(self._variantSetIds):

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
class Dataset(datacontainers.FileSystemDataContainer):
    """
    Class representing a single dataset in the GA4GH data model. A
    dataset corresponds to a single directory in the file system
    that contains a set of reads and variants.
    """
    def __init__(self, dataDir):
        super(Dataset, self).__init__(dataDir, None)
        self._variantSetIdMap = {}
        # Allocate the variantSets
        variantSetsDir = os.path.join(dataDir, "variants")
        for entry in os.listdir(variantSetsDir):
            relativePath = os.path.join(variantSetsDir, entry)
            variantSet = variants.variantSetFactory(relativePath, self)
            variantSetId = variantSet.getId()
            self._variantSetIdMap[variantSetId] = variantSet
        self._variantSetIds = sorted(self._variantSetIdMap.keys())


    def searchVariantSets(self, request, lastId):

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


