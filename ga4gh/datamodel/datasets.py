"""
Dataset objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import random

import ga4gh.protocol as protocol
import ga4gh.datamodel as datamodel
import ga4gh.datamodel.variants as variants
import ga4gh.datamodel.reads as reads


class AbstractDataset(datamodel.DatamodelObject):
    """
    The base class of datasets containing variants and reads
    """
    compoundIdClass = datamodel.DatasetCompoundId

    def __init__(self, localId):
        super(AbstractDataset, self).__init__(None, localId)
        self._variantSetIds = []
        self._variantSetIdMap = {}
        self._readGroupSetIds = []
        self._readGroupSetIdMap = {}
        self._readGroupIds = []
        self._readGroupIdMap = {}

    def toProtocolElement(self):
        dataset = protocol.Dataset()
        dataset.id = self.getId()
        return dataset

    def getVariantSetIds(self):
        """
        Return a list of ids of variant sets that this dataset has
        """
        return self._variantSetIds

    def getVariantSetIdMap(self):
        """
        Return a map of the dataset's variant set ids to variant sets
        """
        return self._variantSetIdMap

    def getVariantSets(self):
        """
        Returns the list of VariantSets in this dataset
        """
        return self._variantSetIdMap.values()

    def getReadGroupSetIds(self):
        """
        Return a list of ids of read group sets that this dataset has
        """
        return self._readGroupSetIds

    def getReadGroupSetIdMap(self):
        """
        Return a map of the dataset's read group set ids to read group sets
        """
        return self._readGroupSetIdMap

    def getReadGroupIds(self):
        """
        Return a list of ids of read groups that this dataset has
        """
        return self._readGroupIds

    def getReadGroupIdMap(self):
        """
        Return a map of the dataset's read group ids to read groups
        """
        return self._readGroupIdMap

    def getReadGroupSets(self):
        """
        Returns the list of ReadGroupSets in this dataset
        """
        return self._readGroupSetIdMap.values()


class SimulatedDataset(AbstractDataset):
    """
    A simulated dataset
    """
    def __init__(
            self, localId, randomSeed=1, numCalls=1,
            variantDensity=1, numVariantSets=1, numAlignments=1):
        super(SimulatedDataset, self).__init__(localId)
        self._randomSeed = randomSeed
        self._randomGenerator = random.Random()
        self._randomGenerator.seed(self._randomSeed)

        # Variants
        for i in range(numVariantSets):
            localId = "simVs{}".format(i)
            seed = self._randomGenerator.randint(0, 2**32 - 1)
            variantSet = variants.SimulatedVariantSet(
                self, localId, seed, numCalls, variantDensity)
            self._variantSetIdMap[variantSet.getId()] = variantSet
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

        # Reads
        localId = 'aReadGroupSet'
        readGroupSet = reads.SimulatedReadGroupSet(
            self, localId, numAlignments)
        self._readGroupSetIdMap[readGroupSet.getId()] = readGroupSet
        for readGroup in readGroupSet.getReadGroups():
            self._readGroupIdMap[readGroup.getId()] = readGroup
        self._readGroupSetIds = sorted(self._readGroupSetIdMap.keys())
        self._readGroupIds = sorted(self._readGroupIdMap.keys())


class FileSystemDataset(AbstractDataset):
    """
    A dataset based on the file system
    """
    def __init__(self, dataDir):
        localId = os.path.basename(os.path.normpath(dataDir))
        super(FileSystemDataset, self).__init__(localId)

        # Variants
        variantSetDir = os.path.join(dataDir, "variants")
        for localId in os.listdir(variantSetDir):
            relativePath = os.path.join(variantSetDir, localId)
            if os.path.isdir(relativePath):
                variantSet = variants.HtslibVariantSet(
                    self, localId, relativePath)
                self._variantSetIdMap[variantSet.getId()] = variantSet
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

        # Reads
        readGroupSetDir = os.path.join(dataDir, "reads")
        for localId in os.listdir(readGroupSetDir):
            relativePath = os.path.join(readGroupSetDir, localId)
            if os.path.isdir(relativePath):
                readGroupSet = reads.HtslibReadGroupSet(
                    self, localId, relativePath)
                self._readGroupSetIdMap[readGroupSet.getId()] = readGroupSet
                for readGroup in readGroupSet.getReadGroups():
                    self._readGroupIdMap[readGroup.getId()] = readGroup
        self._readGroupSetIds = sorted(self._readGroupSetIdMap.keys())
        self._readGroupIds = sorted(self._readGroupIdMap.keys())
