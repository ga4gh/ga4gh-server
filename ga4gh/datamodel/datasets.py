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
    def __init__(self):
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

    def getId(self):
        """
        Return the id of the dataset
        """
        return self._id

    def getDirectory(self):
        """
        Return the path of the directory of this dataset
        """
        return self._datasetDir

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
            self, datasetId, randomSeed, numCalls,
            variantDensity, numVariantSets, numAlignments):
        super(SimulatedDataset, self).__init__()
        self._id = datasetId
        self._randomSeed = randomSeed
        self._randomGenerator = random.Random()
        self._randomGenerator.seed(self._randomSeed)

        # Variants
        for i in range(numVariantSets):
            variantSetId = "simVs{}".format(i)
            seed = self._randomGenerator.randint(0, 2**32 - 1)
            variantSet = variants.SimulatedVariantSet(
                seed, numCalls, variantDensity, variantSetId)
            self._variantSetIdMap[variantSetId] = variantSet
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

        # Reads
        readGroupSetId = "aReadGroupSet"
        readGroupSet = reads.SimulatedReadGroupSet(
            readGroupSetId, numAlignments)
        self._readGroupSetIdMap[readGroupSetId] = readGroupSet
        for readGroup in readGroupSet.getReadGroups():
            self._readGroupIdMap[readGroup.getId()] = readGroup
        self._readGroupSetIds = sorted(self._readGroupSetIdMap.keys())
        self._readGroupIds = sorted(self._readGroupIdMap.keys())


class FileSystemDataset(AbstractDataset):
    """
    A dataset based on the file system
    """
    def __init__(self, datasetDir):
        super(FileSystemDataset, self).__init__()
        self._id = os.path.basename(os.path.normpath(datasetDir))
        self._datasetDir = datasetDir

        # Variants
        variantSetDir = os.path.join(self._datasetDir, "variants")
        for variantSetId in os.listdir(variantSetDir):
            relativePath = os.path.join(variantSetDir, variantSetId)
            if os.path.isdir(relativePath):
                self._variantSetIdMap[variantSetId] = \
                    variants.HtslibVariantSet(
                        variantSetId, relativePath)
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

        # Reads
        readGroupSetDir = os.path.join(self._datasetDir, "reads")
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
