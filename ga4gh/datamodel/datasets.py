"""
Dataset objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.variants as variants
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol


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
        self._readGroupSetNameMap = {}

    def addVariantSet(self, variantSet):
        """
        Adds the specified variantSet to this dataset.
        """
        id_ = variantSet.getId()
        self._variantSetIdMap[id_] = variantSet
        self._variantSetIds.append(id_)

    def addReadGroupSet(self, readGroupSet):
        """
        Adds the specified readGroupSet to this dataset.
        """
        id_ = readGroupSet.getId()
        self._readGroupSetIdMap[id_] = readGroupSet
        self._readGroupSetNameMap[readGroupSet.getLocalId()] = readGroupSet
        self._readGroupSetIds.append(id_)

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
        return [self._variantSetIdMap[id_] for id_ in self._variantSetIds]

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

    def getReadGroupSets(self):
        """
        Returns the list of ReadGroupSets in this dataset
        """
        return self._readGroupSetIdMap.values()

    def getReadGroupSetByName(self, name):
        """
        Returns a ReadGroupSet with the specified name, or raises a
        ReadGroupSetNameNotFoundException if it does not exist.
        """
        if name not in self._readGroupSetNameMap:
            raise exceptions.ReadGroupSetNameNotFoundException(name)
        return self._readGroupSetNameMap[name]

    def getReadGroupSet(self, id_):
        """
        Returns the ReadGroupSet with the specified name, or raises
        a ReadGroupSetNotFoundException otherwise.
        """
        if id_ not in self._readGroupSetIdMap:
            raise exceptions.ReadGroupNotFoundException(id_)
        return self._readGroupSetIdMap[id_]


class SimulatedDataset(AbstractDataset):
    """
    A simulated dataset
    """
    def __init__(
            self, localId, randomSeed=0,
            numVariantSets=1, numCalls=1, variantDensity=0.5,
            numReadGroupSets=1, numReadGroupsPerReadGroupSet=1,
            numAlignments=1):
        super(SimulatedDataset, self).__init__(localId)
        # Variants
        for i in range(numVariantSets):
            localId = "simVs{}".format(i)
            seed = randomSeed + i
            variantSet = variants.SimulatedVariantSet(
                self, localId, seed, numCalls, variantDensity)
            self.addVariantSet(variantSet)
        # Reads
        for i in range(numReadGroupSets):
            localId = 'simRgs{}'.format(i)
            seed = randomSeed + i
            readGroupSet = reads.SimulatedReadGroupSet(
                self, localId, seed, numReadGroupsPerReadGroupSet,
                numAlignments)
            self.addReadGroupSet(readGroupSet)


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
                self.addVariantSet(variantSet)
        # Reads
        readGroupSetDir = os.path.join(dataDir, "reads")
        for localId in os.listdir(readGroupSetDir):
            relativePath = os.path.join(readGroupSetDir, localId)
            if os.path.isdir(relativePath):
                readGroupSet = reads.HtslibReadGroupSet(
                    self, localId, relativePath)
                self.addReadGroupSet(readGroupSet)
