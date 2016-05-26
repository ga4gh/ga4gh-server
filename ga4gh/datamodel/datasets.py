"""
Dataset objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.sequenceAnnotations as sequenceAnnotations
import ga4gh.datamodel.variants as variants
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol
from ga4gh import pb


class Dataset(datamodel.DatamodelObject):
    """
    The base class of datasets containing variants and reads
    """
    compoundIdClass = datamodel.DatasetCompoundId

    def __init__(self, localId):
        super(Dataset, self).__init__(None, localId)
        self._description = None
        self._variantSetIds = []
        self._variantSetIdMap = {}
        self._variantSetNameMap = {}
        self._featureSetIds = []
        self._featureSetIdMap = {}
        self._featureSetNameMap = {}
        self._readGroupSetIds = []
        self._readGroupSetIdMap = {}
        self._readGroupSetNameMap = {}

    def populateFromRow(self, row):
        """
        Populates the instance variables of this Dataset from the
        specified database row.
        """
        self._description = row[b'description']

    def setDescription(self, description):
        """
        Sets the description for this dataset to the specified value.
        """
        self._description = description

    def addVariantSet(self, variantSet):
        """
        Adds the specified variantSet to this dataset.
        """
        id_ = variantSet.getId()
        self._variantSetIdMap[id_] = variantSet
        self._variantSetNameMap[variantSet.getLocalId()] = variantSet
        self._variantSetIds.append(id_)

    def addFeatureSet(self, featureSet):
        """
        Adds the specified variantSet to this dataset.
        """
        id_ = featureSet.getId()
        self._featureSetIdMap[id_] = featureSet
        self._featureSetIds.append(id_)
        name = featureSet.getLocalId()
        self._featureSetNameMap[name] = featureSet

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
        dataset.name = pb.string(self.getLocalId())
        dataset.description = pb.string(self.getDescription())
        return dataset

    def getVariantSets(self):
        """
        Returns the list of VariantSets in this dataset
        """
        return [self._variantSetIdMap[id_] for id_ in self._variantSetIds]

    def getNumVariantSets(self):
        """
        Returns the number of variant sets in this dataset.
        """
        return len(self._variantSetIds)

    def getVariantSet(self, id_):
        """
        Returns the VariantSet with the specified name, or raises a
        VariantSetNotFoundException otherwise.
        """
        if id_ not in self._variantSetIdMap:
            raise exceptions.VariantSetNotFoundException(id_)
        return self._variantSetIdMap[id_]

    def getVariantSetByIndex(self, index):
        """
        Returns the variant set at the specified index in this dataset.
        """
        return self._variantSetIdMap[self._variantSetIds[index]]

    def getVariantSetByName(self, name):
        """
        Returns a VariantSet with the specified name, or raises a
        VariantSetNameNotFoundException if it does not exist.
        """
        if name not in self._variantSetNameMap:
            raise exceptions.VariantSetNameNotFoundException(name)
        return self._variantSetNameMap[name]

    def getFeatureSets(self):
        """
        Returns the list of FeatureSets in this dataset
        """
        return [self._featureSetIdMap[id_] for id_ in self._featureSetIds]

    def getNumFeatureSets(self):
        """
        Returns the number of feature sets in this dataset.
        """
        return len(self._featureSetIds)

    def getFeatureSet(self, id_):
        """
        Returns the FeatureSet with the specified id, or raises a
        FeatureSetNotFoundException otherwise.
        """
        if id_ not in self._featureSetIdMap:
            raise exceptions.FeatureSetNotFoundException(id_)
        return self._featureSetIdMap[id_]

    def getFeatureSetByName(self, name):
        """
        Returns the FeatureSet with the specified name, or raises
        an exception otherwise.
        """
        if name not in self._featureSetNameMap:
            raise exceptions.FeatureSetNameNotFoundException(name)
        return self._featureSetNameMap[name]

    def getFeatureSetByIndex(self, index):
        """
        Returns the feature set at the specified index in this dataset.
        """
        return self._featureSetIdMap[self._featureSetIds[index]]

    def getNumReadGroupSets(self):
        """
        Returns the number of readgroup sets in this dataset.
        """
        return len(self._readGroupSetIds)

    def getReadGroupSets(self):
        """
        Returns the list of ReadGroupSets in this dataset
        """
        return [self._readGroupSetIdMap[id_] for id_ in self._readGroupSetIds]

    def getReadGroupSetByName(self, name):
        """
        Returns a ReadGroupSet with the specified name, or raises a
        ReadGroupSetNameNotFoundException if it does not exist.
        """
        if name not in self._readGroupSetNameMap:
            raise exceptions.ReadGroupSetNameNotFoundException(name)
        return self._readGroupSetNameMap[name]

    def getReadGroupSetByIndex(self, index):
        """
        Returns the readgroup set at the specified index in this dataset.
        """
        return self._readGroupSetIdMap[self._readGroupSetIds[index]]

    def getReadGroupSet(self, id_):
        """
        Returns the ReadGroupSet with the specified name, or raises
        a ReadGroupSetNotFoundException otherwise.
        """
        if id_ not in self._readGroupSetIdMap:
            raise exceptions.ReadGroupNotFoundException(id_)
        return self._readGroupSetIdMap[id_]

    def getDescription(self):
        """
        Returns the free text description of this dataset.
        """
        return self._description


class SimulatedDataset(Dataset):
    """
    A simulated dataset
    """
    def __init__(
            self, localId, referenceSet, randomSeed=0,
            numVariantSets=1, numCalls=1, variantDensity=0.5,
            numReadGroupSets=1, numReadGroupsPerReadGroupSet=1,
            numAlignments=1, numFeatureSets=1):
        super(SimulatedDataset, self).__init__(localId)
        self._description = "Simulated dataset {}".format(localId)
        # TODO create a simulated Ontology
        # Variants
        for i in range(numVariantSets):
            localId = "simVs{}".format(i)
            seed = randomSeed + i
            variantSet = variants.SimulatedVariantSet(
                self, referenceSet, localId, seed, numCalls, variantDensity)
            self.addVariantSet(variantSet)
            variantAnnotationSet = variants.SimulatedVariantAnnotationSet(
                variantSet, "simVas{}".format(i), seed)
            variantSet.addVariantAnnotationSet(variantAnnotationSet)
        # Reads
        for i in range(numReadGroupSets):
            localId = 'simRgs{}'.format(i)
            seed = randomSeed + i
            readGroupSet = reads.SimulatedReadGroupSet(
                self, localId, referenceSet, seed,
                numReadGroupsPerReadGroupSet, numAlignments)
            self.addReadGroupSet(readGroupSet)
        # Features
        for i in range(numFeatureSets):
            localId = "simFs{}".format(i)
            seed = randomSeed + i
            featureSet = sequenceAnnotations.SimulatedFeatureSet(
                self, localId, seed)
            featureSet.setReferenceSet(referenceSet)
            self.addFeatureSet(featureSet)
