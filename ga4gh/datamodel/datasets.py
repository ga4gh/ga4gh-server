"""
Dataset objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import fnmatch
import json
import os

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.variants as variants
import ga4gh.datamodel.biodata as biodata
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
        self._variantAnnotationSetIds = []
        self._variantAnnotationSetIdMap = {}
        self._bioSampleIds = []
        self._bioSampleIdMap = {}
        self._bioSampleNameMap = {}
        self._individualIds = []
        self._individualIdMap = {}
        self._individualNameMap = {}
        self._description = None

    def addBioSample(self, bioSample):
        id_ = bioSample.getId()
        self._bioSampleIdMap[id_] = bioSample
        self._bioSampleIds.append(id_)
        self._bioSampleNameMap[bioSample.getName()] = bioSample

    def addIndividual(self, individual):
        id_ = individual.getId()
        self._individualIdMap[id_] = individual
        self._individualIds.append(id_)
        self._individualNameMap[individual.getName()] = individual

    def addVariantSet(self, variantSet):
        """
        Adds the specified variantSet to this dataset.
        """
        id_ = variantSet.getId()
        self._variantSetIdMap[id_] = variantSet
        self._variantSetIds.append(id_)

    def addVariantAnnotationSet(self, variantAnnotationSet):
        """
        Adds the specified variantAnnotationSet to this dataset.
        """
        id_ = variantAnnotationSet.getId()
        self._variantAnnotationSetIdMap[id_] = variantAnnotationSet
        self._variantAnnotationSetIds.append(id_)

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
        dataset.name = self.getLocalId()
        dataset.description = self.getDescription()
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

    def getVariantAnnotationSets(self):
        """
        Returns the list of VariantAnnotationSets in this dataset
        """
        return [self._variantAnnotationSetIdMap[id_] for id_ in
                self._variantAnnotationSetIds]

    def getVariantAnnotationSet(self, id_):
        """
        Returns the AnnotationSet in this dataset with the specified 'id'
        """
        if id_ not in self._variantAnnotationSetIdMap:
            raise exceptions.AnnotationSetNotFoundException(id_)
        return self._variantAnnotationSetIdMap[id_]

    def getNumVariantAnnotationSets(self):
        """
        Returns the number of variant annotation sets in this dataset.
        """
        return len(self._variantAnnotationSetIds)

    def getVariantSet(self, id_):
        """
        Returns the VariantSet with the specified name, or raises a
        VariantSetNotFoundException otherwise.
        """
        if id_ not in self._variantSetIdMap:
            raise exceptions.VariantSetNotFoundException(id_)
        return self._variantSetIdMap[id_]

    def getBioSample(self, id_):
        """
        Returns the BioSample with the specified ID, or raises a
        BioSampleNotFoundException otherwise.
        """
        if id_ not in self._bioSampleIdMap:
            raise exceptions.BioSampleNotFoundException(id_)
        return self._bioSampleIdMap[id_]

    def getBioSamples(self):
        """
        Returns all BioSamples in this dataset
        """
        return [self._bioSampleIdMap[id_] for id_ in self._bioSampleIds]

    def getIndividual(self, id_):
        """
        Returns the Individual with the specified ID, or raises a
        IndividualNotFoundException otherwise.
        """
        if id_ not in self._individualIdMap:
            raise exceptions.IndividualNotFoundException(id_)
        return self._individualIdMap[id_]

    def getIndividuals(self):
        """
        Returns all Individuals in this dataset
        """
        return [self._individualIdMap[id_] for id_ in self._individualIds]

    def getBioSampleByIndex(self, index):
        """
        Returns the BioSample set at the specified index in this dataset.
        """
        return self._bioSampleIdMap[self._bioSampleIds[index]]

    def getIndividualByIndex(self, index):
        """
        Returns the Individual set at the specified index in this dataset.
        """
        return self._individualIdMap[self._individualIds[index]]

    def getVariantSetByIndex(self, index):
        """
        Returns the variant set at the specified index in this dataset.
        """
        return self._variantSetIdMap[self._variantSetIds[index]]

    def getVariantAnnotationSetByIndex(self, index):
        """
        Returns the variant annotation set at the specified index in this
        dataset.
        """
        return self._variantAnnotationSetIdMap[
            self._variantAnnotationSetIds[index]]

    def getNumReadGroupSets(self):
        """
        Returns the number of readgroup sets in this dataset.
        """
        return len(self._readGroupSetIds)

    def getNumBioSamples(self):
        """
        Returns the number of biosamples in this dataset.
        """
        return len(self._bioSampleIds)

    def getNumIndividuals(self):
        """
        Returns the number of individuals in this dataset.
        """
        return len(self._individualIds)

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


class SimulatedDataset(AbstractDataset):
    """
    A simulated dataset
    """
    def __init__(
            self, localId, referenceSet, randomSeed=0,
            numVariantSets=1, numCalls=1, variantDensity=0.5,
            numReadGroupSets=1, numReadGroupsPerReadGroupSet=1,
            numAlignments=1):
        super(SimulatedDataset, self).__init__(localId)
        self._description = "Simulated dataset {}".format(localId)
        # Variants
        for i in range(numVariantSets):
            localId = "simVs{}".format(i)
            seed = randomSeed + i
            variantSet = variants.SimulatedVariantSet(
                self, localId, seed, numCalls, variantDensity)
            callSets = variantSet.getCallSets()
            # Add biosamples
            for callSet in callSets:
                bioSample = datamodel.biodata.AbstractBioSample(
                    self, callSet.getLocalId(), callSet.getLocalId())
                bioSample2 = datamodel.biodata.AbstractBioSample(
                    self, callSet.getLocalId() + "2", callSet.getLocalId())
                individual = datamodel.biodata.AbstractIndividual(
                    self, callSet.getLocalId())
                self.addIndividual(individual)
                self.addBioSample(bioSample)
                self.addBioSample(bioSample2)
            self.addVariantSet(variantSet)
            variantAnnotationSet = variants.SimulatedVariantAnnotationSet(
                self, "simVas{}".format(i), variantSet)
            self.addVariantAnnotationSet(variantAnnotationSet)
        # Reads
        for i in range(numReadGroupSets):
            localId = 'simRgs{}'.format(i)
            seed = randomSeed + i
            readGroupSet = reads.SimulatedReadGroupSet(
                self, localId, referenceSet, seed,
                numReadGroupsPerReadGroupSet, numAlignments)
            for rg in readGroupSet.getReadGroups():
                bioSample = datamodel.biodata.AbstractBioSample(
                    self, rg.getLocalId(), rg.getLocalId())
                individual = datamodel.biodata.AbstractIndividual(
                    self, rg.getLocalId())
                self.addIndividual(individual)
                self.addBioSample(bioSample)
            self.addReadGroupSet(readGroupSet)


class FileSystemDataset(AbstractDataset):
    """
    A dataset based on the file system
    """
    variantsDirName = "variants"
    readsDirName = "reads"
    biodataDirName = "biodata"
    bioSamplesDirName = biodataDirName + "/biosamples"
    individualsDirName = biodataDirName + "/individuals"

    def __init__(self, localId, dataDir, dataRepository):
        super(FileSystemDataset, self).__init__(localId)
        self._dataDir = dataDir
        self._setMetadata()

        # Variants
        variantSetDir = os.path.join(dataDir, self.variantsDirName)
        for localId in os.listdir(variantSetDir):
            relativePath = os.path.join(variantSetDir, localId)
            if os.path.isdir(relativePath):
                variantSet = variants.HtslibVariantSet(
                    self, localId, relativePath, dataRepository)
                self.addVariantSet(variantSet)
            # Variant annotations sets
                if variantSet.isAnnotated(relativePath):
                    variantAnnotationSet = variants.HtslibVariantAnnotationSet(
                            self, localId, relativePath, dataRepository,
                            variantSet)
                    self.addVariantAnnotationSet(variantAnnotationSet)

        # Reads
        readGroupSetDir = os.path.join(dataDir, self.readsDirName)
        for filename in os.listdir(readGroupSetDir):
            if fnmatch.fnmatch(filename, '*.bam'):
                localId, _ = os.path.splitext(filename)
                bamPath = os.path.join(readGroupSetDir, filename)
                readGroupSet = reads.HtslibReadGroupSet(
                    self, localId, bamPath, dataRepository)
                self.addReadGroupSet(readGroupSet)

        # Biodata
        bioSamplesDir = os.path.join(dataDir, self.bioSamplesDirName)
        if os.path.exists(bioSamplesDir):
            for filename in os.listdir(bioSamplesDir):
                if fnmatch.fnmatch(filename, '*.json'):
                    filepath = os.path.join(bioSamplesDir, filename)
                    localId, _ = os.path.splitext(filename)
                    bioSample = biodata.JsonBioSample(
                        self, localId, filepath)
                    self.addBioSample(bioSample)

        individualsDir = os.path.join(dataDir, self.individualsDirName)
        if os.path.exists(individualsDir):
            for filename in os.listdir(individualsDir):
                if fnmatch.fnmatch(filename, '*.json'):
                    filepath = os.path.join(individualsDir, filename)
                    localId, _ = os.path.splitext(filename)
                    individual = biodata.JsonIndividual(
                        self, localId, filepath)
                    self.addIndividual(individual)

    def _setMetadata(self):
        metadataFileName = '{}.json'.format(self._dataDir)
        if os.path.isfile(metadataFileName):
            with open(metadataFileName) as metadataFile:
                metadata = json.load(metadataFile)
                try:
                    self._description = metadata['description']
                except KeyError as err:
                    raise exceptions.MissingDatasetMetadataException(
                        metadataFileName, str(err))
