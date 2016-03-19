"""
The backing data store for the GA4GH server
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.exceptions as exceptions
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.references as references
import ga4gh.datamodel.ontologies as ontologies


class AbstractDataRepository(object):
    """
    An abstract GA4GH data repository
    """
    def __init__(self):
        self._datasetIdMap = {}
        self._datasetNameMap = {}
        self._datasetIds = []
        self._referenceSetIdMap = {}
        self._referenceSetNameMap = {}
        self._referenceSetIds = []
        self._ontologyNameMap = {}
        self._ontologyNames = []

    def addDataset(self, dataset):
        """
        Adds the specified dataset to this data repository.
        """
        id_ = dataset.getId()
        self._datasetIdMap[id_] = dataset
        self._datasetNameMap[dataset.getLocalId()] = dataset
        self._datasetIds.append(id_)

    def addReferenceSet(self, referenceSet):
        """
        Adds the specified reference set to this data repository.
        """
        id_ = referenceSet.getId()
        self._referenceSetIdMap[id_] = referenceSet
        self._referenceSetNameMap[referenceSet.getLocalId()] = referenceSet
        self._referenceSetIds.append(id_)

    def addOntology(self, ontology):
        """
        Add an ontology to this data repository.
        """
        for name in ontology.keys():
            self._ontologyNameMap[name] = ontology.get(name)
            self._ontologyNames.append(name)

    def getDatasets(self):
        """
        Returns a list of datasets in this data repository
        """
        return [self._datasetIdMap[id_] for id_ in self._datasetIds]

    def getNumDatasets(self):
        """
        Returns the number of datasets in this data repository.
        """
        return len(self._datasetIds)

    def getDataset(self, id_):
        """
        Returns a dataset with the specified ID, or raises a
        DatasetNotFoundException if it does not exist.
        """
        if id_ not in self._datasetIdMap:
            raise exceptions.DatasetNotFoundException(id_)
        return self._datasetIdMap[id_]

    def getDatasetByIndex(self, index):
        """
        Returns the dataset at the specified index.
        """
        return self._datasetIdMap[self._datasetIds[index]]

    def getDatasetByName(self, name):
        """
        Returns the dataset with the specified name.
        """
        if name not in self._datasetNameMap:
            raise exceptions.DatasetNameNotFoundException(name)
        return self._datasetNameMap[name]

    def getReferenceSets(self):
        """
        Returns the list of ReferenceSets in this data repository
        """
        return [self._referenceSetIdMap[id_] for id_ in self._referenceSetIds]

    def getNumReferenceSets(self):
        """
        Returns the number of reference sets in this data repository.
        """
        return len(self._referenceSetIds)

    def getOntology(self, name):
        """
        Returns ontologies
        """
        return self._ontologyNameMap[name]

    def getReferenceSet(self, id_):
        """
        Retuns the ReferenceSet with the specified ID, or raises a
        ReferenceSetNotFoundException if it does not exist.
        """
        if id_ not in self._referenceSetIdMap:
            raise exceptions.ReferenceSetNotFoundException(id_)
        return self._referenceSetIdMap[id_]

    def getReferenceSetByIndex(self, index):
        """
        Returns the reference set at the specified index.
        """
        return self._referenceSetIdMap[self._referenceSetIds[index]]

    def getReferenceSetByName(self, name):
        """
        Returns the reference set with the specified name.
        """
        if name not in self._referenceSetNameMap:
            raise exceptions.ReferenceSetNameNotFoundException(name)
        return self._referenceSetNameMap[name]


class EmptyDataRepository(AbstractDataRepository):
    """
    A data repository that contains no data
    """
    def __init__(self):
        super(EmptyDataRepository, self).__init__()


class SimulatedDataRepository(AbstractDataRepository):
    """
    A data repository that is simulated
    """
    def __init__(
            self, randomSeed=0, numDatasets=2,
            numVariantSets=1, numCalls=2, variantDensity=0.5,
            numReferenceSets=1, numReferencesPerReferenceSet=1,
            numReadGroupSets=1, numReadGroupsPerReadGroupSet=2,
            numAlignments=2):
        super(SimulatedDataRepository, self).__init__()

        # References
        for i in range(numReferenceSets):
            localId = "referenceSet{}".format(i)
            seed = randomSeed + i
            referenceSet = references.SimulatedReferenceSet(
                localId, seed, numReferencesPerReferenceSet)
            self.addReferenceSet(referenceSet)

        # Datasets
        for i in range(numDatasets):
            seed = randomSeed + i
            localId = "simulatedDataset{}".format(i)
            referenceSet = self.getReferenceSetByIndex(i % numReferenceSets)
            dataset = datasets.SimulatedDataset(
                localId, referenceSet=referenceSet, randomSeed=seed,
                numCalls=numCalls, variantDensity=variantDensity,
                numVariantSets=numVariantSets,
                numReadGroupSets=numReadGroupSets,
                numReadGroupsPerReadGroupSet=numReadGroupsPerReadGroupSet,
                numAlignments=numAlignments)
            self.addDataset(dataset)


class FileSystemDataRepository(AbstractDataRepository):
    """
    A data repository based on the file system
    """
    referenceSetsDirName = "referenceSets"
    datasetsDirName = "datasets"
    ontologiesDirName = "ontologies"

    def __init__(self, dataDir):
        super(FileSystemDataRepository, self).__init__()
        self._dataDir = dataDir
        sourceDirNames = [self.referenceSetsDirName,
                          self.ontologiesDirName,
                          self.datasetsDirName]
        constructors = [
            references.HtslibReferenceSet, ontologies.FileSystemOntologies,
            datasets.FileSystemDataset]
        objectAdders = [self.addReferenceSet, self.addOntology,
                        self.addDataset]
        for sourceDirName, constructor, objectAdder in zip(
                sourceDirNames, constructors, objectAdders):
            sourceDir = os.path.join(self._dataDir, sourceDirName)
            for setName in os.listdir(sourceDir):
                relativePath = os.path.join(sourceDir, setName)
                if os.path.isdir(relativePath):
                    objectAdder(constructor(setName, relativePath, self))

    def checkConsistency(self):
        """
        Perform checks that ensure the consistency of the data.
        Factored into a separate method from server init since the
        data repo object can be created on a partially-complete
        data set.
        """
        for dataset in self.getDatasets():
            for readGroupSet in dataset.getReadGroupSets():
                readGroupSet.checkConsistency(self)
            for variantSet in dataset.getVariantSets():
                variantSet.checkConsistency()
