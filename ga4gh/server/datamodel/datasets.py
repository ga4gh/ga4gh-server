"""
Dataset objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.server.datamodel as datamodel
import ga4gh.server.datamodel.reads as reads
import ga4gh.server.datamodel.sequence_annotations as sequence_annotations
import ga4gh.server.datamodel.continuous as continuous
import ga4gh.server.datamodel.variants as variants
import ga4gh.server.exceptions as exceptions
import ga4gh.server.datamodel.bio_metadata as biodata
import ga4gh.server.datamodel.genotype_phenotype as g2p
import ga4gh.server.datamodel.rna_quantification as rnaQuantification

import ga4gh.schemas.pb as pb
import ga4gh.schemas.protocol as protocol


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
        self._continuousSetIds = []
        self._continuousSetIdMap = {}
        self._continuousSetNameMap = {}
        self._readGroupSetIds = []
        self._readGroupSetIdMap = {}
        self._readGroupSetNameMap = {}
        self._biosampleIds = []
        self._biosampleIdMap = {}
        self._biosampleNameMap = {}
        self._individualIds = []
        self._individualIdMap = {}
        self._individualNameMap = {}
        self._phenotypeAssociationSetIdMap = {}
        self._phenotypeAssociationSetNameMap = {}
        self._phenotypeAssociationSetIds = []
        self._rnaQuantificationSetIds = []
        self._rnaQuantificationSetIdMap = {}
        self._rnaQuantificationSetNameMap = {}

    def populateFromRow(self, dataset):
        """
        Populates the instance variables of this Dataset from the
        specified database row.
        """
        self._description = dataset.description
        self.setAttributesJson(dataset.attributes)

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

    def addBiosample(self, biosample):
        """
        Adds the specified biosample to this dataset.
        """
        id_ = biosample.getId()
        self._biosampleIdMap[id_] = biosample
        self._biosampleIds.append(id_)
        self._biosampleNameMap[biosample.getName()] = biosample

    def addIndividual(self, individual):
        """
        Adds the specified individual to this dataset.
        """
        id_ = individual.getId()
        self._individualIdMap[id_] = individual
        self._individualIds.append(id_)
        self._individualNameMap[individual.getName()] = individual

    def addFeatureSet(self, featureSet):
        """
        Adds the specified featureSet to this dataset.
        """
        id_ = featureSet.getId()
        self._featureSetIdMap[id_] = featureSet
        self._featureSetIds.append(id_)
        name = featureSet.getLocalId()
        self._featureSetNameMap[name] = featureSet

    def addContinuousSet(self, continuousSet):
        """
        Adds the specified continuousSet to this dataset.
        """
        id_ = continuousSet.getId()
        self._continuousSetIdMap[id_] = continuousSet
        self._continuousSetIds.append(id_)
        name = continuousSet.getLocalId()
        self._continuousSetNameMap[name] = continuousSet

    def addReadGroupSet(self, readGroupSet):
        """
        Adds the specified readGroupSet to this dataset.
        """
        id_ = readGroupSet.getId()
        self._readGroupSetIdMap[id_] = readGroupSet
        self._readGroupSetNameMap[readGroupSet.getLocalId()] = readGroupSet
        self._readGroupSetIds.append(id_)

    def addRnaQuantificationSet(self, rnaQuantSet):
        """
        Adds the specified rnaQuantification set to this dataset.
        """
        id_ = rnaQuantSet.getId()
        self._rnaQuantificationSetIdMap[id_] = rnaQuantSet
        self._rnaQuantificationSetIds.append(id_)
        name = rnaQuantSet.getLocalId()
        self._rnaQuantificationSetNameMap[name] = rnaQuantSet

    def toProtocolElement(self):
        dataset = protocol.Dataset()
        dataset.id = self.getId()
        dataset.name = pb.string(self.getLocalId())
        dataset.description = pb.string(self.getDescription())
        self.serializeAttributes(dataset)
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

    def addPhenotypeAssociationSet(self, phenotypeAssociationSet):
        """
        Adds the specified g2p association set to this backend.
        """
        id_ = phenotypeAssociationSet.getId()
        self._phenotypeAssociationSetIdMap[id_] = phenotypeAssociationSet
        self._phenotypeAssociationSetNameMap[
            phenotypeAssociationSet.getLocalId()] = phenotypeAssociationSet
        self._phenotypeAssociationSetIds.append(id_)

    def getPhenotypeAssociationSets(self):
        return [self._phenotypeAssociationSetIdMap[id_]
                for id_ in self._phenotypeAssociationSetIdMap]

    def getPhenotypeAssociationSet(self, id_):
        return self._phenotypeAssociationSetIdMap[id_]

    def getPhenotypeAssociationSetByName(self, name):
        if name not in self._phenotypeAssociationSetNameMap:
            raise exceptions.PhenotypeAssociationSetNotFoundException(name)
        return self._phenotypeAssociationSetNameMap[name]

    def getPhenotypeAssociationSetByIndex(self, index):
        return self._phenotypeAssociationSetIdMap[
            self._phenotypeAssociationSetIds[index]]

    def getNumPhenotypeAssociationSets(self):
        """
        Returns the number of reference sets in this data repository.
        """
        return len(self._phenotypeAssociationSetIds)

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

    def getContinuousSets(self):
        """
        Returns the list of ContinuousSets in this dataset
        """
        return [self._continuousSetIdMap[id_]
                for id_ in self._continuousSetIds]

    def getNumContinuousSets(self):
        """
        Returns the number of continuous sets in this dataset.
        """
        return len(self._continuousSetIds)

    def getContinuousSet(self, id_):
        """
        Returns the ContinuousSet with the specified id, or raises a
        ContinuousSetNotFoundException otherwise.
        """
        if id_ not in self._continuousSetIdMap:
            raise exceptions.ContinuousSetNotFoundException(id_)
        return self._continuousSetIdMap[id_]

    def getContinuousSetByName(self, name):
        """
        Returns the ContinuousSet with the specified name, or raises
        an exception otherwise.
        """
        if name not in self._continuousSetNameMap:
            raise exceptions.ContinuousSetNameNotFoundException(name)
        return self._continuousSetNameMap[name]

    def getContinuousSetByIndex(self, index):
        """
        Returns the continuous set at the specified index in this dataset.
        """
        return self._continuousSetIdMap[self._continuousSetIds[index]]

    def getBiosamples(self):
        """
        Returns the list of biosamples in this dataset
        """
        return [self._biosampleIdMap[id_] for id_ in self._biosampleIds]

    def getBiosampleByName(self, name):
        """
        Returns a Biosample with the specified name, or raises a
        BiosampleNameNotFoundException if it does not exist.
        """
        if name not in self._biosampleNameMap:
            raise exceptions.BiosampleNameNotFoundException(name)
        return self._biosampleNameMap[name]

    def getBiosample(self, id_):
        """
        Returns the Biosample with the specified id, or raises
        a BiosampleNotFoundException otherwise.
        """
        if id_ not in self._biosampleIdMap:
            raise exceptions.BiosampleNotFoundException(id_)
        return self._biosampleIdMap[id_]

    def getIndividuals(self):
        """
        Returns the list of individuals in this dataset
        """
        return [self._individualIdMap[id_] for id_ in self._individualIds]

    def getIndividualByName(self, name):
        """
        Returns an individual with the specified name, or raises a
        IndividualNameNotFoundException if it does not exist.
        """
        if name not in self._individualNameMap:
            raise exceptions.IndividualNameNotFoundException(name)
        return self._individualNameMap[name]

    def getIndividual(self, id_):
        """
        Returns the Individual with the specified id, or raises
        a IndividualNotFoundException otherwise.
        """
        if id_ not in self._individualIdMap:
            raise exceptions.IndividualNotFoundException(id_)
        return self._individualIdMap[id_]

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

    def getInfo(self):
        """
        Returns the info of this dataset.
        """
        return self._info

    def getDescription(self):
        """
        Returns the free text description of this dataset.
        """
        return self._description

    def getNumRnaQuantificationSets(self):
        """
        Returns the number of rna quantification sets in this dataset.
        """
        return len(self._rnaQuantificationSetIds)

    def getRnaQuantificationSets(self):
        """
        Returns the list of RnaQuantification sets in this dataset
        """
        return [self._rnaQuantificationSetIdMap[id_] for
                id_ in self._rnaQuantificationSetIds]

    def getRnaQuantificationSetByIndex(self, index):
        """
        Returns the rna quantification set at the specified index in this
        dataset.
        """
        return self._rnaQuantificationSetIdMap[
            self._rnaQuantificationSetIds[index]]

    def getRnaQuantificationSetByName(self, name):
        """
        Returns the RnaQuantification set with the specified name, or raises
        an exception otherwise.
        """
        if name not in self._rnaQuantificationSetNameMap:
            raise exceptions.RnaQuantificationSetNameNotFoundException(name)
        return self._rnaQuantificationSetNameMap[name]

    def getRnaQuantificationSet(self, id_):
        """
        Returns the RnaQuantification set with the specified name, or raises
        a RnaQuantificationSetNotFoundException otherwise.
        """
        if id_ not in self._rnaQuantificationSetIdMap:
            raise exceptions.RnaQuantificationSetNotFoundException(id_)
        return self._rnaQuantificationSetIdMap[id_]


class SimulatedDataset(Dataset):
    """
    A simulated dataset
    """
    def __init__(
            self, localId, referenceSet, randomSeed=0,
            numVariantSets=1, numCalls=1, variantDensity=0.5,
            numReadGroupSets=1, numReadGroupsPerReadGroupSet=1,
            numAlignments=1, numFeatureSets=1, numContinuousSets=1,
            numPhenotypeAssociationSets=1,
            numPhenotypeAssociations=2, numRnaQuantSets=2,
            numExpressionLevels=2):
        super(SimulatedDataset, self).__init__(localId)
        self._description = "Simulated dataset {}".format(localId)

        for i in range(numPhenotypeAssociationSets):
            localId = "simPas{}".format(i)
            seed = randomSeed + i
            phenotypeAssociationSet = g2p.SimulatedPhenotypeAssociationSet(
                self, localId, seed, numPhenotypeAssociations)
            self.addPhenotypeAssociationSet(phenotypeAssociationSet)

        # TODO create a simulated Ontology
        # Variants
        for i in range(numVariantSets):
            localId = "simVs{}".format(i)
            seed = randomSeed + i
            variantSet = variants.SimulatedVariantSet(
                self, referenceSet, localId, seed, numCalls, variantDensity)
            callSets = variantSet.getCallSets()
            # Add biosamples
            for callSet in callSets:
                biosample = biodata.Biosample(
                    self, callSet.getLocalId())
                biosample2 = biodata.Biosample(
                    self, callSet.getLocalId() + "2")
                individual = biodata.Individual(
                    self, callSet.getLocalId())
                biosample.setIndividualId(individual.getId())
                biosample2.setIndividualId(individual.getId())
                self.addIndividual(individual)
                self.addBiosample(biosample)
                self.addBiosample(biosample2)
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
            for rg in readGroupSet.getReadGroups():
                biosample = biodata.Biosample(
                    self, rg.getLocalId())
                individual = biodata.Individual(
                    self, rg.getLocalId())
                biosample.setIndividualId(individual.getId())
                rg.setBiosampleId(biosample.getId())
                self.addIndividual(individual)
                self.addBiosample(biosample)
            self.addReadGroupSet(readGroupSet)
        # Features
        for i in range(numFeatureSets):
            localId = "simFs{}".format(i)
            seed = randomSeed + i
            featureSet = sequence_annotations.SimulatedFeatureSet(
                self, localId, seed)
            featureSet.setReferenceSet(referenceSet)
            self.addFeatureSet(featureSet)
        # Continuous
        for i in range(numContinuousSets):
            localId = "simConts{}".format(i)
            seed = randomSeed + i
            continuousSet = continuous.SimulatedContinuousSet(
                self, localId, seed)
            continuousSet.setReferenceSet(referenceSet)
            self.addContinuousSet(continuousSet)
        # RnaQuantificationSets
        for i in range(numRnaQuantSets):
            localId = 'simRqs{}'.format(i)
            rnaQuantSet = rnaQuantification.SimulatedRnaQuantificationSet(
                self, localId)
            rnaQuantSet.setReferenceSet(referenceSet)
            self.addRnaQuantificationSet(rnaQuantSet)
