"""
repo manager cli
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import glob
import json
import os
import sys
import textwrap
import traceback
import urlparse

import ga4gh.server.cli as cli
import ga4gh.server.datamodel.bio_metadata as bio_metadata
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.datamodel.genotype_phenotype as genotype_phenotype
import ga4gh.server.datamodel.ontologies as ontologies
import ga4gh.server.datamodel.reads as reads
import ga4gh.server.datamodel.references as references
import ga4gh.server.datamodel.rna_quantification as rna_quantification
import ga4gh.server.datamodel.sequence_annotations as sequence_annotations
import ga4gh.server.datamodel.continuous as continuous
import ga4gh.server.datamodel.variants as variants
import ga4gh.server.datamodel.peers as peers
import ga4gh.server.datarepo as datarepo
import ga4gh.server.exceptions as exceptions
import ga4gh.server.repo.rnaseq2ga as rnaseq2ga

import ga4gh.common.cli as common_cli


def getNameFromPath(filePath):
    """
    Returns the filename of the specified path without its extensions.
    This is usually how we derive the default name for a given object.
    """
    if len(filePath) == 0:
        raise ValueError("Cannot have empty path for name")
    fileName = os.path.split(os.path.normpath(filePath))[1]
    # We need to handle things like .fa.gz, so we can't use
    # os.path.splitext
    ret = fileName.split(".")[0]
    assert ret != ""
    return ret


def getRawInput(display):
    """
    Wrapper around raw_input; put into separate function so that it
    can be easily mocked for tests.
    """
    return raw_input(display)


class RepoManager(object):
    """
    Class that provide command line functionality to manage a
    data repository.
    """
    def __init__(self, args):
        self._args = args
        self._registryPath = args.registryPath
        self._repo = datarepo.SqlDataRepository(self._registryPath)

    def _confirmDelete(self, objectType, name, func):
        if self._args.force:
            func()
        else:
            displayString = (
                "Are you sure you want to delete the {} '{}'? "
                "[y|N] ".format(objectType, name))
            userResponse = getRawInput(displayString)
            if userResponse.strip() == 'y':
                func()
            else:
                print("Aborted")

    def _updateRepo(self, func, *args, **kwargs):
        """
        Runs the specified function that updates the repo with the specified
        arguments. This method ensures that all updates are transactional,
        so that if any part of the update fails no changes are made to the
        repo.
        """
        # TODO how do we make this properly transactional?
        self._repo.open(datarepo.MODE_WRITE)
        try:
            func(*args, **kwargs)
            self._repo.commit()
        finally:
            self._repo.close()

    def _openRepo(self):
        if not self._repo.exists():
            raise exceptions.RepoManagerException(
                "Repo '{}' does not exist. Please create a new repo "
                "using the 'init' command.".format(self._registryPath))
        self._repo.open(datarepo.MODE_READ)

    def _checkSequenceOntology(self, ontology):
        so = ontologies.SEQUENCE_ONTOLOGY_PREFIX
        if ontology.getOntologyPrefix() != so:
            raise exceptions.RepoManagerException(
                "Ontology '{}' does not have ontology prefix '{}'".format(
                    ontology.getName(), so))

    def _getFilePath(self, filePath, useRelativePath):
        return filePath if useRelativePath else os.path.abspath(filePath)

    def init(self):
        forceMessage = (
            "Respository '{}' already exists. Use --force to overwrite")
        if self._repo.exists():
            if self._args.force:
                self._repo.delete()
            else:
                raise exceptions.RepoManagerException(
                    forceMessage.format(self._registryPath))
        self._updateRepo(self._repo.initialise)

    def list(self):
        """
        Lists the contents of this repo.
        """
        self._openRepo()
        # TODO this is _very_ crude. We need much more options and detail here.
        self._repo.printSummary()

    def listAnnouncements(self):
        """
        Lists all the announcements the repo has received.
        """
        self._openRepo()
        self._repo.printAnnouncements()

    def clearAnnouncements(self):
        """
        Clears the list of announcements from the repo.
        """
        self._openRepo()
        self._repo.clearAnnouncements()

    def verify(self):
        """
        Checks that the data pointed to in the repository works and
        we don't have any broken URLs, missing files, etc.
        """
        self._openRepo()
        self._repo.verify()

    def addOntology(self):
        """
        Adds a new Ontology to this repo.
        """
        self._openRepo()
        name = self._args.name
        filePath = self._getFilePath(self._args.filePath,
                                     self._args.relativePath)
        if name is None:
            name = getNameFromPath(filePath)
        ontology = ontologies.Ontology(name)
        ontology.populateFromFile(filePath)
        self._updateRepo(self._repo.insertOntology, ontology)

    def addDataset(self):
        """
        Adds a new dataset into this repo.
        """
        self._openRepo()
        dataset = datasets.Dataset(self._args.datasetName)
        dataset.setDescription(self._args.description)
        dataset.setAttributes(json.loads(self._args.attributes))
        self._updateRepo(self._repo.insertDataset, dataset)

    def addReferenceSet(self):
        """
        Adds a new reference set into this repo.
        """
        self._openRepo()
        name = self._args.name
        filePath = self._getFilePath(self._args.filePath,
                                     self._args.relativePath)
        if name is None:
            name = getNameFromPath(self._args.filePath)
        referenceSet = references.HtslibReferenceSet(name)
        referenceSet.populateFromFile(filePath)
        referenceSet.setDescription(self._args.description)
        if self._args.species is not None:
            referenceSet.setSpeciesFromJson(self._args.species)
        referenceSet.setIsDerived(self._args.isDerived)
        referenceSet.setAssemblyId(self._args.assemblyId)
        referenceSet.setAttributes(json.loads(self._args.attributes))
        sourceAccessions = []
        if self._args.sourceAccessions is not None:
            sourceAccessions = self._args.sourceAccessions.split(",")
        referenceSet.setSourceAccessions(sourceAccessions)
        referenceSet.setSourceUri(self._args.sourceUri)
        self._updateRepo(self._repo.insertReferenceSet, referenceSet)

    def addReadGroupSet(self):
        """
        Adds a new ReadGroupSet into this repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        dataUrl = self._args.dataFile
        indexFile = self._args.indexFile
        parsed = urlparse.urlparse(dataUrl)
        # TODO, add https support and others when they have been
        # tested.
        if parsed.scheme in ['http', 'ftp']:
            if indexFile is None:
                raise exceptions.MissingIndexException(dataUrl)
        else:
            if indexFile is None:
                indexFile = dataUrl + ".bai"
            dataUrl = self._getFilePath(self._args.dataFile,
                                        self._args.relativePath)
            indexFile = self._getFilePath(indexFile, self._args.relativePath)
        name = self._args.name
        if self._args.name is None:
            name = getNameFromPath(dataUrl)
        readGroupSet = reads.HtslibReadGroupSet(dataset, name)
        readGroupSet.populateFromFile(dataUrl, indexFile)
        referenceSetName = self._args.referenceSetName
        if referenceSetName is None:
            # Try to find a reference set name from the BAM header.
            referenceSetName = readGroupSet.getBamHeaderReferenceSetName()
        referenceSet = self._repo.getReferenceSetByName(referenceSetName)
        readGroupSet.setReferenceSet(referenceSet)
        readGroupSet.setAttributes(json.loads(self._args.attributes))
        self._updateRepo(self._repo.insertReadGroupSet, readGroupSet)

    def addVariantSet(self):
        """
        Adds a new VariantSet into this repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        dataUrls = self._args.dataFiles
        name = self._args.name
        if len(dataUrls) == 1:
            if self._args.name is None:
                name = getNameFromPath(dataUrls[0])
            if os.path.isdir(dataUrls[0]):
                # Read in the VCF files from the directory.
                # TODO support uncompressed VCF and BCF files
                vcfDir = dataUrls[0]
                pattern = os.path.join(vcfDir, "*.vcf.gz")
                dataUrls = glob.glob(pattern)
                if len(dataUrls) == 0:
                    raise exceptions.RepoManagerException(
                        "Cannot find any VCF files in the directory "
                        "'{}'.".format(vcfDir))
                dataUrls[0] = self._getFilePath(dataUrls[0],
                                                self._args.relativePath)
        elif self._args.name is None:
            raise exceptions.RepoManagerException(
                "Cannot infer the intended name of the VariantSet when "
                "more than one VCF file is provided. Please provide a "
                "name argument using --name.")
        parsed = urlparse.urlparse(dataUrls[0])
        if parsed.scheme not in ['http', 'ftp']:
            dataUrls = map(lambda url: self._getFilePath(
                url, self._args.relativePath), dataUrls)
        # Now, get the index files for the data files that we've now obtained.
        indexFiles = self._args.indexFiles
        if indexFiles is None:
            # First check if all the paths exist locally, as they must
            # if we are making a default index path.
            for dataUrl in dataUrls:
                if not os.path.exists(dataUrl):
                    raise exceptions.MissingIndexException(
                        "Cannot find file '{}'. All variant files must be "
                        "stored locally if the default index location is "
                        "used. If you are trying to create a VariantSet "
                        "based on remote URLs, please download the index "
                        "files to the local file system and provide them "
                        "with the --indexFiles argument".format(dataUrl))
            # We assume that the indexes are made by adding .tbi
            indexSuffix = ".tbi"
            # TODO support BCF input properly here by adding .csi
            indexFiles = [filename + indexSuffix for filename in dataUrls]
        indexFiles = map(lambda url: self._getFilePath(
            url, self._args.relativePath), indexFiles)
        variantSet = variants.HtslibVariantSet(dataset, name)
        variantSet.populateFromFile(dataUrls, indexFiles)
        # Get the reference set that is associated with the variant set.
        referenceSetName = self._args.referenceSetName
        if referenceSetName is None:
            # Try to find a reference set name from the VCF header.
            referenceSetName = variantSet.getVcfHeaderReferenceSetName()
        if referenceSetName is None:
            raise exceptions.RepoManagerException(
                "Cannot infer the ReferenceSet from the VCF header. Please "
                "specify the ReferenceSet to associate with this "
                "VariantSet using the --referenceSetName option")
        referenceSet = self._repo.getReferenceSetByName(referenceSetName)
        variantSet.setReferenceSet(referenceSet)
        variantSet.setAttributes(json.loads(self._args.attributes))
        # Now check for annotations
        annotationSets = []
        if variantSet.isAnnotated() and self._args.addAnnotationSets:
            ontologyName = self._args.ontologyName
            if ontologyName is None:
                raise exceptions.RepoManagerException(
                    "A sequence ontology name must be provided")
            ontology = self._repo.getOntologyByName(ontologyName)
            self._checkSequenceOntology(ontology)
            for annotationSet in variantSet.getVariantAnnotationSets():
                annotationSet.setOntology(ontology)
                annotationSets.append(annotationSet)

        # Add the annotation sets and the variant set as an atomic update
        def updateRepo():
            self._repo.insertVariantSet(variantSet)
            for annotationSet in annotationSets:
                self._repo.insertVariantAnnotationSet(annotationSet)
        self._updateRepo(updateRepo)

    def addPhenotypeAssociationSet(self):
        """
        Adds a new phenotype association set to this repo.
        """
        self._openRepo()
        name = self._args.name
        if name is None:
            name = getNameFromPath(self._args.dirPath)
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        phenotypeAssociationSet = \
            genotype_phenotype.RdfPhenotypeAssociationSet(
                dataset, name, self._args.dirPath)
        phenotypeAssociationSet.setAttributes(
            json.loads(self._args.attributes))
        self._updateRepo(
            self._repo.insertPhenotypeAssociationSet,
            phenotypeAssociationSet)

    def removePhenotypeAssociationSet(self):
        """
        Removes a phenotype association set from the repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        phenotypeAssociationSet = dataset.getPhenotypeAssociationSetByName(
            self._args.name)

        def func():
            self._updateRepo(
                self._repo.removePhenotypeAssociationSet,
                phenotypeAssociationSet)
        self._confirmDelete(
            "PhenotypeAssociationSet",
            phenotypeAssociationSet.getLocalId(),
            func)

    def removeReferenceSet(self):
        """
        Removes a referenceSet from the repo.
        """
        self._openRepo()
        referenceSet = self._repo.getReferenceSetByName(
            self._args.referenceSetName)

        def func():
            self._updateRepo(self._repo.removeReferenceSet, referenceSet)
        self._confirmDelete("ReferenceSet", referenceSet.getLocalId(), func)

    def removeReadGroupSet(self):
        """
        Removes a readGroupSet from the repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        readGroupSet = dataset.getReadGroupSetByName(
            self._args.readGroupSetName)

        def func():
            self._updateRepo(self._repo.removeReadGroupSet, readGroupSet)
        self._confirmDelete("ReadGroupSet", readGroupSet.getLocalId(), func)

    def removeVariantSet(self):
        """
        Removes a variantSet from the repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        variantSet = dataset.getVariantSetByName(self._args.variantSetName)

        def func():
            self._updateRepo(self._repo.removeVariantSet, variantSet)
        self._confirmDelete("VariantSet", variantSet.getLocalId(), func)

    def removeDataset(self):
        """
        Removes a dataset from the repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)

        def func():
            self._updateRepo(self._repo.removeDataset, dataset)
        self._confirmDelete("Dataset", dataset.getLocalId(), func)

    def addFeatureSet(self):
        """
        Adds a new feature set into this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        filePath = self._getFilePath(self._args.filePath,
                                     self._args.relativePath)
        name = getNameFromPath(self._args.filePath)
        featureSet = sequence_annotations.Gff3DbFeatureSet(
            dataset, name)
        referenceSetName = self._args.referenceSetName
        if referenceSetName is None:
            raise exceptions.RepoManagerException(
                "A reference set name must be provided")
        referenceSet = self._repo.getReferenceSetByName(referenceSetName)
        featureSet.setReferenceSet(referenceSet)
        ontologyName = self._args.ontologyName
        if ontologyName is None:
            raise exceptions.RepoManagerException(
                "A sequence ontology name must be provided")
        ontology = self._repo.getOntologyByName(ontologyName)
        self._checkSequenceOntology(ontology)
        featureSet.setOntology(ontology)
        featureSet.populateFromFile(filePath)
        featureSet.setAttributes(json.loads(self._args.attributes))
        self._updateRepo(self._repo.insertFeatureSet, featureSet)

    def removeFeatureSet(self):
        """
        Removes a feature set from this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        featureSet = dataset.getFeatureSetByName(self._args.featureSetName)

        def func():
            self._updateRepo(self._repo.removeFeatureSet, featureSet)
        self._confirmDelete("FeatureSet", featureSet.getLocalId(), func)

    def addContinuousSet(self):
        """
        Adds a new continuous set into this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        filePath = self._getFilePath(self._args.filePath,
                                     self._args.relativePath)
        name = getNameFromPath(self._args.filePath)
        continuousSet = continuous.FileContinuousSet(dataset, name)
        referenceSetName = self._args.referenceSetName
        if referenceSetName is None:
            raise exceptions.RepoManagerException(
                "A reference set name must be provided")
        referenceSet = self._repo.getReferenceSetByName(referenceSetName)
        continuousSet.setReferenceSet(referenceSet)
        continuousSet.populateFromFile(filePath)
        self._updateRepo(self._repo.insertContinuousSet, continuousSet)

    def removeContinuousSet(self):
        """
        Removes a continuous set from this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        continuousSet = dataset.getContinuousSetByName(
                            self._args.continuousSetName)

        def func():
            self._updateRepo(self._repo.removeContinuousSet, continuousSet)
        self._confirmDelete("ContinuousSet", continuousSet.getLocalId(), func)

    def addBiosample(self):
        """
        Adds a new biosample into this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        biosample = bio_metadata.Biosample(
            dataset, self._args.biosampleName)
        biosample.populateFromJson(self._args.biosample)
        self._updateRepo(self._repo.insertBiosample, biosample)

    def removeBiosample(self):
        """
        Removes a biosample from this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        biosample = dataset.getBiosampleByName(self._args.biosampleName)

        def func():
            self._updateRepo(self._repo.removeBiosample, biosample)
        self._confirmDelete("Biosample", biosample.getLocalId(), func)

    def addIndividual(self):
        """
        Adds a new individual into this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        individual = bio_metadata.Individual(
            dataset, self._args.individualName)
        individual.populateFromJson(self._args.individual)
        self._updateRepo(self._repo.insertIndividual, individual)

    def removeIndividual(self):
        """
        Removes an individual from this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        individual = dataset.getIndividualByName(self._args.individualName)

        def func():
            self._updateRepo(self._repo.removeIndividual, individual)
        self._confirmDelete("Individual", individual.getLocalId(), func)

    def addPeer(self):
        """
        Adds a new peer into this repo
        """
        self._openRepo()
        try:
            peer = peers.Peer(
                self._args.url, json.loads(self._args.attributes))
        except exceptions.BadUrlException:
            raise exceptions.RepoManagerException("The URL for the peer was "
                                                  "malformed.")
        except ValueError as e:
            raise exceptions.RepoManagerException(
                "The attributes message "
                "was malformed. {}".format(e))
        self._updateRepo(self._repo.insertPeer, peer)

    def removePeer(self):
        """
        Removes a peer by URL from this repo
        """
        self._openRepo()

        def func():
            self._updateRepo(self._repo.removePeer, self._args.url)
        self._confirmDelete("Peer", self._args.url, func)

    def removeOntology(self):
        """
        Removes an ontology from the repo.
        """
        self._openRepo()
        ontology = self._repo.getOntologyByName(self._args.ontologyName)

        def func():
            self._updateRepo(self._repo.removeOntology, ontology)
        self._confirmDelete("Ontology", ontology.getName(), func)

    def addRnaQuantification(self):
        """
        Adds an rnaQuantification into this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        biosampleId = ""
        if self._args.biosampleName:
            biosample = dataset.getBiosampleByName(self._args.biosampleName)
            biosampleId = biosample.getId()
        if self._args.name is None:
            name = getNameFromPath(self._args.quantificationFilePath)
        else:
            name = self._args.name
        # TODO: programs not fully supported by GA4GH yet
        programs = ""
        featureType = "gene"
        if self._args.transcript:
            featureType = "transcript"
        rnaseq2ga.rnaseq2ga(
            self._args.quantificationFilePath, self._args.filePath, name,
            self._args.format, dataset=dataset, featureType=featureType,
            description=self._args.description, programs=programs,
            featureSetNames=self._args.featureSetNames,
            readGroupSetNames=self._args.readGroupSetName,
            biosampleId=biosampleId)

    def initRnaQuantificationSet(self):
        """
        Initialize an empty RNA quantification set
        """
        store = rnaseq2ga.RnaSqliteStore(self._args.filePath)
        store.createTables()

    def addRnaQuantificationSet(self):
        """
        Adds an rnaQuantificationSet into this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        if self._args.name is None:
            name = getNameFromPath(self._args.filePath)
        else:
            name = self._args.name
        rnaQuantificationSet = rna_quantification.SqliteRnaQuantificationSet(
            dataset, name)
        referenceSetName = self._args.referenceSetName
        if referenceSetName is None:
            raise exceptions.RepoManagerException(
                "A reference set name must be provided")
        referenceSet = self._repo.getReferenceSetByName(referenceSetName)
        rnaQuantificationSet.setReferenceSet(referenceSet)
        rnaQuantificationSet.populateFromFile(self._args.filePath)
        rnaQuantificationSet.setAttributes(json.loads(self._args.attributes))
        self._updateRepo(
            self._repo.insertRnaQuantificationSet, rnaQuantificationSet)

    def removeRnaQuantificationSet(self):
        """
        Removes an rnaQuantificationSet from this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        rnaQuantSet = dataset.getRnaQuantificationSetByName(
            self._args.rnaQuantificationSetName)

        def func():
            self._updateRepo(self._repo.removeRnaQuantificationSet,
                             rnaQuantSet)
        self._confirmDelete(
            "RnaQuantificationSet", rnaQuantSet.getLocalId(), func)

    #
    # Methods to simplify adding common arguments to the parser.
    #

    @classmethod
    def addRepoArgument(cls, subparser):
        subparser.add_argument(
            "registryPath",
            help="the location of the registry database")

    @classmethod
    def addForceOption(cls, subparser):
        subparser.add_argument(
            "-f", "--force", action='store_true',
            default=False, help="do not prompt for confirmation")

    @classmethod
    def addRelativePathOption(cls, subparser):
        subparser.add_argument(
            "-r", "--relativePath", action='store_true',
            default=False, help="store relative path in database")

    @classmethod
    def addDescriptionOption(cls, subparser, objectType):
        subparser.add_argument(
            "-d", "--description", default="",
            help="The human-readable description of the {}.".format(
                objectType))

    @classmethod
    def addDatasetNameArgument(cls, subparser):
        subparser.add_argument(
            "datasetName", help="the name of the dataset")

    @classmethod
    def addAttributesArgument(cls, subparser):
        subparser.add_argument(
            "-A", "--attributes", default="{}",
            help="additional attributes for the message expressed as JSON")

    @classmethod
    def addReferenceSetNameOption(cls, subparser, objectType):
        helpText = (
            "the name of the reference set to associate with this {}"
        ).format(objectType)
        subparser.add_argument(
            "-R", "--referenceSetName", default=None, help=helpText)

    @classmethod
    def addSequenceOntologyNameOption(cls, subparser, objectType):
        helpText = (
            "the name of the sequence ontology instance used to "
            "translate ontology term names to IDs in this {}"
        ).format(objectType)
        subparser.add_argument(
            "-O", "--ontologyName", default=None, help=helpText)

    @classmethod
    def addOntologyNameArgument(cls, subparser):
        subparser.add_argument(
            "ontologyName",
            help="the name of the ontology")

    @classmethod
    def addUrlArgument(cls, subparser):
        subparser.add_argument(
            "url",
            help="The URL of the given resource")

    @classmethod
    def addReadGroupSetNameArgument(cls, subparser):
        subparser.add_argument(
            "readGroupSetName",
            help="the name of the read group set")

    @classmethod
    def addVariantSetNameArgument(cls, subparser):
        subparser.add_argument(
            "variantSetName",
            help="the name of the variant set")

    @classmethod
    def addFeatureSetNameArgument(cls, subparser):
        subparser.add_argument(
            "featureSetName",
            help="the name of the feature set")

    @classmethod
    def addContinuousSetNameArgument(cls, subparser):
        subparser.add_argument(
            "continuousSetName",
            help="the name of the continuous set")

    @classmethod
    def addIndividualNameArgument(cls, subparser):
        subparser.add_argument(
            "individualName",
            help="the name of the individual")

    @classmethod
    def addBiosampleNameArgument(cls, subparser):
        subparser.add_argument(
            "biosampleName",
            help="the name of the biosample")

    @classmethod
    def addBiosampleArgument(cls, subparser):
        subparser.add_argument(
            "biosample",
            help="the JSON of the biosample")

    @classmethod
    def addIndividualArgument(cls, subparser):
        subparser.add_argument(
            "individual",
            help="the JSON of the individual")

    @classmethod
    def addFilePathArgument(cls, subparser, helpText):
        subparser.add_argument("filePath", help=helpText)

    @classmethod
    def addDirPathArgument(cls, subparser, helpText):
        subparser.add_argument("dirPath", help=helpText)

    @classmethod
    def addNameOption(cls, parser, objectType):
        parser.add_argument(
            "-n", "--name", default=None,
            help="The name of the {}".format(objectType))

    @classmethod
    def addNameArgument(cls, parser, objectType):
        parser.add_argument(
            "name", help="The name of the {}".format(objectType))

    @classmethod
    def addRnaQuantificationNameArgument(cls, subparser):
        subparser.add_argument(
            "rnaQuantificationName",
            help="the name of the RNA Quantification")

    @classmethod
    def addClassNameOption(cls, subparser, objectType):
        helpText = (
            "the name of the class used to "
            "fetch features in this {}"
        ).format(objectType)
        subparser.add_argument(
            "-C", "--className",
            default="ga4gh.datamodel.sequence_annotations.Gff3DbFeatureSet",
            help=helpText)

    @classmethod
    def addRnaQuantificationSetNameArgument(cls, subparser):
        subparser.add_argument(
            "rnaQuantificationSetName",
            help="the name of the RNA Quantification Set")

    @classmethod
    def addQuantificationFilePathArgument(cls, subparser, helpText):
        subparser.add_argument("quantificationFilePath", help=helpText)

    @classmethod
    def addRnaFormatArgument(cls, subparser):
        subparser.add_argument(
            "format", help="format of the quantification input data")

    @classmethod
    def addRnaFeatureTypeOption(cls, subparser):
        subparser.add_argument(
            "-t", "--transcript", action="store_true", default=False,
            help="sets the quantification type to transcript")

    @classmethod
    def getParser(cls):
        parser = common_cli.createArgumentParser(
            "GA4GH data repository management tool")
        subparsers = parser.add_subparsers(title='subcommands',)
        cli.addVersionArgument(parser)

        initParser = common_cli.addSubparser(
            subparsers, "init", "Initialize a data repository")
        initParser.set_defaults(runner="init")
        cls.addRepoArgument(initParser)
        cls.addForceOption(initParser)

        verifyParser = common_cli.addSubparser(
            subparsers, "verify",
            "Verifies the repository by examing all data files")
        verifyParser.set_defaults(runner="verify")
        cls.addRepoArgument(verifyParser)

        listParser = common_cli.addSubparser(
            subparsers, "list", "List the contents of the repo")
        listParser.set_defaults(runner="list")
        cls.addRepoArgument(listParser)

        listAnnouncementsParser = common_cli.addSubparser(
            subparsers, "list-announcements", "List the announcements in"
                                              "the repo.")
        listAnnouncementsParser.set_defaults(runner="listAnnouncements")
        cls.addRepoArgument(listAnnouncementsParser)

        clearAnnouncementsParser = common_cli.addSubparser(
            subparsers, "clear-announcements", "List the announcements in"
                                               "the repo.")
        clearAnnouncementsParser.set_defaults(runner="clearAnnouncements")
        cls.addRepoArgument(clearAnnouncementsParser)

        addPeerParser = common_cli.addSubparser(
            subparsers, "add-peer", "Add a peer to the registry by URL.")
        addPeerParser.set_defaults(runner="addPeer")
        cls.addRepoArgument(addPeerParser)
        cls.addUrlArgument(addPeerParser)
        cls.addAttributesArgument(addPeerParser)

        removePeerParser = common_cli.addSubparser(
            subparsers, "remove-peer", "Remove a peer from "
                                       "the registry by URL.")
        removePeerParser.set_defaults(runner="removePeer")
        cls.addRepoArgument(removePeerParser)
        cls.addUrlArgument(removePeerParser)
        cls.addForceOption(removePeerParser)

        addDatasetParser = common_cli.addSubparser(
            subparsers, "add-dataset", "Add a dataset to the data repo")
        addDatasetParser.set_defaults(runner="addDataset")
        cls.addRepoArgument(addDatasetParser)
        cls.addDatasetNameArgument(addDatasetParser)
        cls.addAttributesArgument(addDatasetParser)
        cls.addDescriptionOption(addDatasetParser, "dataset")

        removeDatasetParser = common_cli.addSubparser(
            subparsers, "remove-dataset",
            "Remove a dataset from the data repo")
        removeDatasetParser.set_defaults(runner="removeDataset")
        cls.addRepoArgument(removeDatasetParser)
        cls.addDatasetNameArgument(removeDatasetParser)
        cls.addForceOption(removeDatasetParser)

        objectType = "reference set"
        addReferenceSetParser = common_cli.addSubparser(
            subparsers, "add-referenceset",
            "Add a reference set to the data repo")
        addReferenceSetParser.set_defaults(runner="addReferenceSet")
        cls.addRepoArgument(addReferenceSetParser)
        cls.addFilePathArgument(
            addReferenceSetParser,
            "The path of the FASTA file to use as a reference set. This "
            "file must be bgzipped and indexed.")
        cls.addAttributesArgument(addReferenceSetParser)
        cls.addRelativePathOption(addReferenceSetParser)
        cls.addNameOption(addReferenceSetParser, objectType)
        cls.addDescriptionOption(addReferenceSetParser, objectType)
        addReferenceSetParser.add_argument(
            "--species", default=None,
            help="The species ontology term as a JSON string")
        addReferenceSetParser.add_argument(
            "--isDerived", default=False, type=bool,
            help="Indicates if this reference set is derived from another")
        addReferenceSetParser.add_argument(
            "--assemblyId", default=None,
            help="The assembly id")
        addReferenceSetParser.add_argument(
            "--sourceAccessions", default=None,
            help="The source accessions (pass as comma-separated list)")
        addReferenceSetParser.add_argument(
            "--sourceUri", default=None,
            help="The source URI")

        removeReferenceSetParser = common_cli.addSubparser(
            subparsers, "remove-referenceset",
            "Remove a reference set from the repo")
        removeReferenceSetParser.set_defaults(runner="removeReferenceSet")
        cls.addRepoArgument(removeReferenceSetParser)
        removeReferenceSetParser.add_argument(
            "referenceSetName",
            help="the name of the reference set")
        cls.addForceOption(removeReferenceSetParser)

        objectType = "ReadGroupSet"
        addReadGroupSetParser = common_cli.addSubparser(
            subparsers, "add-readgroupset",
            "Add a read group set to the data repo")
        addReadGroupSetParser.set_defaults(runner="addReadGroupSet")
        cls.addRepoArgument(addReadGroupSetParser)
        cls.addDatasetNameArgument(addReadGroupSetParser)
        cls.addNameOption(addReadGroupSetParser, objectType)
        cls.addReferenceSetNameOption(addReadGroupSetParser, "ReadGroupSet")
        cls.addAttributesArgument(addReadGroupSetParser)
        cls.addRelativePathOption(addReadGroupSetParser)
        addReadGroupSetParser.add_argument(
            "dataFile",
            help="The file path or URL of the BAM file for this ReadGroupSet")
        addReadGroupSetParser.add_argument(
            "-I", "--indexFile", default=None,
            help=(
                "The file path of the BAM index for this ReadGroupSet. "
                "If the dataFile argument is a local file, this will "
                "be automatically inferred by appending '.bai' to the "
                "file name. If the dataFile is a remote URL the path to "
                "a local file containing the BAM index must be provided"))

        addOntologyParser = common_cli.addSubparser(
            subparsers, "add-ontology",
            "Adds an ontology in OBO format to the repo. Currently, "
            "a sequence ontology (SO) instance is required to translate "
            "ontology term names held in annotations to ontology IDs. "
            "Sequence ontology files can be found at "
            "https://github.com/The-Sequence-Ontology/SO-Ontologies")
        addOntologyParser.set_defaults(runner="addOntology")
        cls.addRepoArgument(addOntologyParser)
        cls.addFilePathArgument(
            addOntologyParser,
            "The path of the OBO file defining this ontology.")
        cls.addRelativePathOption(addOntologyParser)
        cls.addNameOption(addOntologyParser, "ontology")

        removeOntologyParser = common_cli.addSubparser(
            subparsers, "remove-ontology",
            "Remove an ontology from the repo")
        removeOntologyParser.set_defaults(runner="removeOntology")
        cls.addRepoArgument(removeOntologyParser)
        cls.addOntologyNameArgument(removeOntologyParser)
        cls.addForceOption(removeOntologyParser)

        removeReadGroupSetParser = common_cli.addSubparser(
            subparsers, "remove-readgroupset",
            "Remove a read group set from the repo")
        removeReadGroupSetParser.set_defaults(runner="removeReadGroupSet")
        cls.addRepoArgument(removeReadGroupSetParser)
        cls.addDatasetNameArgument(removeReadGroupSetParser)
        cls.addReadGroupSetNameArgument(removeReadGroupSetParser)
        cls.addForceOption(removeReadGroupSetParser)

        objectType = "VariantSet"
        addVariantSetParser = common_cli.addSubparser(
            subparsers, "add-variantset",
            "Add a variant set to the data repo based on one or "
            "more VCF files. ")
        addVariantSetParser.set_defaults(runner="addVariantSet")
        cls.addRepoArgument(addVariantSetParser)
        cls.addDatasetNameArgument(addVariantSetParser)
        cls.addRelativePathOption(addVariantSetParser)
        addVariantSetParser.add_argument(
            "dataFiles", nargs="+",
            help=(
                "The VCF/BCF files representing the new VariantSet. "
                "These may be specified either one or more paths "
                "to local files or remote URLS, or as a path to "
                "a local directory containing VCF files. Either "
                "a single directory argument may be passed or a "
                "list of file paths/URLS, but not a mixture of "
                "directories and paths.")
            )
        addVariantSetParser.add_argument(
            "-I", "--indexFiles", nargs="+", metavar="indexFiles",
            help=(
                "The index files for the VCF/BCF files provided in "
                "the dataFiles argument. These must be provided in the "
                "same order as the data files."))
        cls.addNameOption(addVariantSetParser, objectType)
        cls.addReferenceSetNameOption(addVariantSetParser, objectType)
        cls.addSequenceOntologyNameOption(addVariantSetParser, objectType)
        cls.addAttributesArgument(addVariantSetParser)
        addVariantSetParser.add_argument(
            "-a", "--addAnnotationSets", action="store_true",
            help=(
                "If the supplied VCF file contains annotations, create the "
                "corresponding VariantAnnotationSet."))

        removeVariantSetParser = common_cli.addSubparser(
            subparsers, "remove-variantset",
            "Remove a variant set from the repo")
        removeVariantSetParser.set_defaults(runner="removeVariantSet")
        cls.addRepoArgument(removeVariantSetParser)
        cls.addDatasetNameArgument(removeVariantSetParser)
        cls.addVariantSetNameArgument(removeVariantSetParser)
        cls.addForceOption(removeVariantSetParser)

        addFeatureSetParser = common_cli.addSubparser(
            subparsers, "add-featureset", "Add a feature set to the data repo")
        addFeatureSetParser.set_defaults(runner="addFeatureSet")
        cls.addRepoArgument(addFeatureSetParser)
        cls.addDatasetNameArgument(addFeatureSetParser)
        cls.addAttributesArgument(addFeatureSetParser)
        cls.addRelativePathOption(addFeatureSetParser)
        cls.addFilePathArgument(
            addFeatureSetParser,
            "The path to the converted SQLite database containing Feature "
            "data")
        cls.addReferenceSetNameOption(addFeatureSetParser, "feature set")
        cls.addSequenceOntologyNameOption(addFeatureSetParser, "feature set")
        cls.addClassNameOption(addFeatureSetParser, "feature set")

        removeFeatureSetParser = common_cli.addSubparser(
            subparsers, "remove-featureset",
            "Remove a feature set from the repo")
        removeFeatureSetParser.set_defaults(runner="removeFeatureSet")
        cls.addRepoArgument(removeFeatureSetParser)
        cls.addDatasetNameArgument(removeFeatureSetParser)
        cls.addFeatureSetNameArgument(removeFeatureSetParser)
        cls.addForceOption(removeFeatureSetParser)

        addContinuousSetParser = common_cli.addSubparser(
            subparsers, "add-continuousset",
            "Add a continuous set to the data repo")
        addContinuousSetParser.set_defaults(runner="addContinuousSet")
        cls.addRepoArgument(addContinuousSetParser)
        cls.addDatasetNameArgument(addContinuousSetParser)
        cls.addRelativePathOption(addContinuousSetParser)
        cls.addFilePathArgument(
            addContinuousSetParser,
            "The path to the file contianing the continuous data ")
        cls.addReferenceSetNameOption(addContinuousSetParser, "continuous set")
        cls.addClassNameOption(addContinuousSetParser, "continuous set")

        removeContinuousSetParser = common_cli.addSubparser(
            subparsers, "remove-continuousset",
            "Remove a continuous set from the repo")
        removeContinuousSetParser.set_defaults(runner="removeContinuousSet")
        cls.addRepoArgument(removeContinuousSetParser)
        cls.addDatasetNameArgument(removeContinuousSetParser)
        cls.addContinuousSetNameArgument(removeContinuousSetParser)
        cls.addForceOption(removeContinuousSetParser)

        addBiosampleParser = common_cli.addSubparser(
            subparsers, "add-biosample", "Add a Biosample to the dataset")
        addBiosampleParser.set_defaults(runner="addBiosample")
        cls.addRepoArgument(addBiosampleParser)
        cls.addDatasetNameArgument(addBiosampleParser)
        cls.addBiosampleNameArgument(addBiosampleParser)
        cls.addBiosampleArgument(addBiosampleParser)

        removeBiosampleParser = common_cli.addSubparser(
            subparsers, "remove-biosample",
            "Remove a Biosample from the repo")
        removeBiosampleParser.set_defaults(runner="removeBiosample")
        cls.addRepoArgument(removeBiosampleParser)
        cls.addDatasetNameArgument(removeBiosampleParser)
        cls.addBiosampleNameArgument(removeBiosampleParser)
        cls.addForceOption(removeBiosampleParser)

        addIndividualParser = common_cli.addSubparser(
            subparsers, "add-individual", "Add an Individual to the dataset")
        addIndividualParser.set_defaults(runner="addIndividual")
        cls.addRepoArgument(addIndividualParser)
        cls.addDatasetNameArgument(addIndividualParser)
        cls.addIndividualNameArgument(addIndividualParser)
        cls.addIndividualArgument(addIndividualParser)

        removeIndividualParser = common_cli.addSubparser(
            subparsers, "remove-individual",
            "Remove an Individual from the repo")
        removeIndividualParser.set_defaults(runner="removeIndividual")
        cls.addRepoArgument(removeIndividualParser)
        cls.addDatasetNameArgument(removeIndividualParser)
        cls.addIndividualNameArgument(removeIndividualParser)
        cls.addForceOption(removeIndividualParser)

        objectType = "RnaQuantification"
        addRnaQuantificationParser = common_cli.addSubparser(
            subparsers, "add-rnaquantification",
            "Add an RNA quantification to the data repo")
        addRnaQuantificationParser.set_defaults(
            runner="addRnaQuantification")
        cls.addFilePathArgument(
            addRnaQuantificationParser,
            "The path to the RNA SQLite database to create or modify")
        cls.addQuantificationFilePathArgument(
            addRnaQuantificationParser, "The path to the expression file.")
        cls.addRnaFormatArgument(addRnaQuantificationParser)
        cls.addRepoArgument(addRnaQuantificationParser)
        cls.addDatasetNameArgument(addRnaQuantificationParser)
        addRnaQuantificationParser.add_argument(
            "--biosampleName", default=None, help="Biosample Name")
        addRnaQuantificationParser.add_argument(
            "--readGroupSetName", default=None, help="Read Group Set Name")
        addRnaQuantificationParser.add_argument(
            "--featureSetNames", default=None, help="Comma separated list")
        cls.addNameOption(addRnaQuantificationParser, "rna quantification")
        cls.addDescriptionOption(addRnaQuantificationParser, objectType)
        cls.addRnaFeatureTypeOption(addRnaQuantificationParser)
        cls.addAttributesArgument(addRnaQuantificationParser)

        objectType = "RnaQuantificationSet"
        initRnaQuantificationSetParser = common_cli.addSubparser(
            subparsers, "init-rnaquantificationset",
            "Initializes an RNA quantification set")
        initRnaQuantificationSetParser.set_defaults(
            runner="initRnaQuantificationSet")
        cls.addRepoArgument(initRnaQuantificationSetParser)
        cls.addFilePathArgument(
            initRnaQuantificationSetParser,
            "The path to the resulting Quantification Set")

        addRnaQuantificationSetParser = common_cli.addSubparser(
            subparsers, "add-rnaquantificationset",
            "Add an RNA quantification set to the data repo")
        addRnaQuantificationSetParser.set_defaults(
            runner="addRnaQuantificationSet")
        cls.addRepoArgument(addRnaQuantificationSetParser)
        cls.addDatasetNameArgument(addRnaQuantificationSetParser)
        cls.addFilePathArgument(
            addRnaQuantificationSetParser,
            "The path to the converted SQLite database containing RNA data")
        cls.addReferenceSetNameOption(
            addRnaQuantificationSetParser, objectType)
        cls.addNameOption(addRnaQuantificationSetParser, objectType)
        cls.addAttributesArgument(addRnaQuantificationSetParser)

        removeRnaQuantificationSetParser = common_cli.addSubparser(
            subparsers, "remove-rnaquantificationset",
            "Remove an RNA quantification set from the repo")
        removeRnaQuantificationSetParser.set_defaults(
            runner="removeRnaQuantificationSet")
        cls.addRepoArgument(removeRnaQuantificationSetParser)
        cls.addDatasetNameArgument(removeRnaQuantificationSetParser)
        cls.addRnaQuantificationSetNameArgument(
            removeRnaQuantificationSetParser)
        cls.addForceOption(removeRnaQuantificationSetParser)

        addPhenotypeAssociationSetParser = common_cli.addSubparser(
            subparsers, "add-phenotypeassociationset",
            "Adds phenotypes in ttl format to the repo.")
        addPhenotypeAssociationSetParser.set_defaults(
            runner="addPhenotypeAssociationSet")
        cls.addRepoArgument(addPhenotypeAssociationSetParser)
        cls.addDatasetNameArgument(addPhenotypeAssociationSetParser)
        cls.addDirPathArgument(
            addPhenotypeAssociationSetParser,
            "The path of the directory containing ttl files.")
        cls.addNameOption(
            addPhenotypeAssociationSetParser,
            "PhenotypeAssociationSet")
        cls.addAttributesArgument(addPhenotypeAssociationSetParser)

        removePhenotypeAssociationSetParser = common_cli.addSubparser(
            subparsers, "remove-phenotypeassociationset",
            "Remove an phenotypes from the repo")
        removePhenotypeAssociationSetParser.set_defaults(
            runner="removePhenotypeAssociationSet")
        cls.addRepoArgument(removePhenotypeAssociationSetParser)
        cls.addDatasetNameArgument(removePhenotypeAssociationSetParser)
        cls.addNameArgument(
            removePhenotypeAssociationSetParser,
            "phenotype association set")
        cls.addForceOption(removePhenotypeAssociationSetParser)

        return parser

    @classmethod
    def runCommand(cls, args):
        parser = cls.getParser()
        parsedArgs = parser.parse_args(args)
        if "runner" not in parsedArgs:
            parser.print_help()
        manager = RepoManager(parsedArgs)
        runMethod = getattr(manager, parsedArgs.runner)
        runMethod()


def getRepoManagerParser():
    """
    Used by sphinx.argparse.
    """
    return RepoManager.getParser()


def repoExitError(message):
    """
    Exits the repo manager with error status.
    """
    wrapper = textwrap.TextWrapper(
        break_on_hyphens=False, break_long_words=False)
    formatted = wrapper.fill("{}: error: {}".format(sys.argv[0], message))
    sys.exit(formatted)


def repo_main(args=None):
    try:
        RepoManager.runCommand(args)
    except exceptions.RepoManagerException as exception:
        # These are exceptions that happen throughout the CLI, and are
        # used to communicate back to the user
        repoExitError(str(exception))
    except exceptions.NotFoundException as exception:
        # We expect NotFoundExceptions to occur when we're looking for
        # datasets, readGroupsets, etc.
        repoExitError(str(exception))
    except exceptions.DataException as exception:
        # We expect DataExceptions to occur when a file open fails,
        # a URL cannot be reached, etc.
        repoExitError(str(exception))
    except Exception as exception:
        # Uncaught exception: this is a bug
        message = """
An internal error has occurred.  Please file a bug report at
https://github.com/ga4gh/ga4gh-server/issues
with all the relevant details, and the following stack trace.
"""
        print("{}: error:".format(sys.argv[0]), file=sys.stderr)
        print(message, file=sys.stderr)
        traceback.print_exception(*sys.exc_info())
        sys.exit(1)
