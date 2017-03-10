"""
The backing data store for the GA4GH server
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import os
import datetime

import ga4gh.server.datamodel as datamodel
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.datamodel.ontologies as ontologies
import ga4gh.server.datamodel.reads as reads
import ga4gh.server.datamodel.references as references
import ga4gh.server.datamodel.variants as variants
import ga4gh.server.datamodel.sequence_annotations as sequence_annotations
import ga4gh.server.datamodel.continuous as continuous
import ga4gh.server.datamodel.bio_metadata as biodata
import ga4gh.server.datamodel.genotype_phenotype as genotype_phenotype
import ga4gh.server.datamodel.genotype_phenotype_featureset as g2pFeatureset
import ga4gh.server.datamodel.rna_quantification as rna_quantification
import ga4gh.server.datamodel.peers as peers
import ga4gh.server.exceptions as exceptions
import ga4gh.server.repo.models as models

import ga4gh.schemas.protocol as protocol

MODE_READ = 'r'
MODE_WRITE = 'w'


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
        self._ontologyIdMap = {}
        self._ontologyIds = []
        self._peers = []

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
        Add an ontology map to this data repository.
        """
        self._ontologyNameMap[ontology.getName()] = ontology
        self._ontologyIdMap[ontology.getId()] = ontology
        self._ontologyIds.append(ontology.getId())

    def getDatasets(self):
        """
        Returns a list of datasets in this data repository
        """
        return [self._datasetIdMap[id_] for id_ in self._datasetIds]

    def insertPeer(self, peer):
        """
        Adds a peer to the list of peers in the repository. Used only in
        testing.
        """
        self._peers.append(peer)

    def getPeer(self, url):
        """
        Select the first peer in the datarepo with the given url simulating
        the behavior of selecting by URL. This is only used during testing.
        """
        peers = filter(lambda x: x.getUrl() == url, self.getPeers())
        if len(peers) == 0:
            raise exceptions.PeerNotFoundException(url)
        return peers[0]

    def getPeers(self, offset=0, limit=100):
        """
        Returns the list of peers with an optional offset and limit
        simulating the SQL registry for testing.
        """
        return self._peers[offset:offset + limit]

    def insertAnnouncement(self, announcement):
        """
        A placeholder function to simulate receiving an announcement used
        in testing. It will throw an exception if the URL is invalid.
        """
        peers.Peer(announcement.get('url'))

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

    def getOntology(self, id_):
        """
        Returns the ontology with the specified ID.
        """
        if id_ not in self._ontologyIdMap:
            raise exceptions.OntologyNotFoundException(id_)
        return self._ontologyIdMap[id_]

    def getOntologyByName(self, name):
        """
        Returns an ontology by name
        """
        if name not in self._ontologyNameMap:
            raise exceptions.OntologyNameNotFoundException(name)
        return self._ontologyNameMap[name]

    def getOntologys(self):
        """
        Returns all ontologys in the repo
        """
        return [self._ontologyIdMap[id_] for id_ in self._ontologyIds]

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

    def getReadGroupSet(self, id_):
        """
        Returns the readgroup set with the specified ID.
        """
        compoundId = datamodel.ReadGroupSetCompoundId.parse(id_)
        dataset = self.getDataset(compoundId.dataset_id)
        return dataset.getReadGroupSet(id_)

    def getVariantSet(self, id_):
        """
        Returns the readgroup set with the specified ID.
        """
        compoundId = datamodel.VariantSetCompoundId.parse(id_)
        dataset = self.getDataset(compoundId.dataset_id)
        return dataset.getVariantSet(id_)

    def printSummary(self):
        """
        Prints a summary of this data repository to stdout.
        """
        print("Ontologies:")
        for ontology in self.getOntologys():
            print(
                "",
                ontology.getOntologyPrefix(),
                ontology.getName(),
                ontology.getDataUrl(),
                sep="\t")
        print("ReferenceSets:")
        for referenceSet in self.getReferenceSets():
            print(
                "", referenceSet.getLocalId(), referenceSet.getId(),
                referenceSet.getDescription(), referenceSet.getDataUrl(),
                sep="\t")
            for reference in referenceSet.getReferences():
                print(
                    "\t", reference.getLocalId(), reference.getId(),
                    sep="\t")
        print("Datasets:")
        for dataset in self.getDatasets():
            print(
                "", dataset.getLocalId(), dataset.getId(),
                dataset.getDescription(), sep="\t")
            print("\tReadGroupSets:")
            for readGroupSet in dataset.getReadGroupSets():
                print(
                    "\t", readGroupSet.getLocalId(),
                    readGroupSet.getReferenceSet().getLocalId(),
                    readGroupSet.getId(),
                    readGroupSet.getDataUrl(), sep="\t")
                for readGroup in readGroupSet.getReadGroups():
                    print(
                        "\t\t", readGroup.getId(), readGroup.getLocalId(),
                        sep="\t")
            print("\tVariantSets:")
            for variantSet in dataset.getVariantSets():
                print(
                    "\t", variantSet.getLocalId(),
                    variantSet.getReferenceSet().getLocalId(),
                    variantSet.getId(),
                    sep="\t")
                if variantSet.getNumVariantAnnotationSets() > 0:
                    print("\t\tVariantAnnotationSets:")
                    for vas in variantSet.getVariantAnnotationSets():
                        print(
                            "\t\t", vas.getLocalId(),
                            vas.getAnnotationType(),
                            vas.getOntology().getName(), sep="\t")
            print("\tFeatureSets:")
            for featureSet in dataset.getFeatureSets():
                print(
                    "\t", featureSet.getLocalId(),
                    featureSet.getReferenceSet().getLocalId(),
                    featureSet.getOntology().getName(),
                    featureSet.getId(),
                    sep="\t")
            print("\tContinuousSets:")
            for continuousSet in dataset.getContinuousSets():
                print(
                    "\t", continuousSet.getLocalId(),
                    continuousSet.getReferenceSet().getLocalId(),
                    continuousSet.getId(),
                    sep="\t")
            print("\tPhenotypeAssociationSets:")
            for phenotypeAssociationSet in \
                    dataset.getPhenotypeAssociationSets():
                print(
                    "\t", phenotypeAssociationSet.getLocalId(),
                    phenotypeAssociationSet.getParentContainer().getId(),
                    sep="\t")
                # TODO -  please improve this listing
            print("\tRnaQuantificationSets:")
            for rna_quantification_set in dataset.getRnaQuantificationSets():
                print(
                    "\t", rna_quantification_set.getLocalId(),
                    rna_quantification_set.getId(), sep="\t")
                for quant in rna_quantification_set.getRnaQuantifications():
                        print(
                            "\t\t", quant.getLocalId(),
                            quant._description,
                            ",".join(quant._readGroupIds),
                            ",".join(quant._featureSetIds), sep="\t")

    def allReferences(self):
        """
        Return an iterator over all references in the data repo
        """
        for referenceSet in self.getReferenceSets():
            for reference in referenceSet.getReferences():
                yield reference

    def allBiosamples(self):
        """
        Return an iterator over all biosamples in the data repo
        """
        for dataset in self.getDatasets():
            for biosample in dataset.getBiosamples():
                yield biosample

    def allIndividuals(self):
        """
        Return an iterator over all individuals in the data repo
        """
        for dataset in self.getDatasets():
            for individual in dataset.getIndividuals():
                yield individual

    def allReadGroupSets(self):
        """
        Return an iterator over all read group sets in the data repo
        """
        for dataset in self.getDatasets():
            for readGroupSet in dataset.getReadGroupSets():
                yield readGroupSet

    def allReadGroups(self):
        """
        Return an iterator over all read groups in the data repo
        """
        for dataset in self.getDatasets():
            for readGroupSet in dataset.getReadGroupSets():
                for readGroup in readGroupSet.getReadGroups():
                    yield readGroup

    def allVariantSets(self):
        """
        Return an iterator over all read variant sets in the data repo
        """
        for dataset in self.getDatasets():
            for variantSet in dataset.getVariantSets():
                yield variantSet

    def allFeatureSets(self):
        """
        Return an iterator over all feature sets in the data repo
        """
        for dataset in self.getDatasets():
            for featureSet in dataset.getFeatureSets():
                yield featureSet

    def allFeatures(self):
        """
        Return an iterator over all features in the data repo
        """
        for dataset in self.getDatasets():
            for featureSet in dataset.getFeatureSets():
                for feature in featureSet.getFeatures():
                    yield feature

    def allContinuousSets(self):
        """
        Return an iterator over all continuous sets in the data repo
        """
        for dataset in self.getDatasets():
            for continuousSet in dataset.getContinuousSets():
                yield continuousSet

    def allCallSets(self):
        """
        Return an iterator over all call sets in the data repo
        """
        for dataset in self.getDatasets():
            for variantSet in dataset.getVariantSets():
                for callSet in variantSet.getCallSets():
                    yield callSet

    def allVariantAnnotationSets(self):
        """
        Return an iterator over all variant annotation sets
        in the data repo
        """
        for dataset in self.getDatasets():
            for variantSet in dataset.getVariantSets():
                for vaSet in variantSet.getVariantAnnotationSets():
                    yield vaSet

    def allPhenotypeAssociationSets(self):
        """
        Return an iterator over all phenotype association sets
        in the data repo
        """
        for dataset in self.getDatasets():
            for paSet in dataset.getPhenotypeAssociationSets():
                yield paSet

    def allRnaQuantificationSets(self):
        """
        Return an iterator over all rna quantification sets
        """
        for dataset in self.getDatasets():
            for rnaQuantificationSet in dataset.getRnaQuantificationSets():
                yield rnaQuantificationSet

    def allRnaQuantifications(self):
        """
        Return an iterator over all rna quantifications
        """
        for dataset in self.getDatasets():
            for rnaQuantificationSet in dataset.getRnaQuantificationSets():
                for rnaQuantification in \
                        rnaQuantificationSet.getRnaQuantifications():
                    yield rnaQuantification

    def allExpressionLevels(self):
        """
        Return an iterator over all expression levels
        """
        for dataset in self.getDatasets():
            for rnaQuantificationSet in dataset.getRnaQuantificationSets():
                for rnaQuantification in \
                        rnaQuantificationSet.getRnaQuantifications():
                    for expressionLevel in \
                            rnaQuantification.getExpressionLevels():
                        yield expressionLevel


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
            numVariantSets=1, numCalls=1, variantDensity=0.5,
            numReferenceSets=1, numReferencesPerReferenceSet=1,
            numReadGroupSets=1, numReadGroupsPerReadGroupSet=1,
            numPhenotypeAssociations=2,
            numPhenotypeAssociationSets=1,
            numAlignments=2, numRnaQuantSets=2, numExpressionLevels=2,
            numPeers=200):
        super(SimulatedDataRepository, self).__init__()

        for i in xrange(numPeers):
            peer = peers.Peer("http://test{}.org".format(i))
            self.insertPeer(peer)

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
                numAlignments=numAlignments,
                numPhenotypeAssociations=numPhenotypeAssociations,
                numPhenotypeAssociationSets=numPhenotypeAssociationSets,
                numRnaQuantSets=numRnaQuantSets,
                numExpressionLevels=numExpressionLevels)
            self.addDataset(dataset)


class SqlDataRepository(AbstractDataRepository):
    """
    A data repository based on a SQL database.
    """
    class SchemaVersion(object):
        """
        The version of the data repository SQL schema
        """
        def __init__(self, versionString):
            splits = versionString.split('.')
            assert len(splits) == 2
            self.major = splits[0]
            self.minor = splits[1]

        def __str__(self):
            return "{}.{}".format(self.major, self.minor)

    version = SchemaVersion("2.1")
    systemKeySchemaVersion = "schemaVersion"
    systemKeyCreationTimeStamp = "creationTimeStamp"

    def __init__(self, fileName):
        super(SqlDataRepository, self).__init__()
        self._dbFilename = fileName
        # We open the repo in either read or write mode. When we want to
        # update the repo we open it in write mode. For normal online
        # server use, we open it in read mode.
        self._openMode = None
        # Values filled in using the DB. These will all be None until
        # we have called load()
        self._schemaVersion = None
        # Connection to the DB.
        self.database = models.SqliteDatabase(self._dbFilename, **{})
        models.databaseProxy.initialize(self.database)

    def _checkWriteMode(self):
        if self._openMode != MODE_WRITE:
            raise ValueError("Repo must be opened in write mode")

    def getPeer(self, url):
        """
        Finds a peer by URL and return the first peer record with that URL.
        """
        peers = list(models.Peer.select().where(models.Peer.url == url))
        if len(peers) == 0:
            raise exceptions.PeerNotFoundException(url)
        return peers[0]

    def getPeers(self, offset=0, limit=1000):
        """
        Get the list of peers using an SQL offset and limit. Returns a list
        of peer datamodel objects in a list.
        """
        select = models.Peer.select().order_by(
            models.Peer.url).limit(limit).offset(offset)
        return [peers.Peer(p.url, record=p) for p in select]

    def tableToTsv(self, model):
        """
        Takes a model class and attempts to create a table in TSV format
        that can be imported into a spreadsheet program.
        """
        first = True
        for item in model.select():
            if first:
                header = "".join(
                    ["{}\t".format(x) for x in model._meta.fields.keys()])
                print(header)
                first = False
            row = "".join(
                ["{}\t".format(
                    getattr(item, key)) for key in model._meta.fields.keys()])
            print(row)

    def printAnnouncements(self):
        """
        Prints the announcement table to the log in tsv format.
        """
        self.tableToTsv(models.Announcement)

    def clearAnnouncements(self):
        """
        Flushes the announcement table.
        """
        try:
            q = models.Announcement.delete().where(
                models.Announcement.id > 0)
            q.execute()
        except Exception as e:
            raise exceptions.RepoManagerException(e)

    def insertAnnouncement(self, announcement):
        """
        Adds an announcement to the registry for later analysis.
        """
        url = announcement.get('url', None)
        try:
            peers.Peer(url)
        except:
            raise exceptions.BadUrlException(url)
        try:
            # TODO get more details about the user agent
            models.Announcement.create(
                url=announcement.get('url'),
                attributes=json.dumps(announcement.get('attributes', {})),
                remote_addr=announcement.get('remote_addr', None),
                user_agent=announcement.get('user_agent', None))
        except Exception as e:
            raise exceptions.RepoManagerException(e)

    def open(self, mode=MODE_READ):
        """
        Opens this repo in the specified mode.

        TODO: figure out the correct semantics of this and document
        the intended future behaviour as well as the current
        transitional behaviour.
        """
        if mode not in [MODE_READ, MODE_WRITE]:
            error = "Open mode must be '{}' or '{}'".format(
                MODE_READ, MODE_WRITE)
            raise ValueError(error)
        self._openMode = mode
        if mode == MODE_READ:
            self.assertExists()
        if mode == MODE_READ:
            # This is part of the transitional behaviour where
            # we load the whole DB into memory to get access to
            # the data model.
            self.load()

    def commit(self):
        """
        Commits any changes made to the repo. It is an error to call
        this function if the repo is not opened in write-mode.
        """
        self._checkWriteMode()

    def close(self):
        """
        Closes this repo.
        """
        if self._openMode is None:
            raise ValueError("Repo already closed")
        self._openMode = None

    def verify(self):
        """
        Verifies that the data in the repository is consistent.
        """
        # TODO this should emit to a log that we can configure so we can
        # have verbosity levels. We should provide a way to configure
        # where we look at various chromosomes and so on. This will be
        # an important debug tool for administrators.
        for ontology in self.getOntologys():
            print(
                "Verifying Ontology", ontology.getName(),
                "@", ontology.getDataUrl())
            # TODO how do we verify this? Check some well-know SO terms?
        for referenceSet in self.getReferenceSets():
            print(
                "Verifying ReferenceSet", referenceSet.getLocalId(),
                "@", referenceSet.getDataUrl())
            for reference in referenceSet.getReferences():
                length = min(reference.getLength(), 1000)
                bases = reference.getBases(0, length)
                assert len(bases) == length
                print(
                    "\tReading", length, "bases from",
                    reference.getLocalId())
        for dataset in self.getDatasets():
            print("Verifying Dataset", dataset.getLocalId())
            for featureSet in dataset.getFeatureSets():
                for referenceSet in self.getReferenceSets():
                    # TODO cycle through references?
                    reference = referenceSet.getReferences()[0]
                    print(
                        "\tVerifying FeatureSet",
                        featureSet.getLocalId(),
                        "with reference", reference.getLocalId())
                    length = min(reference.getLength(), 1000)
                    features = featureSet.getFeatures(
                        reference.getLocalId(), 0, length, None, 3)
                    for feature in features:
                        print("\t{}".format(feature))
            # for continuousSet in dataset.getContinuousSets():
            # -- there is no getContinuous
            for readGroupSet in dataset.getReadGroupSets():
                print(
                    "\tVerifying ReadGroupSet", readGroupSet.getLocalId(),
                    "@", readGroupSet.getDataUrl())
                references = readGroupSet.getReferenceSet().getReferences()
                # TODO should we cycle through the references? Should probably
                # be an option.
                reference = references[0]
                max_alignments = 10
                for readGroup in readGroupSet.getReadGroups():
                    alignments = readGroup.getReadAlignments(reference)
                    for i, alignment in enumerate(alignments):
                        if i == max_alignments:
                            break
                    print(
                        "\t\tRead", i, "alignments from",
                        readGroup.getLocalId())
            for variantSet in dataset.getVariantSets():
                print("\tVerifying VariantSet", variantSet.getLocalId())
                max_variants = 10
                max_annotations = 10
                refMap = variantSet.getReferenceToDataUrlIndexMap()
                for referenceName, (dataUrl, indexFile) in refMap.items():
                    variants = variantSet.getVariants(referenceName, 0, 2**31)
                    for i, variant in enumerate(variants):
                        if i == max_variants:
                            break
                    print(
                        "\t\tRead", i, "variants from reference",
                        referenceName, "@", dataUrl)
                for annotationSet in variantSet.getVariantAnnotationSets():
                    print(
                        "\t\tVerifying VariantAnnotationSet",
                        annotationSet.getLocalId())
                    for referenceName in refMap.keys():
                        annotations = annotationSet.getVariantAnnotations(
                            referenceName, 0, 2**31)
                        for i, annotation in enumerate(annotations):
                            if i == max_annotations:
                                break
                    print(
                        "\t\t\tRead", i, "annotations from reference",
                        referenceName)
            for phenotypeAssociationSet \
                    in dataset.getPhenotypeAssociationSets():
                print("\t\tVerifying PhenotypeAssociationSet")
                print(
                    "\t\t\t", phenotypeAssociationSet.getLocalId(),
                    phenotypeAssociationSet.getParentContainer().getId(),
                    sep="\t")
                # TODO - please improve this verification,
                #        print out number of tuples in graph

    def _createSystemTable(self):
        self.database.create_table(models.System)
        models.System.create(
            key=self.systemKeySchemaVersion, value=self.version)
        models.System.create(
            key=self.systemKeyCreationTimeStamp, value=datetime.datetime.now())

    def _readSystemTable(self):
        if not self.exists():
            raise exceptions.RepoNotFoundException(
                self._dbFilename)
        try:
            self._schemaVersion = models.System.get(
                models.System.key == self.systemKeySchemaVersion).value
        except Exception:
            raise exceptions.RepoInvalidDatabaseException(self._dbFilename)
        schemaVersion = self.SchemaVersion(self._schemaVersion)
        if schemaVersion.major != self.version.major:
            raise exceptions.RepoSchemaVersionMismatchException(
                schemaVersion, self.version)

    def _createOntologyTable(self):
        self.database.create_table(models.Ontology)

    def insertOntology(self, ontology):
        """
        Inserts the specified ontology into this repository.
        """
        try:
            models.Ontology.create(
                    id=ontology.getName(),
                    name=ontology.getName(),
                    dataurl=ontology.getDataUrl(),
                    ontologyprefix=ontology.getOntologyPrefix())
        except Exception:
            raise exceptions.DuplicateNameException(
                ontology.getName())
        # TODO we need to create a proper ID when we're doing ID generation
        # for the rest of the container objects.

    def _readOntologyTable(self):
        for ont in models.Ontology.select():
            ontology = ontologies.Ontology(ont.name)
            ontology.populateFromRow(ont)
            self.addOntology(ontology)

    def removeOntology(self, ontology):
        """
        Removes the specified ontology term map from this repository.
        """
        q = models.Ontology.delete().where(id == ontology.getId())
        q.execute()

    def _createReferenceTable(self):
        self.database.create_table(models.Reference)

    def insertReference(self, reference):
        """
        Inserts the specified reference into this repository.
        """
        models.Reference.create(
            id=reference.getId(),
            referencesetid=reference.getParentContainer().getId(),
            name=reference.getLocalId(),
            length=reference.getLength(),
            isderived=reference.getIsDerived(),
            species=json.dumps(reference.getSpecies()),
            md5checksum=reference.getMd5Checksum(),
            sourceaccessions=json.dumps(reference.getSourceAccessions()),
            sourceuri=reference.getSourceUri())

    def _readReferenceTable(self):
        for referenceRecord in models.Reference.select():
            referenceSet = self.getReferenceSet(
                referenceRecord.referencesetid.id)
            reference = references.HtslibReference(
                referenceSet, referenceRecord.name)
            reference.populateFromRow(referenceRecord)
            assert reference.getId() == referenceRecord.id
            referenceSet.addReference(reference)

    def _createReferenceSetTable(self):
        self.database.create_table(models.Referenceset)

    def insertReferenceSet(self, referenceSet):
        """
        Inserts the specified referenceSet into this repository.
        """
        try:
            models.Referenceset.create(
                id=referenceSet.getId(),
                name=referenceSet.getLocalId(),
                description=referenceSet.getDescription(),
                assemblyid=referenceSet.getAssemblyId(),
                isderived=referenceSet.getIsDerived(),
                species=json.dumps(referenceSet.getSpecies()),
                md5checksum=referenceSet.getMd5Checksum(),
                sourceaccessions=json.dumps(
                    referenceSet.getSourceAccessions()),
                sourceuri=referenceSet.getSourceUri(),
                dataurl=referenceSet.getDataUrl())
            for reference in referenceSet.getReferences():
                self.insertReference(reference)
        except Exception:
            raise exceptions.DuplicateNameException(
                referenceSet.getLocalId())

    def _readReferenceSetTable(self):
        for referenceSetRecord in models.Referenceset.select():
            referenceSet = references.HtslibReferenceSet(
                referenceSetRecord.name)
            referenceSet.populateFromRow(referenceSetRecord)
            assert referenceSet.getId() == referenceSetRecord.id
            # Insert the referenceSet into the memory-based object model.
            self.addReferenceSet(referenceSet)

    def _createDatasetTable(self):
        self.database.create_table(models.Dataset)

    def insertDataset(self, dataset):
        """
        Inserts the specified dataset into this repository.
        """
        try:
            models.Dataset.create(
                id=dataset.getId(),
                name=dataset.getLocalId(),
                description=dataset.getDescription(),
                attributes=json.dumps(dataset.getAttributes()))
        except Exception:
            raise exceptions.DuplicateNameException(
                dataset.getLocalId())

    def removeDataset(self, dataset):
        """
        Removes the specified dataset from this repository. This performs
        a cascading removal of all items within this dataset.
        """
        for datasetRecord in models.Dataset.select().where(
                        models.Dataset.id == dataset.getId()):
            datasetRecord.delete_instance(recursive=True)

    def removePhenotypeAssociationSet(self, phenotypeAssociationSet):
        """
        Remove a phenotype association set from the repo
        """
        q = models.Phenotypeassociationset.delete().where(
            models.Phenotypeassociationset.id ==
            phenotypeAssociationSet.getId())
        q.execute()

    def removeFeatureSet(self, featureSet):
        """
        Removes the specified featureSet from this repository.
        """
        q = models.Featureset.delete().where(
            models.Featureset.id == featureSet.getId())
        q.execute()

    def removeContinuousSet(self, continuousSet):
        """
        Removes the specified continuousSet from this repository.
        """
        q = models.ContinuousSet.delete().where(
            models.ContinuousSet.id == continuousSet.getId())
        q.execute()

    def _readDatasetTable(self):
        for datasetRecord in models.Dataset.select():
            dataset = datasets.Dataset(datasetRecord.name)
            dataset.populateFromRow(datasetRecord)
            assert dataset.getId() == datasetRecord.id
            # Insert the dataset into the memory-based object model.
            self.addDataset(dataset)

    def _createReadGroupTable(self):
        self.database.create_table(models.Readgroup)

    def insertReadGroup(self, readGroup):
        """
        Inserts the specified readGroup into the DB.
        """
        statsJson = json.dumps(protocol.toJsonDict(readGroup.getStats()))
        experimentJson = json.dumps(
            protocol.toJsonDict(readGroup.getExperiment()))
        try:
            models.Readgroup.create(
                id=readGroup.getId(),
                readgroupsetid=readGroup.getParentContainer().getId(),
                name=readGroup.getLocalId(),
                predictedinsertedsize=readGroup.getPredictedInsertSize(),
                samplename=readGroup.getSampleName(),
                description=readGroup.getDescription(),
                stats=statsJson,
                experiment=experimentJson,
                biosampleid=readGroup.getBiosampleId(),
                attributes=json.dumps(readGroup.getAttributes()))
        except Exception as e:
            raise exceptions.RepoManagerException(e)

    def removeReadGroupSet(self, readGroupSet):
        """
        Removes the specified readGroupSet from this repository. This performs
        a cascading removal of all items within this readGroupSet.
        """
        for readGroupSetRecord in models.Readgroupset.select().where(
                        models.Readgroupset.id == readGroupSet.getId()):
            readGroupSetRecord.delete_instance(recursive=True)

    def removeVariantSet(self, variantSet):
        """
        Removes the specified variantSet from this repository. This performs
        a cascading removal of all items within this variantSet.
        """
        for variantSetRecord in models.Variantset.select().where(
                        models.Variantset.id == variantSet.getId()):
            variantSetRecord.delete_instance(recursive=True)

    def removeBiosample(self, biosample):
        """
        Removes the specified biosample from this repository.
        """
        q = models.Biosample.delete().where(
            models.Biosample.id == biosample.getId())
        q.execute()

    def removeIndividual(self, individual):
        """
        Removes the specified individual from this repository.
        """
        q = models.Individual.delete().where(
            models.Individual.id == individual.getId())
        q.execute()

    def _readReadGroupTable(self):
        for readGroupRecord in models.Readgroup.select():
            readGroupSet = self.getReadGroupSet(
                readGroupRecord.readgroupsetid.id)
            readGroup = reads.HtslibReadGroup(
                readGroupSet, readGroupRecord.name)
            # TODO set the reference set.
            readGroup.populateFromRow(readGroupRecord)
            assert readGroup.getId() == readGroupRecord.id
            # Insert the readGroupSet into the memory-based object model.
            readGroupSet.addReadGroup(readGroup)

    def _createReadGroupSetTable(self):
        self.database.create_table(models.Readgroupset)

    def insertReadGroupSet(self, readGroupSet):
        """
        Inserts a the specified readGroupSet into this repository.
        """
        programsJson = json.dumps(
            [protocol.toJsonDict(program) for program in
             readGroupSet.getPrograms()])
        statsJson = json.dumps(protocol.toJsonDict(readGroupSet.getStats()))
        try:
            models.Readgroupset.create(
                id=readGroupSet.getId(),
                datasetid=readGroupSet.getParentContainer().getId(),
                referencesetid=readGroupSet.getReferenceSet().getId(),
                name=readGroupSet.getLocalId(),
                programs=programsJson,
                stats=statsJson,
                dataurl=readGroupSet.getDataUrl(),
                indexfile=readGroupSet.getIndexFile(),
                attributes=json.dumps(readGroupSet.getAttributes()))
            for readGroup in readGroupSet.getReadGroups():
                self.insertReadGroup(readGroup)
        except Exception as e:
            raise exceptions.RepoManagerException(e)

    def removeReferenceSet(self, referenceSet):
        """
        Removes the specified referenceSet from this repository. This performs
        a cascading removal of all references within this referenceSet.
        However, it does not remove any of the ReadGroupSets or items that
        refer to this ReferenceSet. These must be deleted before the
        referenceSet can be removed.
        """
        try:
            q = models.Reference.delete().where(
                    models.Reference.referencesetid == referenceSet.getId())
            q.execute()
            q = models.Referenceset.delete().where(
                    models.Referenceset.id == referenceSet.getId())
            q.execute()
        except Exception:
            msg = ("Unable to delete reference set.  "
                   "There are objects currently in the registry which are "
                   "aligned against it.  Remove these objects before removing "
                   "the reference set.")
            raise exceptions.RepoManagerException(msg)

    def _readReadGroupSetTable(self):
        for readGroupSetRecord in models.Readgroupset.select():
            dataset = self.getDataset(readGroupSetRecord.datasetid.id)
            readGroupSet = reads.HtslibReadGroupSet(
                dataset, readGroupSetRecord.name)
            referenceSet = self.getReferenceSet(
                readGroupSetRecord.referencesetid.id)
            readGroupSet.setReferenceSet(referenceSet)
            readGroupSet.populateFromRow(readGroupSetRecord)
            assert readGroupSet.getId() == readGroupSetRecord.id
            # Insert the readGroupSet into the memory-based object model.
            dataset.addReadGroupSet(readGroupSet)

    def _createVariantAnnotationSetTable(self):
        self.database.create_table(models.Variantannotationset)

    def insertVariantAnnotationSet(self, variantAnnotationSet):
        """
        Inserts a the specified variantAnnotationSet into this repository.
        """
        analysisJson = json.dumps(
            protocol.toJsonDict(variantAnnotationSet.getAnalysis()))
        try:
            models.Variantannotationset.create(
                id=variantAnnotationSet.getId(),
                variantsetid=variantAnnotationSet.getParentContainer().getId(),
                ontologyid=variantAnnotationSet.getOntology().getId(),
                name=variantAnnotationSet.getLocalId(),
                analysis=analysisJson,
                annotationtype=variantAnnotationSet.getAnnotationType(),
                created=variantAnnotationSet.getCreationTime(),
                updated=variantAnnotationSet.getUpdatedTime(),
                attributes=json.dumps(variantAnnotationSet.getAttributes()))
        except Exception as e:
            raise exceptions.RepoManagerException(e)

    def _readVariantAnnotationSetTable(self):
        for annotationSetRecord in models.Variantannotationset.select():
            variantSet = self.getVariantSet(
                annotationSetRecord.variantsetid.id)
            ontology = self.getOntology(annotationSetRecord.ontologyid.id)
            variantAnnotationSet = variants.HtslibVariantAnnotationSet(
                variantSet, annotationSetRecord.name)
            variantAnnotationSet.setOntology(ontology)
            variantAnnotationSet.populateFromRow(annotationSetRecord)
            assert variantAnnotationSet.getId() == annotationSetRecord.id
            # Insert the variantAnnotationSet into the memory-based model.
            variantSet.addVariantAnnotationSet(variantAnnotationSet)

    def _createCallSetTable(self):
        self.database.create_table(models.Callset)

    def insertCallSet(self, callSet):
        """
        Inserts a the specified callSet into this repository.
        """
        try:
            models.Callset.create(
                id=callSet.getId(),
                name=callSet.getLocalId(),
                variantsetid=callSet.getParentContainer().getId(),
                biosampleid=callSet.getBiosampleId(),
                attributes=json.dumps(callSet.getAttributes()))
        except Exception as e:
            raise exceptions.RepoManagerException(e)

    def _readCallSetTable(self):
        for callSetRecord in models.Callset.select():
            variantSet = self.getVariantSet(callSetRecord.variantsetid.id)
            callSet = variants.CallSet(variantSet, callSetRecord.name)
            callSet.populateFromRow(callSetRecord)
            assert callSet.getId() == callSetRecord.id
            # Insert the callSet into the memory-based object model.
            variantSet.addCallSet(callSet)

    def _createVariantSetTable(self):
        self.database.create_table(models.Variantset)

    def insertVariantSet(self, variantSet):
        """
        Inserts a the specified variantSet into this repository.
        """
        # We cheat a little here with the VariantSetMetadata, and encode these
        # within the table as a JSON dump. These should really be stored in
        # their own table
        metadataJson = json.dumps(
            [protocol.toJsonDict(metadata) for metadata in
             variantSet.getMetadata()])
        urlMapJson = json.dumps(variantSet.getReferenceToDataUrlIndexMap())
        try:
            models.Variantset.create(
                id=variantSet.getId(),
                datasetid=variantSet.getParentContainer().getId(),
                referencesetid=variantSet.getReferenceSet().getId(),
                name=variantSet.getLocalId(),
                created=datetime.datetime.now(),
                updated=datetime.datetime.now(),
                metadata=metadataJson,
                dataurlindexmap=urlMapJson,
                attributes=json.dumps(variantSet.getAttributes()))
        except Exception as e:
            raise exceptions.RepoManagerException(e)
        for callSet in variantSet.getCallSets():
            self.insertCallSet(callSet)

    def _readVariantSetTable(self):
        for variantSetRecord in models.Variantset.select():
            dataset = self.getDataset(variantSetRecord.datasetid.id)
            referenceSet = self.getReferenceSet(
                variantSetRecord.referencesetid.id)
            variantSet = variants.HtslibVariantSet(
                dataset, variantSetRecord.name)
            variantSet.setReferenceSet(referenceSet)
            variantSet.populateFromRow(variantSetRecord)
            assert variantSet.getId() == variantSetRecord.id
            # Insert the variantSet into the memory-based object model.
            dataset.addVariantSet(variantSet)

    def _createFeatureSetTable(self):
        self.database.create_table(models.Featureset)

    def insertFeatureSet(self, featureSet):
        """
        Inserts a the specified featureSet into this repository.
        """
        # TODO add support for info and sourceUri fields.
        try:
            models.Featureset.create(
                id=featureSet.getId(),
                datasetid=featureSet.getParentContainer().getId(),
                referencesetid=featureSet.getReferenceSet().getId(),
                ontologyid=featureSet.getOntology().getId(),
                name=featureSet.getLocalId(),
                dataurl=featureSet.getDataUrl(),
                attributes=json.dumps(featureSet.getAttributes()))
        except Exception as e:
            raise exceptions.RepoManagerException(e)

    def _readFeatureSetTable(self):
        for featureSetRecord in models.Featureset.select():
            dataset = self.getDataset(featureSetRecord.datasetid.id)
            # FIXME this should be handled elsewhere
            if 'cgd' in featureSetRecord.name:
                featureSet = \
                    g2pFeatureset \
                    .PhenotypeAssociationFeatureSet(
                        dataset, featureSetRecord.name)
            else:
                featureSet = sequence_annotations.Gff3DbFeatureSet(
                    dataset, featureSetRecord.name)
            featureSet.setReferenceSet(
                self.getReferenceSet(
                    featureSetRecord.referencesetid.id))
            featureSet.setOntology(
                self.getOntology(featureSetRecord.ontologyid.id))
            featureSet.populateFromRow(featureSetRecord)
            assert featureSet.getId() == featureSetRecord.id
            dataset.addFeatureSet(featureSet)

    def _createContinuousSetTable(self):
        self.database.create_table(models.ContinuousSet)

    def insertContinuousSet(self, continuousSet):
        """
        Inserts a the specified continuousSet into this repository.
        """
        # TODO add support for info and sourceUri fields.
        try:
            models.ContinuousSet.create(
                id=continuousSet.getId(),
                datasetid=continuousSet.getParentContainer().getId(),
                referencesetid=continuousSet.getReferenceSet().getId(),
                name=continuousSet.getLocalId(),
                dataurl=continuousSet.getDataUrl(),
                attributes=json.dumps(continuousSet.getAttributes()))
        except Exception as e:
            raise exceptions.RepoManagerException(e)

    def _readContinuousSetTable(self):
        for continuousSetRecord in models.ContinuousSet.select():
            dataset = self.getDataset(continuousSetRecord.datasetid.id)
            continuousSet = continuous.FileContinuousSet(
                    dataset, continuousSetRecord.name)
            continuousSet.setReferenceSet(
                self.getReferenceSet(
                    continuousSetRecord.referencesetid.id))
            continuousSet.populateFromRow(continuousSetRecord)
            assert continuousSet.getId() == continuousSetRecord.id
            dataset.addContinuousSet(continuousSet)

    def _createBiosampleTable(self):
        self.database.create_table(models.Biosample)

    def insertBiosample(self, biosample):
        """
        Inserts the specified Biosample into this repository.
        """
        try:
            models.Biosample.create(
                id=biosample.getId(),
                datasetid=biosample.getParentContainer().getId(),
                name=biosample.getLocalId(),
                description=biosample.getDescription(),
                disease=json.dumps(biosample.getDisease()),
                created=biosample.getCreated(),
                updated=biosample.getUpdated(),
                individualid=biosample.getIndividualId(),
                attributes=json.dumps(biosample.getAttributes()),
                individualAgeAtCollection=json.dumps(
                        biosample.getIndividualAgeAtCollection()))
        except Exception:
            raise exceptions.DuplicateNameException(
                biosample.getLocalId(),
                biosample.getParentContainer().getLocalId())

    def _readBiosampleTable(self):
        for biosampleRecord in models.Biosample.select():
            dataset = self.getDataset(biosampleRecord.datasetid.id)
            biosample = biodata.Biosample(
                dataset, biosampleRecord.name)
            biosample.populateFromRow(biosampleRecord)
            assert biosample.getId() == biosampleRecord.id
            dataset.addBiosample(biosample)

    def _createIndividualTable(self):
        self.database.create_table(models.Individual)

    def insertIndividual(self, individual):
        """
        Inserts the specified individual into this repository.
        """
        try:
            models.Individual.create(
                id=individual.getId(),
                datasetId=individual.getParentContainer().getId(),
                name=individual.getLocalId(),
                description=individual.getDescription(),
                created=individual.getCreated(),
                updated=individual.getUpdated(),
                species=json.dumps(individual.getSpecies()),
                sex=json.dumps(individual.getSex()),
                attributes=json.dumps(individual.getAttributes()))
        except Exception:
            raise exceptions.DuplicateNameException(
                individual.getLocalId(),
                individual.getParentContainer().getLocalId())

    def _readIndividualTable(self):
        for individualRecord in models.Individual.select():
            dataset = self.getDataset(individualRecord.datasetid.id)
            individual = biodata.Individual(
                dataset, individualRecord.name)
            individual.populateFromRow(individualRecord)
            assert individual.getId() == individualRecord.id
            dataset.addIndividual(individual)

    def _createPhenotypeAssociationSetTable(self):
        self.database.create_table(models.Phenotypeassociationset)

    def _createRnaQuantificationSetTable(self):
        self.database.create_table(models.Rnaquantificationset)

    def insertPhenotypeAssociationSet(self, phenotypeAssociationSet):
        """
        Inserts the specified phenotype annotation set into this repository.
        """
        datasetId = phenotypeAssociationSet.getParentContainer().getId()
        attributes = json.dumps(phenotypeAssociationSet.getAttributes())
        try:
            models.Phenotypeassociationset.create(
                id=phenotypeAssociationSet.getId(),
                name=phenotypeAssociationSet.getLocalId(),
                datasetid=datasetId,
                dataurl=phenotypeAssociationSet._dataUrl,
                attributes=attributes)
        except Exception:
            raise exceptions.DuplicateNameException(
                phenotypeAssociationSet.getParentContainer().getId())

    def _readPhenotypeAssociationSetTable(self):
        for associationSetRecord in models.Phenotypeassociationset.select():
            dataset = self.getDataset(associationSetRecord.datasetid.id)
            phenotypeAssociationSet = \
                genotype_phenotype.RdfPhenotypeAssociationSet(
                    dataset,
                    associationSetRecord.name,
                    associationSetRecord.dataurl)
            dataset.addPhenotypeAssociationSet(phenotypeAssociationSet)

    def insertRnaQuantificationSet(self, rnaQuantificationSet):
        """
        Inserts a the specified rnaQuantificationSet into this repository.
        """
        try:
            models.Rnaquantificationset.create(
                id=rnaQuantificationSet.getId(),
                datasetid=rnaQuantificationSet.getParentContainer().getId(),
                referencesetid=rnaQuantificationSet.getReferenceSet().getId(),
                name=rnaQuantificationSet.getLocalId(),
                dataurl=rnaQuantificationSet.getDataUrl(),
                attributes=json.dumps(rnaQuantificationSet.getAttributes()))
        except Exception:
            raise exceptions.DuplicateNameException(
                rnaQuantificationSet.getLocalId(),
                rnaQuantificationSet.getParentContainer().getLocalId())

    def _readRnaQuantificationSetTable(self):
        for quantificationSetRecord in models.Rnaquantificationset.select():
            dataset = self.getDataset(quantificationSetRecord.datasetid.id)
            referenceSet = self.getReferenceSet(
                quantificationSetRecord.referencesetid.id)
            rnaQuantificationSet = \
                rna_quantification.SqliteRnaQuantificationSet(
                    dataset, quantificationSetRecord.name)
            rnaQuantificationSet.setReferenceSet(referenceSet)
            rnaQuantificationSet.populateFromRow(quantificationSetRecord)
            assert rnaQuantificationSet.getId() == quantificationSetRecord.id
            dataset.addRnaQuantificationSet(rnaQuantificationSet)

    def removeRnaQuantificationSet(self, rnaQuantificationSet):
        """
        Removes the specified rnaQuantificationSet from this repository. This
        performs a cascading removal of all items within this
        rnaQuantificationSet.
        """
        q = models.Rnaquantificationset.delete().where(
            models.Rnaquantificationset.id == rnaQuantificationSet.getId())
        q.execute()

    def insertPeer(self, peer):
        """
        Accepts a peer datamodel object and adds it to the registry.
        """
        try:
            models.Peer.create(
                url=peer.getUrl(),
                attributes=json.dumps(peer.getAttributes()))
        except Exception as e:
            raise exceptions.RepoManagerException(e)

    def removePeer(self, url):
        """
        Remove peers by URL.
        """
        q = models.Peer.delete().where(
            models.Peer.url == url)
        q.execute()

    def _createNetworkTables(self):
        """"""
        self.database.create_table(models.Peer)
        self.database.create_table(models.Announcement)

    def initialise(self):
        """
        Initialise this data repository, creating any necessary directories
        and file paths.
        """
        self._checkWriteMode()
        self._createSystemTable()
        self._createNetworkTables()
        self._createOntologyTable()
        self._createReferenceSetTable()
        self._createReferenceTable()
        self._createDatasetTable()
        self._createReadGroupSetTable()
        self._createReadGroupTable()
        self._createCallSetTable()
        self._createVariantSetTable()
        self._createVariantAnnotationSetTable()
        self._createFeatureSetTable()
        self._createContinuousSetTable()
        self._createBiosampleTable()
        self._createIndividualTable()
        self._createPhenotypeAssociationSetTable()
        self._createRnaQuantificationSetTable()

    def exists(self):
        """
        Checks that this data repository exists in the file system and has the
        required structure.
        """
        # TODO should this invoke a full load operation or just check the DB
        # exists?
        return os.path.exists(self._dbFilename)

    def assertExists(self):
        if not self.exists():
            raise exceptions.RepoNotFoundException(self._dbFilename)

    def delete(self):
        """
        Delete this data repository by recursively removing all directories.
        This will delete ALL data stored within the repository!!
        """
        os.unlink(self._dbFilename)

    def load(self):
        """
        Loads this data repository into memory.
        """
        self._readSystemTable()
        self._readOntologyTable()
        self._readReferenceSetTable()
        self._readReferenceTable()
        self._readDatasetTable()
        self._readReadGroupSetTable()
        self._readReadGroupTable()
        self._readVariantSetTable()
        self._readCallSetTable()
        self._readVariantAnnotationSetTable()
        self._readFeatureSetTable()
        self._readContinuousSetTable()
        self._readBiosampleTable()
        self._readIndividualTable()
        self._readPhenotypeAssociationSetTable()
        self._readRnaQuantificationSetTable()
