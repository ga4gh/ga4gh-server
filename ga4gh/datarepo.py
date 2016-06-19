"""
The backing data store for the GA4GH server
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import os
import sqlite3

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.ontologies as ontologies
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.references as references
import ga4gh.datamodel.variants as variants
import ga4gh.datamodel.sequenceAnnotations as sequenceAnnotations
import ga4gh.exceptions as exceptions
from ga4gh import protocol

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

    version = SchemaVersion("2.0")
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
        self._creationTimeStamp = None
        # Connection to the DB.
        self._dbConnection = None

    def _checkWriteMode(self):
        if self._openMode != MODE_WRITE:
            raise ValueError("Repo must be opened in write mode")

    def _checkReadMode(self):
        if self._openMode != MODE_READ:
            raise ValueError("Repo must be opened in read mode")

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
        self._safeConnect()
        # Turn on foreign key constraints.
        cursor = self._dbConnection.cursor()
        cursor.execute("PRAGMA foreign_keys = ON;")
        self._dbConnection.commit()
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
        self._dbConnection.commit()

    def close(self):
        """
        Closes this repo.
        """
        if self._openMode is None:
            raise ValueError("Repo already closed")
        self._openMode = None
        self._dbConnection.close()
        self._dbConnection = None

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

    def _safeConnect(self):
        try:
            # The next line creates the file if it did not exist previously
            self._dbConnection = sqlite3.connect(self._dbFilename)
        except sqlite3.OperationalError:
            # raised e.g. when directory passed as dbFilename
            raise exceptions.RepoInvalidDatabaseException(self._dbFilename)

    def _createSystemTable(self, cursor):
        sql = """
            CREATE TABLE System (
                key TEXT NOT NULL PRIMARY KEY,
                value TEXT NOT NULL
            );
        """
        cursor.execute(sql)
        cursor.execute(
            "INSERT INTO System VALUES "
            "('{}', '{}')".format(
                self.systemKeySchemaVersion, self.version))
        cursor.execute(
            "INSERT INTO System VALUES ('{}', datetime('now'))".format(
                self.systemKeyCreationTimeStamp))

    def _readSystemTable(self, cursor):
        sql = "SELECT key, value FROM System;"
        cursor.execute(sql)
        config = {}
        for row in cursor:
            config[row[0]] = row[1]
        row = cursor.fetchone()
        self._schemaVersion = config[self.systemKeySchemaVersion]
        self._creationTimeStamp = config[self.systemKeyCreationTimeStamp]
        schemaVersion = self.SchemaVersion(self._schemaVersion)
        if schemaVersion.major != self.version.major:
            raise exceptions.RepoSchemaVersionMismatchException(
                schemaVersion, self.version)

    def _createOntologyTable(self, cursor):
        sql = """
            CREATE TABLE Ontology(
                id TEXT NOT NULL PRIMARY KEY,
                name TEXT NOT NULL,
                dataUrl TEXT NOT NULL,
                ontologyPrefix TEXT NOT NULL,
                UNIQUE (name)
            );
        """
        cursor.execute(sql)

    def insertOntology(self, ontology):
        """
        Inserts the specified ontology into this repository.
        """
        sql = """
            INSERT INTO Ontology(id, name, dataUrl, ontologyPrefix)
            VALUES (?, ?, ?, ?);
        """
        cursor = self._dbConnection.cursor()
        # TODO we need to create a proper ID when we're doing ID generation
        # for the rest of the container objects.
        try:
            cursor.execute(sql, (
                ontology.getName(),
                ontology.getName(),
                ontology.getDataUrl(),
                ontology.getOntologyPrefix()))
        except sqlite3.IntegrityError:
            raise exceptions.DuplicateNameException(ontology.getName())

    def _readOntologyTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM Ontology;")
        for row in cursor:
            ontology = ontologies.Ontology(row[b'name'])
            ontology.populateFromRow(row)
            self.addOntology(ontology)

    def removeOntology(self, ontology):
        """
        Removes the specified ontology term map from this repository.
        """
        sql = "DELETE FROM Ontology WHERE name=?"
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (ontology.getName(),))

    def _createReferenceTable(self, cursor):
        sql = """
            CREATE TABLE Reference (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                referenceSetId TEXT NOT NULL,
                length INTEGER,
                isDerived INTEGER,
                md5checksum TEXT,
                ncbiTaxonId INTEGER,
                sourceAccessions TEXT,
                sourceDivergence REAL,
                sourceUri TEXT,
                UNIQUE (referenceSetId, name),
                FOREIGN KEY(referenceSetId) REFERENCES ReferenceSet(id)
                    ON DELETE CASCADE
            );
        """
        cursor.execute(sql)

    def insertReference(self, reference):
        """
        Inserts the specified reference into this repository.
        """
        sql = """
            INSERT INTO Reference (
                id, referenceSetId, name, length, isDerived, md5checksum,
                ncbiTaxonId, sourceAccessions, sourceUri)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);
        """
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (
            reference.getId(), reference.getParentContainer().getId(),
            reference.getLocalId(), reference.getLength(),
            reference.getIsDerived(), reference.getMd5Checksum(),
            reference.getNcbiTaxonId(),
            # We store the list of sourceAccessions as a JSON string. Perhaps
            # this should be another table?
            json.dumps(reference.getSourceAccessions()),
            reference.getSourceUri()))

    def _readReferenceTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM Reference;")
        for row in cursor:
            referenceSet = self.getReferenceSet(row[b'referenceSetId'])
            reference = references.HtslibReference(referenceSet, row[b'name'])
            reference.populateFromRow(row)
            assert reference.getId() == row[b"id"]
            referenceSet.addReference(reference)

    def _createReferenceSetTable(self, cursor):
        sql = """
            CREATE TABLE ReferenceSet (
                id TEXT NOT NULL PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT,
                assemblyId TEXT,
                isDerived INTEGER,
                md5checksum TEXT,
                ncbiTaxonId INTEGER,
                sourceAccessions TEXT,
                sourceUri TEXT,
                dataUrl TEXT NOT NULL,
                UNIQUE (name)
            );
        """
        cursor.execute(sql)

    def insertReferenceSet(self, referenceSet):
        """
        Inserts the specified referenceSet into this repository.
        """
        sql = """
            INSERT INTO ReferenceSet (
                id, name, description, assemblyId, isDerived, md5checksum,
                ncbiTaxonId, sourceAccessions, sourceUri, dataUrl)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
        """
        cursor = self._dbConnection.cursor()
        try:
            cursor.execute(sql, (
                referenceSet.getId(), referenceSet.getLocalId(),
                referenceSet.getDescription(), referenceSet.getAssemblyId(),
                referenceSet.getIsDerived(), referenceSet.getMd5Checksum(),
                referenceSet.getNcbiTaxonId(),
                # We store the list of sourceAccessions as a JSON string.
                # Perhaps this should be another table?
                json.dumps(referenceSet.getSourceAccessions()),
                referenceSet.getSourceUri(), referenceSet.getDataUrl()))
        except sqlite3.IntegrityError:
            raise exceptions.DuplicateNameException(referenceSet.getLocalId())
        for reference in referenceSet.getReferences():
            self.insertReference(reference)

    def _readReferenceSetTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM ReferenceSet;")
        for row in cursor:
            referenceSet = references.HtslibReferenceSet(row[b'name'])
            referenceSet.populateFromRow(row)
            assert referenceSet.getId() == row[b"id"]
            # Insert the referenceSet into the memory-based object model.
            self.addReferenceSet(referenceSet)

    def _createDatasetTable(self, cursor):
        sql = """
            CREATE TABLE Dataset (
                id TEXT NOT NULL PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT,
                UNIQUE (name)
            );
        """
        cursor.execute(sql)

    def insertDataset(self, dataset):
        """
        Inserts the specified dataset into this repository.
        """
        sql = """
            INSERT INTO Dataset (id, name, description)
            VALUES (?, ?, ?);
        """
        cursor = self._dbConnection.cursor()
        try:
            cursor.execute(sql, (
                dataset.getId(), dataset.getLocalId(),
                dataset.getDescription()))
        except sqlite3.IntegrityError:
            raise exceptions.DuplicateNameException(dataset.getLocalId())

    def removeDataset(self, dataset):
        """
        Removes the specified dataset from this repository. This performs
        a cascading removal of all items within this dataset.
        """
        sql = "DELETE FROM Dataset WHERE id=?"
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (dataset.getId(),))

    def removeFeatureSet(self, featureSet):
        """
        Removes the specified featureSet from this repository.
        """
        sql = "DELETE FROM FeatureSet WHERE id=?"
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (featureSet.getId(),))

    def _readDatasetTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM Dataset;")
        for row in cursor:
            dataset = datasets.Dataset(row[b'name'])
            dataset.populateFromRow(row)
            assert dataset.getId() == row[b"id"]
            # Insert the dataset into the memory-based object model.
            self.addDataset(dataset)

    def _createReadGroupTable(self, cursor):
        sql = """
            CREATE TABLE ReadGroup (
                id TEXT NOT NULL PRIMARY KEY,
                readGroupSetId TEXT NOT NULL,
                name TEXT NOT NULL,
                predictedInsertSize INTEGER,
                sampleId TEXT,
                description TEXT,
                stats TEXT NOT NULL,
                experiment TEXT NOT NULL,
                created TEXT,
                updated TEXT,
                UNIQUE (readGroupSetId, name),
                FOREIGN KEY(readGroupSetId) REFERENCES ReadGroupSet(id)
                    ON DELETE CASCADE
            );
        """
        cursor.execute(sql)

    def insertReadGroup(self, readGroup):
        """
        Inserts the specified readGroup into the DB.
        """
        sql = """
            INSERT INTO ReadGroup (
                id, readGroupSetId, name, predictedInsertSize,
                sampleId, description, stats, experiment, created, updated)
            VALUES
                (?, ?, ?, ?, ?, ?, ?, ?, datetime('now'), datetime('now'));
        """
        cursor = self._dbConnection.cursor()
        statsJson = json.dumps(protocol.toJsonDict(readGroup.getStats()))
        experimentJson = json.dumps(
            protocol.toJsonDict(readGroup.getExperiment()))
        cursor.execute(sql, (
            readGroup.getId(), readGroup.getParentContainer().getId(),
            readGroup.getLocalId(), readGroup.getPredictedInsertSize(),
            readGroup.getSampleId(), readGroup.getDescription(),
            statsJson, experimentJson))

    def removeReadGroupSet(self, readGroupSet):
        """
        Removes the specified readGroupSet from this repository. This performs
        a cascading removal of all items within this readGroupSet.
        """
        sql = "DELETE FROM ReadGroupSet WHERE id=?"
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (readGroupSet.getId(),))

    def removeVariantSet(self, variantSet):
        """
        Removes the specified variantSet from this repository. This performs
        a cascading removal of all items within this variantSet.
        """
        sql = "DELETE FROM VariantSet WHERE id=?"
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (variantSet.getId(),))

    def _readReadGroupTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM ReadGroup;")
        for row in cursor:
            readGroupSet = self.getReadGroupSet(row[b'readGroupSetId'])
            readGroup = reads.HtslibReadGroup(readGroupSet, row[b'name'])
            # TODO set the reference set.
            readGroup.populateFromRow(row)
            assert readGroup.getId() == row[b'id']
            # Insert the readGroupSet into the memory-based object model.
            readGroupSet.addReadGroup(readGroup)

    def _createReadGroupSetTable(self, cursor):
        sql = """
            CREATE TABLE ReadGroupSet (
                id TEXT NOT NULL PRIMARY KEY,
                name TEXT NOT NULL,
                datasetId TEXT NOT NULL,
                referenceSetId TEXT NOT NULL,
                programs TEXT,
                stats TEXT NOT NULL,
                dataUrl TEXT NOT NULL,
                indexFile TEXT NOT NULL,
                UNIQUE (datasetId, name),
                FOREIGN KEY(datasetId) REFERENCES Dataset(id)
                    ON DELETE CASCADE,
                FOREIGN KEY(referenceSetId) REFERENCES ReferenceSet(id)
            );
        """
        cursor.execute(sql)

    def insertReadGroupSet(self, readGroupSet):
        """
        Inserts a the specified readGroupSet into this repository.
        """
        sql = """
            INSERT INTO ReadGroupSet (
                id, datasetId, referenceSetId, name, programs, stats,
                dataUrl, indexFile)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?);
        """
        programsJson = json.dumps(
            [protocol.toJsonDict(program) for program in
             readGroupSet.getPrograms()])
        statsJson = json.dumps(protocol.toJsonDict(readGroupSet.getStats()))
        cursor = self._dbConnection.cursor()
        try:
            cursor.execute(sql, (
                readGroupSet.getId(),
                readGroupSet.getParentContainer().getId(),
                readGroupSet.getReferenceSet().getId(),
                readGroupSet.getLocalId(),
                programsJson, statsJson, readGroupSet.getDataUrl(),
                readGroupSet.getIndexFile()))
        except sqlite3.IntegrityError:
            raise exceptions.DuplicateNameException(
                readGroupSet.getLocalId(),
                readGroupSet.getParentContainer().getLocalId())
        for readGroup in readGroupSet.getReadGroups():
            self.insertReadGroup(readGroup)

    def removeReferenceSet(self, referenceSet):
        """
        Removes the specified referenceSet from this repository. This performs
        a cascading removal of all references within this referenceSet.
        However, it does not remove any of the ReadGroupSets or items that
        refer to this ReferenceSet. These must be deleted before the
        referenceSet can be removed.
        """
        sql = "DELETE FROM ReferenceSet WHERE id=?"
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (referenceSet.getId(),))

    def _readReadGroupSetTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM ReadGroupSet;")
        for row in cursor:
            dataset = self.getDataset(row[b'datasetId'])
            readGroupSet = reads.HtslibReadGroupSet(dataset, row[b'name'])
            referenceSet = self.getReferenceSet(row[b'referenceSetId'])
            readGroupSet.setReferenceSet(referenceSet)
            readGroupSet.populateFromRow(row)
            assert readGroupSet.getId() == row[b'id']
            # Insert the readGroupSet into the memory-based object model.
            dataset.addReadGroupSet(readGroupSet)

    def _createVariantAnnotationSetTable(self, cursor):
        sql = """
            CREATE TABLE VariantAnnotationSet (
                id TEXT NOT NULL PRIMARY KEY,
                name TEXT NOT NULL,
                variantSetId TEXT NOT NULL,
                ontologyId TEXT NOT NULL,
                analysis TEXT,
                annotationType TEXT,
                created TEXT,
                updated TEXT,
                UNIQUE (variantSetId, name),
                FOREIGN KEY(variantSetId) REFERENCES VariantSet(id)
                    ON DELETE CASCADE,
                FOREIGN KEY(ontologyId) REFERENCES Ontology(id)
            );
        """
        cursor.execute(sql)

    def insertVariantAnnotationSet(self, variantAnnotationSet):
        """
        Inserts a the specified variantAnnotationSet into this repository.
        """
        sql = """
            INSERT INTO VariantAnnotationSet (
                id, variantSetId, ontologyId, name, analysis, annotationType,
                created, updated)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?);
        """
        analysisJson = json.dumps(
            protocol.toJsonDict(variantAnnotationSet.getAnalysis()))
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (
            variantAnnotationSet.getId(),
            variantAnnotationSet.getParentContainer().getId(),
            variantAnnotationSet.getOntology().getId(),
            variantAnnotationSet.getLocalId(),
            analysisJson,
            variantAnnotationSet.getAnnotationType(),
            variantAnnotationSet.getCreationTime(),
            variantAnnotationSet.getUpdatedTime()))

    def _readVariantAnnotationSetTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM VariantAnnotationSet;")
        for row in cursor:
            variantSet = self.getVariantSet(row[b'variantSetId'])
            ontology = self.getOntology(row[b'ontologyId'])
            variantAnnotationSet = variants.HtslibVariantAnnotationSet(
                variantSet, row[b'name'])
            variantAnnotationSet.setOntology(ontology)
            variantAnnotationSet.populateFromRow(row)
            assert variantAnnotationSet.getId() == row[b'id']
            # Insert the variantAnnotationSet into the memory-based model.
            variantSet.addVariantAnnotationSet(variantAnnotationSet)

    def _createCallSetTable(self, cursor):
        sql = """
            CREATE TABLE CallSet (
                id TEXT NOT NULL PRIMARY KEY,
                name TEXT NOT NULL,
                variantSetId TEXT NOT NULL,
                UNIQUE (variantSetId, name),
                FOREIGN KEY(variantSetId) REFERENCES VariantSet(id)
                    ON DELETE CASCADE
            );
        """
        cursor.execute(sql)

    def insertCallSet(self, callSet):
        """
        Inserts a the specified callSet into this repository.
        """
        sql = """
            INSERT INTO CallSet (
                id, variantSetId, name)
            VALUES (?, ?, ?);
        """
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (
            callSet.getId(), callSet.getParentContainer().getId(),
            callSet.getLocalId()))

    def _readCallSetTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM CallSet;")
        for row in cursor:
            variantSet = self.getVariantSet(row[b'variantSetId'])
            callSet = variants.CallSet(variantSet, row[b'name'])
            callSet.populateFromRow(row)
            assert callSet.getId() == row[b'id']
            # Insert the callSet into the memory-based object model.
            variantSet.addCallSet(callSet)

    def _createVariantSetTable(self, cursor):
        sql = """
            CREATE TABLE VariantSet (
                id TEXT NOT NULL PRIMARY KEY,
                name TEXT NOT NULL,
                datasetId TEXT NOT NULL,
                referenceSetId TEXT NOT NULL,
                created TEXT,
                updated TEXT,
                metadata TEXT,
                dataUrlIndexMap TEXT NOT NULL,
                UNIQUE (datasetID, name),
                FOREIGN KEY(datasetId) REFERENCES Dataset(id)
                    ON DELETE CASCADE,
                FOREIGN KEY(referenceSetId) REFERENCES ReferenceSet(id)
            );
        """
        cursor.execute(sql)

    def insertVariantSet(self, variantSet):
        """
        Inserts a the specified variantSet into this repository.
        """
        sql = """
            INSERT INTO VariantSet (
                id, datasetId, referenceSetId, name, created, updated,
                metadata, dataUrlIndexMap)
            VALUES (?, ?, ?, ?, datetime('now'), datetime('now'), ?, ?);
        """
        cursor = self._dbConnection.cursor()
        # We cheat a little here with the VariantSetMetadata, and encode these
        # within the table as a JSON dump. These should really be stored in
        # their own table
        metadataJson = json.dumps(
            [protocol.toJsonDict(metadata) for metadata in
             variantSet.getMetadata()])
        urlMapJson = json.dumps(variantSet.getReferenceToDataUrlIndexMap())
        try:
            cursor.execute(sql, (
                variantSet.getId(), variantSet.getParentContainer().getId(),
                variantSet.getReferenceSet().getId(), variantSet.getLocalId(),
                metadataJson, urlMapJson))
        except sqlite3.IntegrityError:
            raise exceptions.DuplicateNameException(
                variantSet.getLocalId(),
                variantSet.getParentContainer().getLocalId())
        for callSet in variantSet.getCallSets():
            self.insertCallSet(callSet)

    def _readVariantSetTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM VariantSet;")
        for row in cursor:
            dataset = self.getDataset(row[b'datasetId'])
            referenceSet = self.getReferenceSet(row[b'referenceSetId'])
            variantSet = variants.HtslibVariantSet(dataset, row[b'name'])
            variantSet.setReferenceSet(referenceSet)
            variantSet.populateFromRow(row)
            assert variantSet.getId() == row[b'id']
            # Insert the variantSet into the memory-based object model.
            dataset.addVariantSet(variantSet)

    def _createFeatureSetTable(self, cursor):
        sql = """
            CREATE TABLE FeatureSet (
                id TEXT NOT NULL PRIMARY KEY,
                name TEXT NOT NULL,
                datasetId TEXT NOT NULL,
                referenceSetId TEXT NOT NULL,
                ontologyId TEXT NOT NULL,
                info TEXT,
                sourceUri TEXT,
                dataUrl TEXT NOT NULL,
                UNIQUE (datasetId, name),
                FOREIGN KEY(datasetId) REFERENCES Dataset(id)
                    ON DELETE CASCADE,
                FOREIGN KEY(referenceSetId) REFERENCES ReferenceSet(id)
                FOREIGN KEY(ontologyId) REFERENCES Ontology(id)
            );
        """
        cursor.execute(sql)

    def insertFeatureSet(self, featureSet):
        """
        Inserts a the specified featureSet into this repository.
        """
        # TODO add support for info and sourceUri fields.
        sql = """
            INSERT INTO FeatureSet (
                id, datasetId, referenceSetId, ontologyId, name, dataUrl)
            VALUES (?, ?, ?, ?, ?, ?)
        """
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (
            featureSet.getId(),
            featureSet.getParentContainer().getId(),
            featureSet.getReferenceSet().getId(),
            featureSet.getOntology().getId(),
            featureSet.getLocalId(),
            featureSet.getDataUrl()))

    def _readFeatureSetTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM FeatureSet;")
        for row in cursor:
            dataset = self.getDataset(row[b'datasetId'])
            featureSet = sequenceAnnotations.Gff3DbFeatureSet(
                dataset, row[b'name'])
            featureSet.setReferenceSet(
                self.getReferenceSet(row[b'referenceSetId']))
            featureSet.setOntology(self.getOntology(row[b'ontologyId']))
            featureSet.populateFromRow(row)
            assert featureSet.getId() == row[b'id']
            dataset.addFeatureSet(featureSet)

    def initialise(self):
        """
        Initialise this data repostitory, creating any necessary directories
        and file paths.
        """
        self._checkWriteMode()
        cursor = self._dbConnection
        self._createSystemTable(cursor)
        self._createOntologyTable(cursor)
        self._createReferenceSetTable(cursor)
        self._createReferenceTable(cursor)
        self._createDatasetTable(cursor)
        self._createReadGroupSetTable(cursor)
        self._createReadGroupTable(cursor)
        self._createCallSetTable(cursor)
        self._createVariantSetTable(cursor)
        self._createVariantAnnotationSetTable(cursor)
        self._createFeatureSetTable(cursor)

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
        with sqlite3.connect(self._dbFilename) as db:
            cursor = db.cursor()
            try:
                self._readSystemTable(cursor)
            except (sqlite3.OperationalError, sqlite3.DatabaseError):
                raise exceptions.RepoInvalidDatabaseException(
                    self._dbFilename)
            self._readOntologyTable(cursor)
            self._readReferenceSetTable(cursor)
            self._readReferenceTable(cursor)
            self._readDatasetTable(cursor)
            self._readReadGroupSetTable(cursor)
            self._readReadGroupTable(cursor)
            self._readVariantSetTable(cursor)
            self._readCallSetTable(cursor)
            self._readVariantAnnotationSetTable(cursor)
            self._readFeatureSetTable(cursor)
