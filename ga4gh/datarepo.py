"""
The backing data store for the GA4GH server
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import glob
import os
import sqlite3

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.ontologies as ontologies
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.references as references
import ga4gh.datamodel.variants as variants
import ga4gh.exceptions as exceptions


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
        Add an ontology map to this data repository.
        """
        name = ontology.getLocalId()
        self._ontologyNameMap[name] = ontology
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
        Returns an ontology map by name
        """
        return self._ontologyNameMap.get(name, None)

    def getOntologys(self):
        """
        Returns all ontology maps in the repo
        """
        return [self._ontologyNameMap[name] for name in self._ontologyNames]

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
        dataset = self.getDataset(compoundId.datasetId)
        return dataset.getReadGroupSet(id_)

    def getVariantSet(self, id_):
        """
        Returns the readgroup set with the specified ID.
        """
        compoundId = datamodel.VariantSetCompoundId.parse(id_)
        dataset = self.getDataset(compoundId.datasetId)
        return dataset.getVariantSet(id_)


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
    def __init__(self, fileName):
        super(SqlDataRepository, self).__init__()
        self._dbFilename = fileName
        # Values filled in using the DB. These will all be None until
        # we have called load()
        self._repositoryVersion = None
        self._creationTimeStamp = None
        # Connection to the DB.
        self._dbConnection = None

    def openDb(self):
        """
        Opens the DB connection to make the repository ready for updates.
        Must be called before inserting any objects with the 'insertX'
        methods.
        """
        self._dbConnection = sqlite3.connect(self._dbFilename)

    def closeDb(self):
        """
        Closes the DB connection and commits all changes made.
        """
        self._dbConnection.commit()
        self._dbConnection.close()
        self._dbConnection = None

    def _createSystemTable(self, cursor):
        sql = """
            CREATE TABLE System (
                repositoryVersion TEXT,
                creationTimeStamp TEXT
            );
        """
        cursor.execute(sql)
        cursor.execute("INSERT INTO System VALUES (0.1, datetime('now'));")

    def _readSystemTable(self, cursor):
        sql = "SELECT repositoryVersion, creationTimeStamp FROM System;"
        cursor.execute(sql)
        row = cursor.fetchone()
        self._repositoryVersion = row[0]
        self._creationTimeStamp = row[1]

    def _createOntologyTable(self, cursor):
        sql = """
            CREATE TABLE Ontology (
                name TEXT,
                dataUrl TEXT
            );
        """
        cursor.execute(sql)

    def insertOntology(self, ontology):
        """
        Inserts the specified ontology into this repository.
        """
        sql = """
            INSERT INTO Ontology (name, dataUrl)
            VALUES (?, ?);
        """
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (
            ontology.getLocalId(), ontology.getDataUrl()))

    def _readOntologyTable(self, cursor):
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM Ontology;")
        for row in cursor:
            ontology = ontologies.Ontology(row[b'name'])
            ontology.populateFromRow(row)
            self.addOntology(ontology)

    def _createReferenceTable(self, cursor):
        sql = """
            CREATE TABLE Reference (
                id TEXT,
                referenceSetId TEXT,
                name TEXT,
                length INTEGER,
                isDerived INTEGER,
                md5checksum TEXT,
                ncbiTaxonId TEXT,
                sourceAccessions TEXT,
                sourceDivergence REAL,
                sourceUri TEXT
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
                id TEXT,
                name TEXT,
                description TEXT,
                assemblyId TEXT,
                isDerived INTEGER,
                md5checksum TEXT,
                ncbiTaxonId TEXT,
                sourceAccessions TEXT,
                sourceUri TEXT,
                dataUrl TEXT
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
        cursor.execute(sql, (
            referenceSet.getId(), referenceSet.getLocalId(),
            referenceSet.getDescription(), referenceSet.getAssemblyId(),
            referenceSet.getIsDerived(), referenceSet.getMd5Checksum(),
            referenceSet.getNcbiTaxonId(),
            # We store the list of sourceAccessions as a JSON string. Perhaps
            # this should be another table?
            json.dumps(referenceSet.getSourceAccessions()),
            referenceSet.getSourceUri(), referenceSet.getDataUrl()))
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
                id TEXT,
                name TEXT,
                description TEXT
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
        cursor.execute(sql, (
            dataset.getId(), dataset.getLocalId(), dataset.getDescription()))

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
                id TEXT,
                readGroupSetId TEXT,
                name TEXT,
                predictedInsertSize INTEGER,
                sampleId TEXT,
                description TEXT,
                created TEXT,
                updated TEXT
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
                sampleId, created, updated)
            VALUES (?, ?, ?, ?, ?, datetime('now'), datetime('now'));
        """
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (
            readGroup.getId(), readGroup.getParentContainer().getId(),
            readGroup.getLocalId(), readGroup.getPredictedInsertSize(),
            readGroup.getSampleId()))

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
                id TEXT,
                datasetId TEXT,
                referenceSetId TEXT,
                name TEXT,
                programs TEXT,
                dataUrl TEXT,
                indexFile TEXT
            );
        """
        cursor.execute(sql)

    def insertReadGroupSet(self, readGroupSet):
        """
        Inserts a the specified readGroupSet into this repository.
        """
        sql = """
            INSERT INTO ReadGroupSet (
                id, datasetId, referenceSetId, name, programs,
                dataUrl, indexFile)
            VALUES (?, ?, ?, ?, ?, ?, ?);
        """
        programsJson = json.dumps(
            [program.toJsonDict() for program in readGroupSet.getPrograms()])
        cursor = self._dbConnection.cursor()
        cursor.execute(sql, (
            readGroupSet.getId(), readGroupSet.getParentContainer().getId(),
            readGroupSet.getReferenceSet().getId(), readGroupSet.getLocalId(),
            programsJson, readGroupSet.getDataUrl(),
            readGroupSet.getIndexFile()))
        for readGroup in readGroupSet.getReadGroups():
            self.insertReadGroup(readGroup)

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

    def _createCallSetTable(self, cursor):
        sql = """
            CREATE TABLE CallSet (
                id TEXT,
                variantSetId TEXT,
                name TEXT
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
                id TEXT,
                datasetId TEXT,
                referenceSetId TEXT,
                name TEXT,
                created TEXT,
                updated TEXT,
                metadata TEXT,
                isAnnotated INTEGER,
                dataUrlIndexMap TEXT
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
                metadata, isAnnotated, dataUrlIndexMap)
            VALUES (?, ?, ?, ?, datetime('now'), datetime('now'), ?, ?, ?);
        """
        cursor = self._dbConnection.cursor()
        # We cheat a little here with the VariantSetMetadata, and encode these
        # within the table as a JSON dump. These should really be stored in
        # their own table
        metadataJson = json.dumps(
            [metadata.toJsonDict() for metadata in variantSet.getMetadata()])
        urlMapJson = json.dumps(variantSet.getReferenceToDataUrlIndexMap())
        cursor.execute(sql, (
            variantSet.getId(), variantSet.getParentContainer().getId(),
            variantSet.getReferenceSet().getId(), variantSet.getLocalId(),
            metadataJson, variantSet.isAnnotated(), urlMapJson))
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
            # If it is annotated, also create a VariantAnnotationSet
            if variantSet.isAnnotated():
                sequenceOntology = self.getOntology('sequence_ontology')
                variantAnnotationSet = variants.HtslibVariantAnnotationSet(
                    variantSet, sequenceOntology)
                variantAnnotationSet.populateFromFile()
                dataset.addVariantAnnotationSet(variantAnnotationSet)

    def initialise(self):
        """
        Initialise this data repostitory, creating any necessary directories
        and file paths.
        """
        # Create the SQLite DB
        with sqlite3.connect(self._dbFilename) as db:
            cursor = db.cursor()
            self._createSystemTable(cursor)
            self._createOntologyTable(cursor)
            self._createReferenceSetTable(cursor)
            self._createReferenceTable(cursor)
            self._createDatasetTable(cursor)
            self._createReadGroupSetTable(cursor)
            self._createReadGroupTable(cursor)
            self._createCallSetTable(cursor)
            self._createVariantSetTable(cursor)
            # TODO add ReadStats, Experiment and other ReadGroup tables.

    def exists(self):
        """
        Checks that this data repository exists in the file system and has the
        required structure.
        """
        # TODO should this invoke a full load operation or just check the DB
        # exists?
        return os.path.exists(self._dbFilename)

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
            self._readSystemTable(cursor)
            self._readOntologyTable(cursor)
            self._readReferenceSetTable(cursor)
            self._readReferenceTable(cursor)
            self._readDatasetTable(cursor)
            self._readReadGroupSetTable(cursor)
            self._readReadGroupTable(cursor)
            self._readVariantSetTable(cursor)
            self._readCallSetTable(cursor)

    def printSummary(self):
        """
        Prints a summary of this data repository to stdout.
        """
        print("Repository version {} at path '{}'".format(
            self._repositoryVersion, self._dbFilename))
        print("Created at ", self._creationTimeStamp)
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
                    "\t", readGroupSet.getLocalId(), readGroupSet.getId(),
                    readGroupSet.getDataUrl(), sep="\t")
                for readGroup in readGroupSet.getReadGroups():
                    print(
                        "\t\t", readGroup.getId(), readGroup.getLocalId(),
                        sep="\t")
            print("\tVariantSets:")
            for variantSet in dataset.getVariantSets():
                print(
                    "\t", variantSet.getLocalId(), variantSet.getId(),
                    sep="\t")


class FileSystemDataRepository(AbstractDataRepository):
    """
    Class representing and old-style FileSystem based data repository.
    This is primarily intended to provide an easy way to keep existing
    tests working, and may also provide a smooth upgrade path for any
    users who have data stored in the old file system based repos.

    This is deprecated and should be removed , along with updating
    the necessary tests.
    """
    referenceSetsDirName = "referenceSets"
    datasetsDirName = "datasets"
    ontologiesDirName = "ontologymaps"

    def __init__(self, dataDir):
        super(FileSystemDataRepository, self).__init__()
        self._dataDir = dataDir

        pattern = os.path.join(
            self._dataDir, self.referenceSetsDirName, "*.fa.gz")
        for fastaFile in glob.glob(pattern):
            name = os.path.basename(fastaFile).split(".")[0]
            referenceSet = references.HtslibReferenceSet(name)
            referenceSet.populateFromFile(fastaFile)
            self.addReferenceSet(referenceSet)

        # for ontologies we go into each directory within the
        # main directory.
        sourceDir = os.path.join(self._dataDir, self.ontologiesDirName)
        for setName in os.listdir(sourceDir):
            relativePath = os.path.join(sourceDir, setName)
            if os.path.isdir(relativePath):
                for filename in os.listdir(relativePath):
                    if filename.endswith(".txt"):
                        name = filename.split(".")[0]
                        ontology = ontologies.Ontology(name)
                        ontology.populateFromFile(
                            os.path.join(relativePath, filename))
                        self.addOntology(ontology)

        sourceDir = os.path.join(self._dataDir, self.datasetsDirName)
        for setName in os.listdir(sourceDir):
            relativePath = os.path.join(sourceDir, setName)
            if os.path.isdir(relativePath):
                dataset = datasets.FileSystemDataset(
                    setName, relativePath, self)
                self.addDataset(dataset)

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
