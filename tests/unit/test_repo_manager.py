"""
Tests for the repo manager tool
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import shutil
import tempfile
import unittest

import ga4gh.exceptions as exceptions
import ga4gh.datarepo as datarepo
import ga4gh.cli as cli
import tests.paths as paths


class TestGetNameFromPath(unittest.TestCase):
    """
    Tests the method for deriving the default name of objects from file
    paths.
    """
    def testError(self):
        self.assertRaises(ValueError, cli.getNameFromPath, "")

    def testLocalDirectory(self):
        self.assertEqual(cli.getNameFromPath("no_extension"), "no_extension")
        self.assertEqual(cli.getNameFromPath("x.y"), "x")
        self.assertEqual(cli.getNameFromPath("x.y.z"), "x")

    def testFullPaths(self):
        self.assertEqual(cli.getNameFromPath("/no_ext"), "no_ext")
        self.assertEqual(cli.getNameFromPath("/x.y"), "x")
        self.assertEqual(cli.getNameFromPath("/x.y.z"), "x")
        self.assertEqual(cli.getNameFromPath("/a/no_ext"), "no_ext")
        self.assertEqual(cli.getNameFromPath("/a/x.y"), "x")
        self.assertEqual(cli.getNameFromPath("/a/x.y.z"), "x")

    def testUrls(self):
        self.assertEqual(cli.getNameFromPath("file:///no_ext"), "no_ext")
        self.assertEqual(cli.getNameFromPath("http://example.com/x.y"), "x")
        self.assertEqual(cli.getNameFromPath("ftp://x.y.z"), "x")


class AbstractRepoManagerTest(unittest.TestCase):
    """
    Base class for repo manager tests
    """
    def setUp(self):
        fd, self._repoPath = tempfile.mkstemp(prefix="ga4gh_repoman_test")
        os.unlink(self._repoPath)

    def runCommand(self, cmd):
        cli.RepoManager.runCommand(cmd.split())

    def tearDown(self):
        os.unlink(self._repoPath)

    def readRepo(self):
        repo = datarepo.SqlDataRepository(self._repoPath)
        repo.open(datarepo.MODE_READ)
        return repo

    def init(self):
        self.runCommand("init {}".format(self._repoPath))

    def addOntology(self):
        # Add the sequence ontology
        self._ontologyName = paths.ontologyName
        cmd = "add-ontology {} {}".format(self._repoPath, paths.ontologyPath)
        self.runCommand(cmd)

    def addDataset(self):
        self._datasetName = "test_dataset"
        cmd = "add-dataset {} {}".format(self._repoPath, self._datasetName)
        self.runCommand(cmd)

    def addReferenceSet(self):
        self._referenceSetName = "test_rs"
        fastaFile = paths.ncbi37FaPath
        self.runCommand("add-referenceset {} {} --name={}".format(
            self._repoPath, fastaFile, self._referenceSetName))

    def addReadGroupSet(self):
        bamFile = paths.bamPath
        self._readGroupSetName = "test_rgs"
        cmd = (
            "add-readgroupset {} {} {} --referenceSetName={} "
            "--name={}").format(
            self._repoPath, self._datasetName, bamFile,
            self._referenceSetName, self._readGroupSetName)
        self.runCommand(cmd)

    def addFeatureSet(self):
        featuresPath = paths.featuresPath
        self._featureSetName = paths.featureSetName
        cmd = "add-featureset {} {} {} --referenceSetName={}".format(
            self._repoPath, self._datasetName, featuresPath,
            self._referenceSetName)
        self.runCommand(cmd)

    def getFeatureSet(self):
        repo = self.readRepo()
        dataset = repo.getDatasetByName(self._datasetName)
        featureSet = dataset.getFeatureSetByName(self._featureSetName)
        return featureSet


class TestAddFeatureSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddFeatureSet, self).setUp()
        self.init()
        self.addDataset()
        self.addOntology()
        self.addReferenceSet()
        self.addFeatureSet()

    def testAddFeatureSet(self):
        featureSet = self.getFeatureSet()
        self.assertEqual(featureSet.getLocalId(), self._featureSetName)
        self.assertEqual(
            featureSet._parentContainer.getLocalId(), self._datasetName)
        self.assertEqual(
            featureSet.getReferenceSet().getLocalId(),
            self._referenceSetName)
        # TODO not clear these fields get populated now
        # self.assertEqual(featureSet.getInfo(), "TODO")
        # self.assertEqual(featureSet.getSourceUrl(), "TODO")


class TestRemoveFeatureSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveFeatureSet, self).setUp()
        self.init()
        self.addDataset()
        self.addOntology()
        self.addReferenceSet()
        self.addFeatureSet()

    def testRemoveFeatureSet(self):
        featureSet = self.getFeatureSet()
        cmd = "remove-featureset {} {} {} -f".format(
            self._repoPath, self._datasetName, featureSet.getLocalId())
        self.runCommand(cmd)
        with self.assertRaises(exceptions.FeatureSetNameNotFoundException):
            self.getFeatureSet()


class TestAddDataset(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddDataset, self).setUp()
        self.init()

    def testDefaults(self):
        name = "test_dataset"
        self.runCommand("add-dataset {} {}".format(self._repoPath, name))
        repo = self.readRepo()
        dataset = repo.getDatasetByName(name)
        self.assertEqual(dataset.getLocalId(), name)

    def testSameName(self):
        name = "test_dataset"
        cmd = "add-dataset {} {}".format(self._repoPath, name)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)


class TestAddReferenceSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddReferenceSet, self).setUp()
        self.init()

    def testDefaults(self):
        fastaFile = paths.ncbi37FaPath
        name = os.path.split(fastaFile)[1].split(".")[0]
        self.runCommand("add-referenceset {} {}".format(
            self._repoPath, fastaFile))
        repo = self.readRepo()
        referenceSet = repo.getReferenceSetByName(name)
        self.assertEqual(referenceSet.getLocalId(), name)
        self.assertEqual(referenceSet.getDataUrl(), fastaFile)
        # TODO check that the default values for all fields are set correctly.

    def testWithName(self):
        name = "test_reference_set"
        fastaFile = paths.ncbi37FaPath
        cmd = "add-referenceset {} {} --name={}".format(
            self._repoPath, fastaFile, name)
        self.runCommand(cmd)
        repo = self.readRepo()
        referenceSet = repo.getReferenceSetByName(name)
        self.assertEqual(referenceSet.getLocalId(), name)
        self.assertEqual(referenceSet.getDataUrl(), fastaFile)

    def testWithSameName(self):
        fastaFile = paths.ncbi37FaPath
        # Default name
        cmd = "add-referenceset {} {}".format(self._repoPath, fastaFile)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.RepoManagerException, self.runCommand, cmd)
        # Specified name
        cmd = "add-referenceset {} {} --name=testname".format(
            self._repoPath, fastaFile)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)


class TestAddOntology(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddOntology, self).setUp()
        self.init()

    def testDefaults(self):
        mapFile = paths.ontologyPath
        name = os.path.split(mapFile)[1].split(".")[0]
        self.runCommand("add-ontology {} {}".format(self._repoPath, mapFile))
        repo = self.readRepo()
        ontologyTermMap = repo.getOntologyTermMapByName(name)
        self.assertEqual(ontologyTermMap.getLocalId(), name)
        self.assertEqual(ontologyTermMap.getDataUrl(), mapFile)

    def testWithName(self):
        mapFile = paths.ontologyPath
        name = "test_name"
        self.runCommand("add-ontology {} {} --name={}".format(
            self._repoPath, mapFile, name))
        repo = self.readRepo()
        ontologyTermMap = repo.getOntologyTermMapByName(name)
        self.assertEqual(ontologyTermMap.getLocalId(), name)
        self.assertEqual(ontologyTermMap.getDataUrl(), mapFile)

    def testWithSameName(self):
        mapFile = paths.ontologyPath
        # Default name
        cmd = "add-ontology {} {}".format(self._repoPath, mapFile)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.RepoManagerException, self.runCommand, cmd)
        # Specified name
        cmd = "add-ontology {} {} --name=testname".format(
            self._repoPath, mapFile)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)


class TestRemoveDataset(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveDataset, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()

    def assertDatasetRemoved(self):
        repo = self.readRepo()
        self.assertRaises(
            exceptions.DatasetNameNotFoundException,
            repo.getDatasetByName, self._datasetName)

    def testEmptyDatasetForce(self):
        self.runCommand("remove-dataset {} {} -f".format(
            self._repoPath, self._datasetName))
        self.assertDatasetRemoved()

    def testContainsReadGroupSet(self):
        self.addReadGroupSet()
        self.runCommand("remove-dataset {} {} -f".format(
            self._repoPath, self._datasetName))
        self.assertDatasetRemoved()


class TestRemoveReadGroupSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveReadGroupSet, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()
        self.addReadGroupSet()

    def assertReadGroupSetRemoved(self):
        repo = self.readRepo()
        dataset = repo.getDatasetByName(self._datasetName)
        self.assertRaises(
            exceptions.ReadGroupSetNameNotFoundException,
            dataset.getReadGroupSetByName, self._readGroupSetName)

    def testWithForce(self):
        self.runCommand("remove-readgroupset {} {} {} -f".format(
            self._repoPath, self._datasetName, self._readGroupSetName))
        self.assertReadGroupSetRemoved()


class TestRemoveReferenceSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveReferenceSet, self).setUp()
        self.init()
        self.addReferenceSet()

    def assertReferenceSetRemoved(self):
        repo = self.readRepo()
        self.assertRaises(
            exceptions.ReferenceSetNameNotFoundException,
            repo.getReferenceSetByName, self._referenceSetName)

    def testDefaults(self):
        self.runCommand("remove-referenceset {} {} -f".format(
            self._repoPath, self._referenceSetName))
        self.assertReferenceSetRemoved()


class TestVerify(AbstractRepoManagerTest):

    def setUp(self):
        super(TestVerify, self).setUp()

    @unittest.skip("Skip test until repo manager is updated")
    def testVerify(self):
        self.init()
        self.addDataset()
        self.addReferenceSet()
        self.addReadGroupSet()
        self.addFeatureSet()
        # TODO fill out with other objects
        cmd = "verify {}".format(self._repoPath)
        self.runCommand(cmd)


class TestRemoveOntology(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveOntology, self).setUp()
        self.init()
        self.addOntology()

    def assertOntologyRemoved(self):
        repo = self.readRepo()
        self.assertRaises(
            exceptions.OntologyNameNotFoundException,
            repo.getOntologyTermMapByName, self._ontologyName)

    def testDefaults(self):
        self.runCommand("remove-ontology {} {} -f".format(
            self._repoPath, self._ontologyName))
        self.assertOntologyRemoved()


class TestAddReadGroupSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddReadGroupSet, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()

    def verifyReadGroupSet(self, name, dataUrl, indexFile):
        repo = self.readRepo()
        dataset = repo.getDatasetByName(self._datasetName)
        referenceSet = repo.getReferenceSetByName(self._referenceSetName)
        readGroupSet = dataset.getReadGroupSetByName(name)
        self.assertEqual(readGroupSet.getLocalId(), name)
        self.assertEqual(readGroupSet.getReferenceSet(), referenceSet)
        self.assertEqual(readGroupSet.getDataUrl(), dataUrl)
        self.assertEqual(readGroupSet.getIndexFile(), indexFile)

    def testDefaultsLocalFile(self):
        bamFile = paths.bamPath
        name = os.path.split(bamFile)[1].split(".")[0]
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, bamFile,
                self._referenceSetName)
        self.runCommand(cmd)
        self.verifyReadGroupSet(name, bamFile, bamFile + ".bai")

    def testLocalFileWithIndex(self):
        bamFile = paths.bamPath
        name = os.path.split(bamFile)[1].split(".")[0]
        with tempfile.NamedTemporaryFile() as temp:
            indexFile = temp.name
            shutil.copyfile(bamFile + ".bai", indexFile)
            cmd = "add-readgroupset {} {} {} {} --referenceSetName={}".format(
                    self._repoPath, self._datasetName, bamFile,
                    indexFile, self._referenceSetName)
            self.runCommand(cmd)
            self.verifyReadGroupSet(name, bamFile, indexFile)

    def testLocalFileWithName(self):
        bamFile = paths.bamPath
        name = "test_rgs"
        cmd = (
            "add-readgroupset {} {} {} --referenceSetName={} "
            "--name={}").format(
            self._repoPath, self._datasetName, bamFile,
            self._referenceSetName, name)
        self.runCommand(cmd)
        self.verifyReadGroupSet(name, bamFile, bamFile + ".bai")

    def testAddReadGroupSetWithSameName(self):
        # Default name
        bamFile = paths.bamPath
        name = os.path.split(bamFile)[1].split(".")[0]
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, bamFile,
                self._referenceSetName)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)
        # Specified name
        name = "test_rgs"
        cmd = (
            "add-readgroupset {} {} {} --referenceSetName={} "
            "--name={}").format(
            self._repoPath, self._datasetName, bamFile,
            self._referenceSetName, name)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)

    def testUrlWithMissingIndex(self):
        bamFile = "http://example.com/example.bam"
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, bamFile,
                self._referenceSetName)
        self.assertRaises(
            exceptions.MissingIndexException, self.runCommand, cmd)

    def testMissingDataset(self):
        bamFile = paths.bamPath
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, "not_a_dataset_name", bamFile,
                self._referenceSetName)
        self.assertRaises(
            exceptions.DatasetNameNotFoundException, self.runCommand, cmd)

    def testMissingReferenceSet(self):
        bamFile = paths.bamPath
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, bamFile,
                "not_a_referenceset_name")
        self.assertRaises(
            exceptions.ReferenceSetNameNotFoundException, self.runCommand, cmd)
