"""
Tests for the repo manager tool
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import tempfile
import unittest

import ga4gh.exceptions as exceptions
import ga4gh.datarepo as datarepo
import ga4gh.cli as cli
# import tests.paths as paths


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
        repo.open("r")
        return repo


class TestAddDataset(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddDataset, self).setUp()
        self.runCommand("init {}".format(self._repoPath))

    def testAddDataset(self):
        name = "test_dataset"
        self.runCommand("add-dataset {} {}".format(self._repoPath, name))
        repo = self.readRepo()
        dataset = repo.getDatasetByName(name)
        self.assertEqual(dataset.getLocalId(), name)

    def testAddDatasetWithSameName(self):
        name = "test_dataset"
        cmd = "add-dataset {} {}".format(self._repoPath, name)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.RepoManagerException, self.runCommand, cmd)
    # TODO adapt the code below to test the repo manager code within
    # the above framework.
# class RepoManagerInidividualCommandTest(AbstractRepoManagerTest):
    # """
    # Tests for individiual repo manager commands
    # """
    # def setUp(self):
    #     self.tempdir = self.getTempDirPath()
    #     # self.repoManager = repo_manager.RepoManager(self.tempdir)
    #     self.repoManager.init()

    # def tearDown(self):
    #     try:
    #         self.repoManager.destroy()
    #     except exceptions.RepoManagerException:
    #         pass

    # @unittest.skip("Skip until repo manager completed")
    # def testInit(self):
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.init()

    # @unittest.skip("Skip until repo manager completed")
    # def testDestroy(self):
    #     self.repoManager.destroy()
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.destroy()

    # @unittest.skip("Skip until repo manager completed")
    # def testCheck(self):
    #     self.repoManager.check()

    # @unittest.skip("Skip until repo manager completed")
    # def testList(self):
    #     self.repoManager.list()

    # @unittest.skip("Skip until repo manager completed")
    # def testAddDataset(self):
    #     self.repoManager.addDataset('dataset1')
    #     self.repoManager.addDataset('dataset2')
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.addDataset('dataset2')

    # @unittest.skip("Skip until repo manager completed")
    # def testRemoveDataset(self):
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeDataset('dataset1')
    #     self.repoManager.addDataset('dataset1')
    #     self.repoManager.removeDataset('dataset1')
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeDataset('dataset1')

    # @unittest.skip("Skip until repo manager completed")
    # def testAddReferenceSet(self):
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.addReferenceSet(paths.bamPath, 'link', {})
    #     self.repoManager.addReferenceSet(
    #         paths.faPath, 'link', {'description': 'aDescription'})
    #     self.repoManager.addReferenceSet(
    #         paths.faPath2, 'copy', {})
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.addReferenceSet(paths.faPath, 'link', {})

    # @unittest.skip("Skip until repo manager completed")
    # def testRemoveReferenceSet(self):
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeReferenceSet(paths.referenceSetName)
    #     self.repoManager.addReferenceSet(paths.faPath, 'link', {})
    #     self.repoManager.removeReferenceSet(paths.referenceSetName)
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeReferenceSet(paths.referenceSetName)

    # @unittest.skip("Skip until repo manager completed")
    # def testAddReadGroupSet(self):
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.addReadGroupSet(
    #             'dataset1', paths.bamPath, 'link')
    #     self.repoManager.addDataset('dataset1')
    #     self.repoManager.addDataset('dataset2')
    #     self.repoManager.addReadGroupSet(
    #         'dataset1', paths.bamPath, 'link')
    #     self.repoManager.addReadGroupSet(
    #         'dataset1', paths.bamPath2, 'link')
    #     self.repoManager.addReadGroupSet(
    #         'dataset2', paths.bamPath, 'link')
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.addReadGroupSet(
    #             'dataset1', paths.bamPath, 'link')

    # @unittest.skip("Skip until repo manager completed")
    # def testRemoveReadGroupSet(self):
    #     self.repoManager.addDataset('dataset1')
    #     self.repoManager.addReadGroupSet(
    #         'dataset1', paths.bamPath, 'link')
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeReadGroupSet(
    #             'dataset2', paths.readGroupSetName)
    #     self.repoManager.removeReadGroupSet(
    #         'dataset1', paths.readGroupSetName)
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeReadGroupSet(
    #             'dataset1', paths.readGroupSetName)

    # @unittest.skip("Skip until repo manager completed")
    # def testAddVariantSet(self):
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.addVariantSet(
    #             'dataset1', paths.vcfDirPath, 'link')
    #     self.repoManager.addDataset('dataset1')
    #     self.repoManager.addDataset('dataset2')
    #     self.repoManager.addVariantSet(
    #         'dataset1', paths.vcfDirPath, 'link')
    #     self.repoManager.addVariantSet(
    #         'dataset1', paths.vcfDirPath2, 'link')
    #     self.repoManager.addVariantSet(
    #         'dataset2', paths.vcfDirPath, 'link')
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.addVariantSet(
    #             'dataset1', paths.vcfDirPath, 'link')

    # @unittest.skip("Skip until repo manager completed")
    # def testRemoveVariantSet(self):
    #     self.repoManager.addDataset('dataset1')
    #     self.repoManager.addVariantSet(
    #         'dataset1', paths.vcfDirPath, 'link')
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeVariantSet(
    #             'dataset2', paths.variantSetName)
    #     self.repoManager.removeVariantSet(
    #         'dataset1', paths.variantSetName)
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeVariantSet(
    #             'dataset1', paths.variantSetName)

    # @unittest.skip("Skip until repo manager completed")
    # def testAddOntologyMap(self):
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.addOntologyMap(paths.bamPath, 'link')
    #     self.repoManager.addOntologyMap(
    #         paths.ontologyPath, 'link')
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.addOntologyMap(
    #             paths.ontologyPath, 'link')

    # @unittest.skip("Skip until repo manager completed")
    # def testRemoveOntologyMap(self):
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeOntologyMap(
    #             paths.ontologyName)
    #     self.repoManager.addOntologyMap(
    #         paths.ontologyPath, 'link')
    #     self.repoManager.removeOntologyMap(
    #         paths.ontologyName)
    #     with self.assertRaises(exceptions.RepoManagerException):
    #         self.repoManager.removeOntologyMap(
    #             paths.ontologyName)
