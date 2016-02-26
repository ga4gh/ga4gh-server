"""
Tests for the repo manager tool
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import tempfile
import unittest

import ga4gh.repo_manager as repo_manager
import ga4gh.exceptions as exceptions
import tests.paths as paths


class AbstractRepoManagerTest(unittest.TestCase):
    """
    Base class for repo manager tests
    """
    moveMode = 'link'

    def getTempDirPath(self):
        # tempfile doesn't provide methods for just getting a standalone
        # dir name, so we need to create then delete the dir
        tempdir = tempfile.mkdtemp()
        os.rmdir(tempdir)
        return tempdir


class RepoManagerTest(AbstractRepoManagerTest):
    """
    End to end test for the repo manager
    """
    def setUp(self):
        self.testDirName = self.getTempDirPath()
        self.repoManager = repo_manager.RepoManager(self.testDirName)

    def tearDown(self):
        self.repoManager.destroy()

    def testEndtoEnd(self):
        self.repoManager.init()
        datasetName = 'datasetOne'
        self.repoManager.addDataset(datasetName)
        metadata = {'description': 'aDescription'}
        self.repoManager.addReferenceSet(
            paths.faPath, self.moveMode, metadata)
        self.repoManager.addReadGroupSet(
            datasetName, paths.bamPath, self.moveMode)
        self.repoManager.addVariantSet(
            datasetName, paths.vcfDirPath, self.moveMode)
        with self.assertRaises(exceptions.ReferenceSetNameNotFoundException):
            # ReferenceSet named 'Default' does not exist
            self.repoManager.check(doConsistencyCheck=True)
        self.repoManager.check()
        self.repoManager.list()
        self.repoManager.removeReadGroupSet(
            datasetName, paths.readGroupSetName)
        self.repoManager.removeVariantSet(
            datasetName, paths.variantSetName)
        self.repoManager.removeReferenceSet(paths.referenceSetName)
        self.repoManager.removeDataset(datasetName)


class RepoManagerInidividualCommandTest(AbstractRepoManagerTest):
    """
    Tests for individiual repo manager commands
    """
    def setUp(self):
        self.tempdir = self.getTempDirPath()
        self.repoManager = repo_manager.RepoManager(self.tempdir)
        self.repoManager.init()

    def tearDown(self):
        try:
            self.repoManager.destroy()
        except exceptions.RepoManagerException:
            pass

    def testInit(self):
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.init()

    def testDestroy(self):
        self.repoManager.destroy()
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.destroy()

    def testCheck(self):
        self.repoManager.check()

    def testList(self):
        self.repoManager.list()

    def testAddDataset(self):
        self.repoManager.addDataset('dataset1')
        self.repoManager.addDataset('dataset2')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addDataset('dataset2')

    def testRemoveDataset(self):
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeDataset('dataset1')
        self.repoManager.addDataset('dataset1')
        self.repoManager.removeDataset('dataset1')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeDataset('dataset1')

    def testAddReferenceSet(self):
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addReferenceSet(paths.bamPath, 'link', {})
        self.repoManager.addReferenceSet(
            paths.faPath, 'link', {'description': 'aDescription'})
        self.repoManager.addReferenceSet(
            paths.faPath2, 'copy', {})
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addReferenceSet(paths.faPath, 'link', {})

    def testRemoveReferenceSet(self):
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeReferenceSet(paths.referenceSetName)
        self.repoManager.addReferenceSet(paths.faPath, 'link', {})
        self.repoManager.removeReferenceSet(paths.referenceSetName)
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeReferenceSet(paths.referenceSetName)

    def testAddReadGroupSet(self):
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addReadGroupSet(
                'dataset1', paths.bamPath, 'link')
        self.repoManager.addDataset('dataset1')
        self.repoManager.addDataset('dataset2')
        self.repoManager.addReadGroupSet(
            'dataset1', paths.bamPath, 'link')
        self.repoManager.addReadGroupSet(
            'dataset1', paths.bamPath2, 'link')
        self.repoManager.addReadGroupSet(
            'dataset2', paths.bamPath, 'link')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addReadGroupSet(
                'dataset1', paths.bamPath, 'link')

    def testRemoveReadGroupSet(self):
        self.repoManager.addDataset('dataset1')
        self.repoManager.addReadGroupSet(
            'dataset1', paths.bamPath, 'link')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeReadGroupSet(
                'dataset2', paths.readGroupSetName)
        self.repoManager.removeReadGroupSet(
            'dataset1', paths.readGroupSetName)
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeReadGroupSet(
                'dataset1', paths.readGroupSetName)

    def testAddVariantSet(self):
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addVariantSet(
                'dataset1', paths.vcfDirPath, 'link')
        self.repoManager.addDataset('dataset1')
        self.repoManager.addDataset('dataset2')
        self.repoManager.addVariantSet(
            'dataset1', paths.vcfDirPath, 'link')
        self.repoManager.addVariantSet(
            'dataset1', paths.vcfDirPath2, 'link')
        self.repoManager.addVariantSet(
            'dataset2', paths.vcfDirPath, 'link')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addVariantSet(
                'dataset1', paths.vcfDirPath, 'link')

    def testRemoveVariantSet(self):
        self.repoManager.addDataset('dataset1')
        self.repoManager.addVariantSet(
            'dataset1', paths.vcfDirPath, 'link')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeVariantSet(
                'dataset2', paths.variantSetName)
        self.repoManager.removeVariantSet(
            'dataset1', paths.variantSetName)
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeVariantSet(
                'dataset1', paths.variantSetName)
