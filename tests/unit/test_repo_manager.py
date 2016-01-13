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


class AbstractRepoManagerTest(unittest.TestCase):
    """
    Base class for repo manager tests
    """
    vcfDirPath = os.path.abspath(
        'tests/data/datasets/dataset1/variants/1kgPhase1')
    vcfDirPath2 = os.path.abspath(
        'tests/data/datasets/dataset1/variants/1kgPhase3')
    faPath = os.path.abspath(
        'tests/data/referenceSets/Default/chr17.fa.gz')
    faPath2 = os.path.abspath(
        'tests/data/referenceSets/example_1/simple.fa.gz')
    faPath3 = os.path.abspath(
        'tests/data/referenceSets/example_2/random1.fa.gz')
    bamPath = os.path.abspath(
        'tests/data/datasets/dataset1/reads/chr17.1-250.bam')
    bamPath2 = os.path.abspath(
        'tests/data/datasets/dataset1/reads/'
        'wgEncodeUwRepliSeqBg02esG1bAlnRep1_sample.bam')
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
            self.faPath, self.moveMode, metadata)
        self.repoManager.addReadGroupSet(
            datasetName, self.bamPath, self.moveMode)
        self.repoManager.addVariantSet(
            datasetName, self.vcfDirPath, self.moveMode)
        self.repoManager.check()
        self.repoManager.list()
        self.repoManager.removeReadGroupSet(datasetName, 'chr17.1-250')
        self.repoManager.removeVariantSet(datasetName, '1kgPhase1')
        self.repoManager.removeReferenceSet('chr17')
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
            self.repoManager.addReferenceSet(self.bamPath, 'link', {})
        self.repoManager.addReferenceSet(
            self.faPath, 'link', {'description': 'aDescription'})
        self.repoManager.addReferenceSet(
            self.faPath2, 'copy', {})
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addReferenceSet(self.faPath, 'link', {})

    def testRemoveReferenceSet(self):
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeReferenceSet('chr17')
        self.repoManager.addReferenceSet(self.faPath, 'link', {})
        self.repoManager.removeReferenceSet('chr17')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeReferenceSet('chr17')

    def testAddReadGroupSet(self):
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addReadGroupSet(
                'dataset1', self.bamPath, 'link')
        self.repoManager.addDataset('dataset1')
        self.repoManager.addDataset('dataset2')
        self.repoManager.addReadGroupSet(
            'dataset1', self.bamPath, 'link')
        self.repoManager.addReadGroupSet(
            'dataset1', self.bamPath2, 'link')
        self.repoManager.addReadGroupSet(
            'dataset2', self.bamPath, 'link')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addReadGroupSet(
                'dataset1', self.bamPath, 'link')

    def testRemoveReadGroupSet(self):
        self.repoManager.addDataset('dataset1')
        self.repoManager.addReadGroupSet(
            'dataset1', self.bamPath, 'link')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeReadGroupSet('dataset2', 'chr17.1-250')
        self.repoManager.removeReadGroupSet('dataset1', 'chr17.1-250')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeReadGroupSet('dataset1', 'chr17.1-250')

    def testAddVariantSet(self):
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addVariantSet(
                'dataset1', self.vcfDirPath, 'link')
        self.repoManager.addDataset('dataset1')
        self.repoManager.addDataset('dataset2')
        self.repoManager.addVariantSet(
            'dataset1', self.vcfDirPath, 'link')
        self.repoManager.addVariantSet(
            'dataset1', self.vcfDirPath2, 'link')
        self.repoManager.addVariantSet(
            'dataset2', self.vcfDirPath, 'link')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.addVariantSet(
                'dataset1', self.vcfDirPath, 'link')

    def testRemoveVariantSet(self):
        self.repoManager.addDataset('dataset1')
        self.repoManager.addVariantSet(
            'dataset1', self.vcfDirPath, 'link')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeVariantSet('dataset2', '1kgPhase1')
        self.repoManager.removeVariantSet('dataset1', '1kgPhase1')
        with self.assertRaises(exceptions.RepoManagerException):
            self.repoManager.removeVariantSet('dataset1', '1kgPhase1')
