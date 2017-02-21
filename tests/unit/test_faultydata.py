"""
Unit tests for faulty data sets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import unittest

import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.datamodel.variants as variants
import ga4gh.server.exceptions as exceptions

import ga4gh.schemas.protocol as protocol


class FaultyVariantDataTest(unittest.TestCase):
    """
    Superclass of faulty variant data tests.
    """
    def setUp(self):
        self.testDataDir = "tests/faultydata/variants"
        self.dataset = datasets.Dataset('dataset1')

    def getFullPath(self, localId):
        return os.path.join(self.testDataDir, localId)


class TestVariantSetNoIndexedVcf(FaultyVariantDataTest):
    localIds = ["no_indexed_vcf"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            variantSet = variants.HtslibVariantSet(self.dataset, localId)
            with self.assertRaises(exceptions.NotIndexedException):
                variantSet.populateFromDirectory(path)


class TestInconsistentMetaData(FaultyVariantDataTest):
    localIds = ["inconsist_meta"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            variantSet = variants.HtslibVariantSet(self.dataset, localId)
            variantSet.populateFromDirectory(path)
            with self.assertRaises(exceptions.InconsistentMetaDataException):
                variantSet.checkConsistency()


class TestInconsistentCallSetId(FaultyVariantDataTest):
    localIds = ["inconsist_sampleid", "inconsist_sampleid2"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            variantSet = variants.HtslibVariantSet(self.dataset, localId)
            variantSet.populateFromDirectory(path)
            with self.assertRaises(exceptions.InconsistentCallSetIdException):
                variantSet.checkConsistency()


class TestOverlappingVcfVariants(FaultyVariantDataTest):
    localIds = ["overlapping_vcf"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            variantSet = variants.HtslibVariantSet(self.dataset, localId)
            with self.assertRaises(exceptions.OverlappingVcfException):
                variantSet.populateFromDirectory(path)


class TestDuplicateCallSetId(FaultyVariantDataTest):
    """
    THIS SECTION IS CURRENTLY NOT WORKING
    It returns the following error:

    [E::bcf_hdr_add_sample] Duplicated sample name 'S1'
    Aborted (core dumped)

    which is coming from:
    htslib/vcf.c function bcf_hdr_add_sample

    UNABLE TO CAPTURE EXCEPTION
    """
    localIds = ["duplicated_sampleid"]

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            variantSet = variants.HtslibVariantSet(self.dataset, localId)
            with self.assertRaises(exceptions.DuplicateCallSetIdException):
                variantSet.populateFromDirectory(path)


class FaultyReferenceDataTest(unittest.TestCase):
    """
    Superclass of faulty reference data tests
    """
    def getFullPath(self, localId):
        testDataDir = "tests/faultydata/references"
        return os.path.join(testDataDir, localId)
