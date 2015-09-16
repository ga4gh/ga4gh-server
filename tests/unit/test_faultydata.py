"""
Unit tests for faulty data sets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import unittest

import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.references as references
import ga4gh.datamodel.variants as variants
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol


class FaultyVariantDataTest(unittest.TestCase):
    """
    Superclass of faulty variant data tests.
    """
    def setUp(self):
        self.testDataDir = "tests/faultydata/variants"
        self.dataset = datasets.AbstractDataset('dataset1')

    def getFullPath(self, localId):
        return os.path.join(self.testDataDir, localId)


class TestVariantSetNoIndexedVcf(FaultyVariantDataTest):
    localIds = ["no_indexed_vcf"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            self.assertRaises(
                exceptions.NotIndexedException,
                variants.HtslibVariantSet, self.dataset, localId, path)


class TestInconsistentMetaData(FaultyVariantDataTest):
    localIds = ["inconsist_meta"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.InconsistentMetaDataException):
                variants.HtslibVariantSet(self.dataset, localId, path)


class TestInconsistentCallSetId(FaultyVariantDataTest):
    localIds = ["inconsist_sampleid", "inconsist_sampleid2"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.InconsistentCallSetIdException):
                variants.HtslibVariantSet(self.dataset, localId, path)


class TestOverlappingVcfVariants(FaultyVariantDataTest):
    localIds = ["overlapping_vcf"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.OverlappingVcfException):
                variants.HtslibVariantSet(self.dataset, localId, path)


class TestEmptyDirException(FaultyVariantDataTest):
    localIds = ["empty_dir"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            self.assertRaises(
                exceptions.EmptyDirException,
                variants.HtslibVariantSet, self.dataset, localId, path)


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
            self.assertRaises(
                exceptions.DuplicateCallSetIdException,
                variants.HtslibVariantSet, self.dataset, localId, path)


class FaultyReferenceDataTest(unittest.TestCase):
    """
    Superclass of faulty reference data tests
    """

    def getFullPath(self, localId):
        testDataDir = "tests/faultydata/references"
        return os.path.join(testDataDir, localId)


class TestTwoReferences(FaultyReferenceDataTest):
    """
    Tests for FASTA files with more than one reference.
    """

    def testInstantiation(self):
        localId = "two_references"
        path = self.getFullPath(localId)
        self.assertRaises(
            exceptions.NotExactlyOneReferenceException,
            references.HtslibReferenceSet, localId, path)


class TestInconsistentReferenceName(FaultyReferenceDataTest):
    """
    Tests the case in which we have a reference file with a different
    name to the ID in the fasta file.
    """

    def testInstantiation(self):
        localId = "inconsistent_reference_name"
        path = self.getFullPath(localId)
        self.assertRaises(
            exceptions.InconsistentReferenceNameException,
            references.HtslibReferenceSet, localId, path)
