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
import ga4gh.backend as backend


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
                variants.HtslibVariantSet, self.dataset, localId, path,
                None)


class TestInconsistentMetaData(FaultyVariantDataTest):
    localIds = ["inconsist_meta"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.InconsistentMetaDataException):
                variants.HtslibVariantSet(self.dataset, localId, path, None)


class TestInconsistentCallSetId(FaultyVariantDataTest):
    localIds = ["inconsist_sampleid", "inconsist_sampleid2"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.InconsistentCallSetIdException):
                variants.HtslibVariantSet(self.dataset, localId, path, None)


class TestOverlappingVcfVariants(FaultyVariantDataTest):
    localIds = ["overlapping_vcf"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.OverlappingVcfException):
                variants.HtslibVariantSet(self.dataset, localId, path, None)


class TestEmptyDirException(FaultyVariantDataTest):
    localIds = ["empty_dir"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            self.assertRaises(
                exceptions.EmptyDirException,
                variants.HtslibVariantSet, self.dataset, localId, path, None)


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
                variants.HtslibVariantSet, self.dataset, localId, path,
                None)


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
            references.HtslibReferenceSet, localId, path, None)


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
            references.HtslibReferenceSet, localId, path, None)


class FaultyReferenceSetDataTest(unittest.TestCase):
    """
    Superclass of faulty reference set data tests
    """
    def getFullPath(self, localId):
        testDataDir = "tests/faultydata/references"
        return os.path.join(testDataDir, localId)


class TestNoReferenceSetMetadata(FaultyReferenceSetDataTest):
    """
    Tests an error is thrown with a missing reference set metadata file
    """
    def testNoReferenceSetMetadata(self):
        localId = "no_refset_meta"
        path = self.getFullPath(localId)
        with self.assertRaises(IOError):
            references.HtslibReferenceSet(localId, path, None)


class TestMissingReferenceSetMetadata(FaultyReferenceSetDataTest):
    """
    Tests an error is thrown with a reference set metadata file that
    is missing entries
    """
    def testMissingReferenceSetMetadata(self):
        localId = "missing_refset_meta"
        path = self.getFullPath(localId)
        with self.assertRaises(exceptions.MissingReferenceSetMetadata):
            references.HtslibReferenceSet(localId, path, None)


class TestInvalidReferenceSetMetadata(FaultyReferenceSetDataTest):
    """
    Tests an error is thrown with a reference set metadata file that
    can not be parsed
    """
    def testMissingReferenceSetMetadata(self):
        localId = "invalid_refset_meta"
        path = self.getFullPath(localId)
        with self.assertRaises(ValueError):
            references.HtslibReferenceSet(localId, path, None)


class FaultyDatasetTest(unittest.TestCase):
    """
    Superclass of faulty dataset tests.
    """
    def setUp(self):
        self.testDataDir = "tests/faultydata/datasets"

    def getFullPath(self, localId):
        return os.path.join(self.testDataDir, localId)


class TestBadDatasetMetadata(FaultyDatasetTest):
    """
    Tests that we raise an expcetion if the metadata is not correct.
    """
    def testBadReferenceDatasetMetadata(self):
        localId = "bad_metadata"
        path = self.getFullPath(localId)
        localBackend = backend.EmptyBackend()
        with self.assertRaises(exceptions.MissingDatasetMetadataException):
            datasets.FileSystemDataset(localId, path, localBackend)
