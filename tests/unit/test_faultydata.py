"""
Unit tests for faulty data sets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import unittest

import ga4gh.datamodel.variants as variants
import ga4gh.datamodel.references as references
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol


class FaultyVariantDataTest(unittest.TestCase):
    def setUp(self):
        self.testDataDir = "tests/faultydata/variants"

    def getFullPath(self, variant_set_id):
        return os.path.join(self.testDataDir, variant_set_id)


class TestVariantSetNoIndexedVcf(FaultyVariantDataTest):
    setIds = ["no_indexed_vcf"]

    def testInstantiation(self):
        for variant_set_id in self.setIds:
            path = self.getFullPath(variant_set_id)
            self.assertRaises(
                exceptions.NotIndexedException,
                variants.HtslibVariantSet, variant_set_id, path)


class TestInconsistentMetaData(FaultyVariantDataTest):
    setIds = ["inconsist_meta"]

    def testInstantiation(self):
        for variant_set_id in self.setIds:
            path = self.getFullPath(variant_set_id)
            self.assertRaises(
                exceptions.InconsistentMetaDataException,
                variants.HtslibVariantSet, variant_set_id, path)


class TestInconsistentCallSetId(FaultyVariantDataTest):
    setIds = ["inconsist_sampleid", "inconsist_sampleid2"]

    def testInstantiation(self):
        for variant_set_id in self.setIds:
            path = self.getFullPath(variant_set_id)
            self.assertRaises(
                exceptions.InconsistentCallSetIdException,
                variants.HtslibVariantSet, variant_set_id, path)


class TestOverlappingVcfVariants(FaultyVariantDataTest):
    setIds = ["overlapping_vcf"]

    def testInstantiation(self):
        for variant_set_id in self.setIds:
            path = self.getFullPath(variant_set_id)
            self.assertRaises(
                exceptions.OverlappingVcfException,
                variants.HtslibVariantSet, variant_set_id, path)


class TestEmptyDirException(FaultyVariantDataTest):
    setIds = ["empty_dir"]

    def testInstantiation(self):
        for variant_set_id in self.setIds:
            path = self.getFullPath(variant_set_id)
            self.assertRaises(
                exceptions.EmptyDirException,
                variants.HtslibVariantSet, variant_set_id, path)


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
    setIds = ["duplicated_sampleid"]

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testInstantiation(self):
        for variant_set_id in self.setIds:
            path = self.getFullPath(variant_set_id)
            self.assertRaises(
                exceptions.DuplicateCallSetIdException,
                variants.HtslibVariantSet, variant_set_id, path)


class FaultyReferenceDataTest(unittest.TestCase):
    def setUp(self):
        self.testDataDir = "tests/faultydata/references"

    def getFullPath(self, reference_set_id):
        return os.path.join(self.testDataDir, reference_set_id)


class TestTwoReferences(FaultyReferenceDataTest):
    setIds = ["two_references"]

    def testInstantiation(self):
        for reference_set_id in self.setIds:
            path = self.getFullPath(reference_set_id)
            with self.assertRaises(exceptions.NotExactlyOneReferenceException):
                references.HtslibReferenceSet(reference_set_id, path)
