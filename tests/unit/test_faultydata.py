"""
Unit tests for faulty data sets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import unittest

import ga4gh.datamodel.variants as variants
import ga4gh.exceptions as exceptions


class FaultyDataTest(unittest.TestCase):
    def setUp(self):
        self.testDataDir = "tests/faultydata/variants"

    def getFullPath(self, variantSetId):
        return os.path.join(self.testDataDir, variantSetId)


class TestVariantSetNoIndexedVcf(FaultyDataTest):
    setIds = ["no_indexed_vcf"]

    def testInstantiation(self):
        for variantSetId in self.setIds:
            path = self.getFullPath(variantSetId)
            self.assertRaises(
                exceptions.NotIndexedException,
                variants.HtslibVariantSet, variantSetId, path)


class TestInconsistentMetaData(FaultyDataTest):
    setIds = ["inconsist_meta"]

    def testInstantiation(self):
        for variantSetId in self.setIds:
            path = self.getFullPath(variantSetId)
            self.assertRaises(
                exceptions.InconsistentMetaDataException,
                variants.HtslibVariantSet, variantSetId, path)


class TestInconsistentCallSetId(FaultyDataTest):
    setIds = ["inconsist_sampleid", "inconsist_sampleid2"]

    def testInstantiation(self):
        for variantSetId in self.setIds:
            path = self.getFullPath(variantSetId)
            self.assertRaises(
                exceptions.InconsistentCallSetIdException,
                variants.HtslibVariantSet, variantSetId, path)


class TestOverlappingVcfVariants(FaultyDataTest):
    setIds = ["overlapping_vcf"]

    def testInstantiation(self):
        for variantSetId in self.setIds:
            path = self.getFullPath(variantSetId)
            self.assertRaises(
                exceptions.OverlappingVcfException,
                variants.HtslibVariantSet, variantSetId, path)


class TestEmptyDirException(FaultyDataTest):
    setIds = ["empty_dir"]

    def testInstantiation(self):
        for variantSetId in self.setIds:
            path = self.getFullPath(variantSetId)
            self.assertRaises(
                exceptions.EmptyDirException,
                variants.HtslibVariantSet, variantSetId, path)


class TestDuplicateCallSetId(FaultyDataTest):
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

    def testInstantiation(self):
        """
        for variantSetId in self.setIds:
            path = self.getFullPath(variantSetId)
            self.assertRaises(
                exceptions.DuplicateCallSetIdException,
                variants.HtslibVariantSet, variantSetId, path)
        """
        pass
