"""
Unit tests for faulty data sets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import glob
import os
import unittest

import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol
import ga4gh.datasource.htslib as htslib


def getVariantSetFromDirectory(name, directory):
    """
    Returns a variant set populated with the VCF files in the specified
    directory. Convenience method.
    """
    pattern = os.path.join(directory, "*.vcf.gz")
    files = glob.glob(pattern)
    indexes = [vcf + ".tbi" for vcf in files]
    return htslib.HtslibVariantSet(name, files, indexes)


class FaultyVariantDataTest(unittest.TestCase):
    """
    Superclass of faulty variant data tests.
    """
    def setUp(self):
        self.testDataDir = "tests/faultydata/variants"
        # self.dataset = datasets.Dataset('dataset1')

    def getFullPath(self, localId):
        return os.path.join(self.testDataDir, localId)


class TestVariantSetNoIndexedVcf(FaultyVariantDataTest):
    localIds = ["no_indexed_vcf"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.NotIndexedException):
                getVariantSetFromDirectory(localId, path)


@unittest.skip("FIX VCF import consistencty check")
class TestInconsistentMetaData(FaultyVariantDataTest):
    localIds = ["inconsist_meta"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.InconsistentMetaDataException):
                getVariantSetFromDirectory(localId, path)


@unittest.skip("FIX VCF import consistencty check")
class TestInconsistentCallSetId(FaultyVariantDataTest):
    localIds = ["inconsist_sampleid", "inconsist_sampleid2"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.InconsistentCallSetIdException):
                getVariantSetFromDirectory(localId, path)


class TestOverlappingVcfVariants(FaultyVariantDataTest):
    localIds = ["overlapping_vcf"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.OverlappingVcfException):
                getVariantSetFromDirectory(localId, path)


class TestEmptyDirException(FaultyVariantDataTest):
    localIds = ["empty_dir"]

    def testInstantiation(self):
        for localId in self.localIds:
            path = self.getFullPath(localId)
            with self.assertRaises(exceptions.EmptyVariantSetException):
                getVariantSetFromDirectory(localId, path)


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
            with self.assertRaises(exceptions.DuplicateCallSetIdException):
                getVariantSetFromDirectory(localId, path)
