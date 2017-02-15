"""
Unit tests for variant objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import hashlib
import unittest

import ga4gh.server.datarepo as datarepo
import ga4gh.server.datamodel.variants as variants
import ga4gh.server.datamodel.datasets as datasets

import tests.paths as paths

import ga4gh.schemas.protocol as protocol


class TestHtslibVariantAnnotationSet(unittest.TestCase):
    """
    Unit tests for the abstract variant set.
    """
    def _createVariantAnnotationSet(self, vcfDir):
        """
        Creates a VariantAnnotationSet from the specified directory of
        VCF files.
        """
        self._variantSetName = "testVariantSet"
        self._repo = datarepo.SqlDataRepository(paths.testDataRepo)
        self._repo.open(datarepo.MODE_READ)
        self._dataset = datasets.Dataset("testDs")
        self._variantSet = variants.HtslibVariantSet(
            self._dataset, self._variantSetName)
        self._variantSet.populateFromDirectory(vcfDir)
        self._variantAnnotationSet = variants.HtslibVariantAnnotationSet(
            self._variantSet, "testVAs")
        self._variantAnnotationSet.setOntology(
            self._repo.getOntologyByName(paths.ontologyName))

    def setUp(self):
        vcfDir = "tests/data/datasets/dataset1/variants/WASH7P_annotation"
        self._createVariantAnnotationSet(vcfDir)

    def testConvertLocation(self):
        loc = protocol.AlleleLocation()
        loc.start = 150
        pos = "151/305"
        testLoc = self._variantAnnotationSet.convertLocation(pos)
        self.assertEqual(testLoc, loc)

    def testThousandGenomesAnnotation(self):
        vcfDir = "tests/data/datasets/dataset1/variants/1kg.3.annotations"
        self._createVariantAnnotationSet(vcfDir)
        self.assertTrue(self._variantSet.isAnnotated())

    def testConvertLocationHgvsC(self):
        loc = protocol.AlleleLocation()
        loc.start = 430
        loc.reference_sequence = "T"
        loc.alternate_sequence = "A"
        hgvsC = "NM_001005484.1:c.431T>A"
        testLoc = self._variantAnnotationSet.convertLocationHgvsC(hgvsC)
        self.assertEqual(testLoc, loc)

    def testConvertLocationHgvsP(self):
        loc = protocol.AlleleLocation()
        loc.start = 143
        loc.alternate_sequence = "Asn"
        loc.reference_sequence = "Ile"
        hgvsP = "NM_001005484.1:p.Ile144Asn"
        testLoc = self._variantAnnotationSet.convertLocationHgvsP(hgvsP)
        self.assertEqual(testLoc, loc)

    def testAddLocations(self):
        effect = protocol.TranscriptEffect()
        effect.hgvs_annotation.protein = "NM_001005484.1:p.Ile144Asn"
        effect.hgvs_annotation.transcript = "NM_001005484.1:c.431T>A"
        effect.protein_location.alternate_sequence = "Asn"
        effect.protein_location.reference_sequence = "Ile"
        effect.protein_location.start = 143
        effect.cds_location.alternate_sequence = "A"
        effect.cds_location.reference_sequence = "T"
        effect.cds_location.start = 430
        effect.cdna_location.start = 430
        protPos = "144/305"
        cdnaPos = "431/918"
        testEffect = self._variantAnnotationSet.addLocations(
            effect, protPos, cdnaPos)
        self.assertEqual(testEffect, effect)

    def testHashVariantAnnotation(self):
        annotation = protocol.VariantAnnotation()
        variant = protocol.Variant()
        expected = hashlib.md5('\t()\t[]\t').hexdigest()
        hashed = self._variantAnnotationSet.hashVariantAnnotation(
            variant, annotation)
        self.assertEqual(hashed, expected)

    def testGetTranscriptEffectId(self):
        effect = protocol.TranscriptEffect()
        expected = hashlib.md5("\t\t[]\t").hexdigest()
        hashed = self._variantAnnotationSet.getTranscriptEffectId(effect)
        self.assertEqual(hashed, expected)
