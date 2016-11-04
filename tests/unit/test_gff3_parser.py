"""
GFF3 parser unit tests.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.server.gff3 as gff3
import unittest

_testDataDir = "tests/data/datasets/dataset1/sequenceAnnotations/"


class TestGff3ParserOnTypicalFile(unittest.TestCase):
    """
    Data driven unit tests for the GFF3 parser
    """
    def setUp(self):
        testDataFile = _testDataDir + "gencodeV21Set1.gff3"
        self.gff3Parser = gff3.Gff3Parser(testDataFile)
        self.gff3Data = self.gff3Parser.parse()

    def testFileParsedHasSomeRootFeatures(self):
        self.assertIsNotNone(self.gff3Data.roots, "No root features")
        self.assertNotEqual(len(self.gff3Data.roots), 0, "No root features")

    def testSomeFeatureIsWellFormed(self):
        featId = self.gff3Data.byFeatureName.keys()[0]
        feat = self.gff3Data.byFeatureName[featId][0]
        self.assertEqual(featId, feat.featureName, "featureName mismatch")
        self.assertIsNotNone(feat.seqname, "sequence name is not populated")
        self.assertGreaterEqual(feat.end, feat.start, "end less than start")
        self.assertIn(feat.strand, u"+-", "strand is neither + nor -")
        self.assertIsNotNone(feat.source, "source is unspecified")
        self.assertIsNotNone(feat.type, "feature type is unspecified")
        self.assertIsInstance(feat.parents, set, "parents not a set")
        self.assertIsInstance(feat.children, set, "children not a set")

    def testRootFeaturesHaveNoParents(self):
        for root in self.gff3Data.roots:
            self.assertEqual(
                len(root.parents), 0, "root feature has a parent")

    def testAllFeaturesContainAllRootFeatures(self):
        for root in self.gff3Data.roots:
            feat = self.gff3Data.byFeatureName[root.featureName]
            self.assertGreaterEqual(
                len(feat), 1,
                "root feature not in list of all features")

    def testInvalidFeatureNameKeyQueryFails(self):
        badFeatureName = "987654"
        badFeat = self.gff3Data.byFeatureName[badFeatureName]
        self.assertEqual(
            len(badFeat), 0,
            "invalid feature ID returned valid object")

    def testAllChildrenFeaturesArePresentInSet(self):
        for featList in self.gff3Data.byFeatureName.values():
            for feat in featList:
                for child in feat.children:
                    childLookup = self.gff3Data.byFeatureName[
                        child.featureName]
                    self.assertGreaterEqual(
                        len(childLookup), 1,
                        "child feature not in set")


class TestGff3ParserOnDiscontinuousFeatureFile(TestGff3ParserOnTypicalFile):
    """
    Data driven parser test on file with discontinuous features.
    The tests here rely on specific data in the file being parsed.
    """
    def setUp(self):
        testDataFile = _testDataDir + "discontinuous.gff3"
        self.gff3Parser = gff3.Gff3Parser(testDataFile)
        self.gff3Data = self.gff3Parser.parse()

    def testDiscontinuousFeature(self):
        feat = self.gff3Data.byFeatureName['apidb|cds_MAL13P1.103-1']
        self.assertEqual(
            len(feat), 10,
            "not all parts of discontinuous feature parsed")


class TestGff3ParserOnSacCerFile(TestGff3ParserOnTypicalFile):
    """
    Data driven parser test on file from Saccharomyces cerevisiae S288C genome.
    """
    def setUp(self):
        testDataFile = _testDataDir + "sacCerTest.gff3"
        self.gff3Parser = gff3.Gff3Parser(testDataFile)
        self.gff3Data = self.gff3Parser.parse()


class TestGff3ParserOnSpecialCasesFile(TestGff3ParserOnTypicalFile):
    """
    Data driven parser test on a GFF3 file representing edge cases.
    """
    def setUp(self):
        testDataFile = _testDataDir + "specialCasesTest.gff3"
        self.gff3Parser = gff3.Gff3Parser(testDataFile)
        self.gff3Data = self.gff3Parser.parse()
