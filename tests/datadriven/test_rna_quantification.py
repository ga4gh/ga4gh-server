"""
Data-driven tests for rna quantification.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import random
import ga4gh.protocol as protocol
import ga4gh.datamodel.rna_quantification as rna_quantification
import tests.datadriven as datadriven

# TODO: tests are a bit silly now as it is reading the flat file in the same
# way as the module.  Will have more meaning once rna_quantification is
# switched over to use the parsers


def testReferenceSets():
    testDataDir = "tests/data/dataset1/rnaQuant"
    for test in datadriven.makeTests(testDataDir, RnaQuantificationTest):
        yield test


class RnaQuantificationTest(datadriven.DataDrivenTest):
    """
    Data driven test class for rna quantification. Builds an alternative model
    of a rna quantification, and verifies that it is consistent with the model
    built by the rna_quantification.RNASeqResult object.
    """
    class RnaQuantInfo(object):
        """
        Container class for information about a quantification.
        test data is a tab file with no header.  Columns are:
        Id, annotations, description, name, readGroupId
        where annotation is a comma separated list
        """
        def __init__(self, rnaQuantFileName):
            rnaQuantificationData = open(rnaQuantFileName, "r")
            quantData = rnaQuantificationData.readline()
            fields = quantData.strip().split('\t')
            self.id = fields[0]
            self.annotationIds = fields[1].split(',')
            self.description = fields[2]
            self.name = fields[3]
            self.readGroupId = fields[4]

    class CharacterizationInfo(object):
        """
        Container class for characterization information related to a
        quantification.
        test data is tab file with no header.  Columns are:
        analysisId, complexity, exonicFraction, fractionMapped,
        intergenicFraction, intronicFraction
        """
        def __init__(self, characterizationFileName):
            characterizationData = open(characterizationFileName, "r")
            characterization = characterizationData.readline()
            fields = characterization.strip().split('\t')
            self.analysisId = fields[0]
            self.complexity = float(fields[1])
            self.exonicFraction = float(fields[2])
            self.fractionMapped = float(fields[3])
            self.intergenicFraction = float(fields[4])
            self.intronicFraction = float(fields[5])

    class ReadCountInfo(object):
        """
        Container class for characterization information related to a
        quantification.
        test data is tab file with no header.  Columns are:
        analysisId, multiCount, multiSpliceCount, totalReadCount, uniqueCount,
        uniqueSpliceCount
        """
        def __init__(self, readCountFileName):
            countData = open(readCountFileName, "r")
            counts = countData.readline()
            fields = counts.strip().split('\t')
            self.analysisId = fields[0]
            self.multiCount = int(fields[1])
            self.multiSpliceCount = int(fields[2])
            self.totalReadCount = int(fields[3])
            self.uniqueCount = int(fields[4])
            self.uniqueSpliceCount = int(fields[5])

    class ExpressionLevelInfo(object):
        """
        Container class for expression information
        test data is tab file with no header.  Columns are:
        id, annotationId, expression, featureGroupId,
        isNormalized, rawReadCount, score, units
        """
        def __init__(self, expressionFileName):
            expressionLevelData = open(expressionFileName, "r")
            self.expressionLevel = {}
            for expressionData in expressionLevelData.readlines():
                fields = expressionData.strip().split('\t')
                id = fields[0]
                self.expressionLevel[id] = {}
                self.expressionLevel[id]["annotationId"] = fields[1]
                self.expressionLevel[id]["expression"] = fields[2]
                self.expressionLevel[id]["featureGroupId"] = fields[3]
                self.expressionLevel[id]["isNormalized"] = fields[4]
                self.expressionLevel[id]["rawReadCount"] = fields[5]
                self.expressionLevel[id]["score"] = fields[6]
                self.expressionLevel[id]["units"] = fields[7]

    def __init__(self, rnaQuantificationId, baseDir):
        super(RnaQuantificationTest, self).__init__(rnaQuantificationId,
                                                    baseDir)
        self._rnaQuantFileName = os.path.join(baseDir, "rnaseq.table")
        self._rnaQuantInfo = self.RnaQuantInfo(self._rnaQuantFileName)
        self._characterizationFileName = os.path.join(baseDir, "dist.table")
        self._characterizationInfo = self.CharacterizationInfo(
            self._characterizationFileName)
        self._readCountFileName = os.path.join(baseDir, "counts.table")
        self._readCountInfo = self.ReadCountInfo(self._readCountFileName)
        self._expressionFileName = os.path.join(baseDir, "expression.table")
        self._expressionLevelInfo = self.ExpressionLevelInfo(
            self._expressionFileName)
        self._featureGroupInfo = self._getFeatureGroupInfo(
            self._expressionLevelInfo, rnaQuantificationId)

    def _getFeatureGroupInfo(self, expressionInfo, rnaQuantId):
        featureGroupInfo = {}
        for expressionId in expressionInfo.keys():
            id = expressionInfo[expressionId]["featureGroupId"]
            featureGroupInfo[id] = {"analysisId": rnaQuantId}
            featureGroupInfo[id] = {"name": id}

        return featureGroupInfo

    def getDataModelClass(self):
        return rna_quantification.RNASeqResult

    def getProtocolClass(self):
        return protocol.RnaQuantification

    def testValidateObjects(self):
        # test that validation works on reference sets and references
        rnaQuantification = self._gaObject
        rnaQuantificationPe = rnaQuantification.toProtocolElement()
        self.assertValid(
            protocol.RnaQuantification, rnaQuantificationPe.toJsonDict())

    def testGetRnaQuantification(self):
        # test searching with no arguments succeeds
        rnaQuantification = self._gaObject
        for gaRnaQuant in rnaQuantification.getRnaQuantification(None):
            self.assertRnaQuantsEqual(gaRnaQuant, self._rnaQuantInfo)

    def assertRnaQuantsEqual(self, gaRnaQuant, rnaQuant):
        self.assertEqual(gaRnaQuant.id, rnaQuant.id)
        self.assertEqual(gaRnaQuant.annotationIds, rnaQuant.annotationIds)
        self.assertEqual(gaRnaQuant.description, rnaQuant.description)
        self.assertEqual(gaRnaQuant.name, rnaQuant.name)
        self.assertEqual(gaRnaQuant.readGroupId, rnaQuant.readGroupId)

    def testGetCharacterization(self):
        rnaQuantification = self._gaObject
        for gaCharacterization in rnaQuantification.getCharacterization(None):
            self.assertCharacterizationEqual(gaCharacterization,
                                             self._characterizationInfo)

    def assertCharacterizationEqual(self, gaCharacterization,
                                    characterization):
        self.assertEqual(gaCharacterization.analysisId,
                         characterization.analysisId)
        self.assertEqual(gaCharacterization.complexity,
                         characterization.complexity)
        self.assertEqual(gaCharacterization.exonicFraction,
                         characterization.exonicFraction)
        self.assertEqual(gaCharacterization.fractionMapped,
                         characterization.fractionMapped)
        self.assertEqual(gaCharacterization.intergenicFraction,
                         characterization.intergenicFraction)
        self.assertEqual(gaCharacterization.intronicFraction,
                         characterization.intronicFraction)

    def testGetReadCounts(self):
        rnaQuantification = self._gaObject
        for gaReadCounts in rnaQuantification.getReadCounts(None):
            self.assertReadCountsEqual(gaReadCounts, self._readCountInfo)

    def assertReadCountsEqual(self, gaReadCounts, readCounts):
        self.assertEqual(gaReadCounts.analysisId, readCounts.analysisId)
        self.assertEqual(gaReadCounts.multiCount, readCounts.multiCount)
        self.assertEqual(gaReadCounts.multiSpliceCount,
                         readCounts.multiSpliceCount)
        self.assertEqual(gaReadCounts.totalReadCount,
                         readCounts.totalReadCount)
        self.assertEqual(gaReadCounts.uniqueCount, readCounts.uniqueCount)
        self.assertEqual(gaReadCounts.uniqueSpliceCount,
                         readCounts.uniqueSpliceCount)

    def testGetExpressionLevel(self):
        # test searching with no arguments succeeds
        rnaQuantification = self._gaObject
        for gaExpression in rnaQuantification.getExpressionLevel(None, None):
            self.assertExpressionEqual(gaExpression,
                                       self._expressionLevelInfo)

    def assertExpressionEqual(self, gaExpression, expressionLevel):
        id = gaExpression.id
        self.assertIn(id, expressionLevel.keys())
        self.assertEqual(gaExpression.annotationId,
                         expressionLevel[id]["annotationId"])
        self.assertEqual(gaExpression.expression,
                         expressionLevel[id]["expression"])
        self.assertEqual(gaExpression.featureGroupId,
                         expressionLevel[id]["featureGroupId"])
        self.assertEqual(gaExpression.isNormalized,
                         expressionLevel[id]["isNormalized"])
        self.assertEqual(gaExpression.rawReadCount,
                         expressionLevel[id]["rawReadCount"])
        self.assertEqual(gaExpression.score, expressionLevel[id]["score"])
        self.assertEqual(gaExpression.units, expressionLevel[id]["units"])

    def testGetFeatureGroup(self):
        # test searching with no arguments succeeds
        rnaQuantification = self._gaObject
        groupId = self.getRandomFeatureGroup()
        for gaFeatureGroup in rnaQuantification.getFeatureGroup(groupId):
            self.assertFeatureGroupEqual(gaFeatureGroup,
                                         self._featureGroupInfo)

    def getRandomFeatureGroup(self):
        expressionId = random.choice(self._expressionLevelInfo.keys())
        return self._expressionLevelInfo[expressionId]["featureGroupId"]

    def assertFeatureGroupEqual(self, gaFeatureGroup, featureGroupInfo):
        id = gaFeatureGroup.id
        self.assertIn(id, featureGroupInfo.keys())
        self.assertEqual(gaFeatureGroup.annotationId,
                         featureGroupInfo[id]["analysisId"])
        self.assertEqual(gaFeatureGroup.FeatureGroup,
                         featureGroupInfo[id]["name"])
