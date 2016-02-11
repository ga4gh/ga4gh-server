"""
Data-driven tests for rna quantification.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import random

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.protocol as protocol
import ga4gh.datamodel.rna_quantification as rna_quantification
import tests.datadriven as datadriven
import ga4gh.exceptions as exceptions


def testRnaQuantification():
    testDataDir = "tests/data/datasets/dataset1/rnaQuant"
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
                self.expressionLevel[id] = {
                    "annotationId": fields[1],
                    "expression": fields[2],
                    "featureGroupId": fields[3],
                    "isNormalized": fields[4],
                    "rawReadCount": fields[5],
                    "score": fields[6],
                    "units": fields[7]
                }

    def __init__(self, rnaQuantificationId, baseDir):
        self._dataset = datasets.AbstractDataset("ds")
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
            self._expressionLevelInfo, self._gaObject.getId())

    def _getFeatureGroupInfo(self, expressionInfo, rnaQuantId):
        featureGroupInfo = {}
        for expressionId in expressionInfo.expressionLevel.keys():
            id = expressionInfo.expressionLevel[expressionId]["featureGroupId"]
            featureGroupInfo[id] = {"analysisId": rnaQuantId, "name": id}

        return featureGroupInfo

    def getDataModelInstance(self, localId, dataPath):
        return rna_quantification.RNASeqResult(self._dataset, localId,
                                               dataPath)

    def getProtocolClass(self):
        return protocol.RnaQuantification

    def testValidateObjects(self):
        # test that validation works on reference sets and references
        rnaQuantification = self._gaObject
        rnaQuantificationPe = rnaQuantification.toProtocolElement()
        self.assertValid(
            protocol.RnaQuantification, rnaQuantificationPe.toJsonDict())

    def assertRnaQuantsEqual(self, gaRnaQuant, rnaQuant):
        self.assertEqual(gaRnaQuant.id, rnaQuant.id)
        self.assertEqual(gaRnaQuant.annotationIds, rnaQuant.annotationIds)
        self.assertEqual(gaRnaQuant.description, rnaQuant.description)
        self.assertEqual(gaRnaQuant.name, rnaQuant.name)
        self.assertEqual(gaRnaQuant.readGroupId, rnaQuant.readGroupId)

    def testGetCharacterization(self):
        rnaQuantification = self._gaObject
        gaChar = rnaQuantification.getCharacterization()
        self.assertCharacterizationEqual(gaChar.toProtocolElement(),
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
        gaReadCounts = rnaQuantification.getReadCounts()
        self.assertReadCountsEqual(gaReadCounts.toProtocolElement(),
                                   self._readCountInfo)

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
        rnaQuantification = self._gaObject
        # positive test: get the expected expression level
        expressionID = self._expressionLevelInfo.expressionLevel.keys()[0]
        compoundId = datamodel.ExpressionLevelCompoundId(
            rnaQuantification.getCompoundId(), expressionID)
        gaExpression = rnaQuantification.getExpressionLevel(str(compoundId),
                                                            None)
        self.assertExpressionEqual(gaExpression,
                                   self._expressionLevelInfo.expressionLevel)

        # negative test: bad name
        with self.assertRaises(exceptions.ExpressionLevelNotFoundException):
            rnaQuantification.getExpressionLevel("not here", None)

    def assertExpressionEqual(self, gaExpressionObj,
                              expressionLevel):
        id = gaExpressionObj.getName()
        gaExpression = gaExpressionObj.toProtocolElement()
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
        rnaQuantification = self._gaObject
        # positive test: get the expected feature group
        groupId = self.getRandomFeatureGroupId()
        compoundId = datamodel.FeatureGroupCompoundId(
            rnaQuantification.getCompoundId(), groupId)
        gaFeatureGroup = rnaQuantification.getFeatureGroup(str(compoundId))
        self.assertFeatureGroupEqual(gaFeatureGroup, self._featureGroupInfo)

        # negative test: bad name
        with self.assertRaises(exceptions.FeatureGroupNotFoundException):
            rnaQuantification.getFeatureGroup("not here")

    def getRandomFeatureGroupId(self):
        idList = self._expressionLevelInfo.expressionLevel.keys()
        expressionId = random.choice(idList)
        expression = self._expressionLevelInfo.expressionLevel[expressionId]
        return expression["featureGroupId"]

    def assertFeatureGroupEqual(self, gaFeatureGroupObj,
                                featureGroupInfo):
        gaFeatureGroup = gaFeatureGroupObj.toProtocolElement()
        id = gaFeatureGroupObj.getLocalId()
        self.assertIn(id, featureGroupInfo.keys())
        self.assertEqual(gaFeatureGroup.analysisId,
                         featureGroupInfo[id]["analysisId"])
        self.assertEqual(gaFeatureGroup.name,
                         featureGroupInfo[id]["name"])
