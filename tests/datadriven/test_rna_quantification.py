"""
Data-driven tests for rna quantification.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import random
import sqlite3

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.protocol as protocol
import ga4gh.datamodel.rna_quantification as rna_quantification
import tests.datadriven as datadriven
import ga4gh.exceptions as exceptions


def testRnaQuantification():
    testDataDir = "tests/data/datasets/dataset1/rnaQuant"
    for test in datadriven.makeTests(testDataDir, RnaQuantificationTest,
                                     '*.db'):
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
        def __init__(self, dbConn):
            self._dbConn = dbConn
            sql = ("SELECT id FROM RNAQUANTIFICATION")
            query = self._dbconn.execute(sql)
            self.id = query.fetchone()
            sql = ("SELECT annotation_ids FROM RNAQUANTIFICATION")
            query = self._dbconn.execute(sql)
            self.annotationIds = query.fetchone().split(',')
            sql = ("SELECT description FROM RNAQUANTIFICATION")
            query = self._dbconn.execute(sql)
            self.description = query.fetchone()
            sql = ("SELECT name FROM RNAQUANTIFICATION")
            query = self._dbconn.execute(sql)
            self.name = query.fetchone()
            sql = ("SELECT read_group_id FROM RNAQUANTIFICATION")
            query = self._dbconn.execute(sql)
            self.readGroupId = query.fetchone()

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
        self._dbFile = os.path.join(baseDir, "ENCFF305LZB.db")
        self._dbconn = sqlite3.connect(self._dbFile)
        self._rnaQuantInfo = self.RnaQuantInfo(self._dbconn)

        # self._expressionLevelInfo = self.ExpressionLevelInfo(
        #     self._expressionFileName)
        # self._featureGroupInfo = self._getFeatureGroupInfo(
        #     self._expressionLevelInfo, self._gaObject.getId())

    def _getFeatureGroupInfo(self, expressionInfo, rnaQuantId):
        featureGroupInfo = {}
        for expressionId in expressionInfo.expressionLevel.keys():
            id = expressionInfo.expressionLevel[expressionId]["featureGroupId"]
            featureGroupInfo[id] = {"analysisId": rnaQuantId, "name": id}

        return featureGroupInfo

    def getDataModelInstance(self, localId, dataPath):
        return rna_quantification.RNASeqResult(self._dataset, localId,
                                               dataPath, None)

    def getProtocolClass(self):
        return protocol.RnaQuantification

    def testValidateObjects(self):
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
