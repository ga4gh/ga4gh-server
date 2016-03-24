"""
Data-driven tests for rna quantification.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# import os
import random
# import sqlite3

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.protocol as protocol
import ga4gh.datamodel.rna_quantification as rna_quantification
import tests.datadriven as datadriven
import ga4gh.exceptions as exceptions


_rnaQuantTestData = {
    "id": "ZHM6RU5DRkYzMDVMWkI",
    "annotation_ids": ["Gencodev16"],
    "description": "RNAseq data from ENCODE evaluation",
    "name": "ENCFF305LZB",
    "read_group_id": ""
}


_expressionTestData = {
    "id": "ZHM6RU5DRkYzMDVMWkI6RU5TRzAwMDAwMDc2OTg0LjEz",
    "feature_compound_id": "ZHM6RU5DRkYzMDVMWkI6RU5TRzAwMDAwMDc2OTg0LjEz",
    "name": "ENSG00000076984.13",
    "annotation_id": "Gencodev16",
    "expression": 24.52,
    "feature_group_id": "ENSG00000076984.13",
    "is_normalized": True,
    "raw_read_count": 4317.0,
    "score": 23.34315,
    "units": "TPM"

}


_featureGroupTestData = {
    "id": "ZHM6RU5DRkYzMDVMWkI6RU5TRzAwMDAwMDc2OTg0LjEz",
    "analysisId": "ENCFF305LZB",
    "name": "ENSG00000076984.13"
}


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
    def __init__(self, rnaQuantificationLocalId, baseDir):
        self._dataset = datasets.AbstractDataset("ds")
        # self._datarepo = datarepo.FileSystemDataRepository("tests/data")
        # self._backend = backend.Backend(self._datarepo)
        rnaQuantificationId = rnaQuantificationLocalId[:-3]  # remove '.db'
        super(RnaQuantificationTest, self).__init__(rnaQuantificationId,
                                                    baseDir)
        # self._dbFile = os.path.join(baseDir, "ENCFF305LZB.db")
        # self._dbconn = sqlite3.connect(self._dbFile)
        # self._rnaQuantInfo = self.RnaQuantInfo(self._dbconn)

        # self._expressionLevelInfo = self.ExpressionLevelInfo(
        #     self._expressionFileName)
        # self._featureGroupInfo = self._getFeatureGroupInfo(
        #     self._expressionLevelInfo, self._gaObject.getId())

    def getDataModelInstance(self, localId, dataPath):
        return rna_quantification.RNASeqResult(self._dataset, localId,
                                               dataPath, None)

    def getProtocolClass(self):
        return protocol.RnaQuantification

    def _getFeatureGroupInfo(self, expressionInfo, rnaQuantId):
        featureGroupInfo = {}
        for expressionId in expressionInfo.expressionLevel.keys():
            id = expressionInfo.expressionLevel[expressionId]["featureGroupId"]
            featureGroupInfo[id] = {"analysisId": rnaQuantId, "name": id}

        return featureGroupInfo

    def testValidateObjects(self):
        rnaQuantification = self._gaObject
        rnaQuantificationPe = rnaQuantification.toProtocolElement()
        self.assertValid(
            protocol.RnaQuantification, rnaQuantificationPe.toJsonDict())

    def testRnaQuantificationObject(self):
        gaRnaQuant = self._gaObject.toProtocolElement()
        self.assertEqual(gaRnaQuant.id, _rnaQuantTestData["id"])
        self.assertEqual(gaRnaQuant.annotationIds,
                         _rnaQuantTestData["annotation_ids"])
        self.assertEqual(gaRnaQuant.description,
                         _rnaQuantTestData["description"])
        self.assertEqual(gaRnaQuant.name, _rnaQuantTestData["name"])
        self.assertEqual(gaRnaQuant.readGroupId,
                         _rnaQuantTestData["read_group_id"])

    def testGetExpressionLevelById(self):
        rnaQuantification = self._gaObject
        expressionId = _expressionTestData["id"]
        gaExpression = rnaQuantification.getExpressionLevel(expressionId)
        self.assertExpressionEqual(gaExpression, _expressionTestData)

    def testGetExpressionLevelByBadIdFails(self):
        rnaQuantification = self._gaObject
        rawId = datamodel.ExpressionLevelCompoundId.deobfuscate(
            _expressionTestData["id"])
        badId = rawId + "not_in_database"
        expressionId = datamodel.ExpressionLevelCompoundId.obfuscate(badId)
        with self.assertRaises(exceptions.ExpressionLevelNotFoundException):
            rnaQuantification.getExpressionLevel(expressionId)

    def assertExpressionEqual(self, gaExpressionObj, testData):
        gaExpression = gaExpressionObj.toProtocolElement()
        self.assertEqual(gaExpression.id, testData["id"])
        self.assertEqual(gaExpression.annotationId,
                         testData["annotation_id"])
        self.assertEqual(gaExpression.expression,
                         testData["expression"])
        self.assertEqual(gaExpression.featureGroupId,
                         testData["feature_group_id"])
        self.assertEqual(gaExpression.isNormalized,
                         testData["is_normalized"])
        self.assertEqual(gaExpression.rawReadCount,
                         testData["raw_read_count"])
        self.assertEqual(gaExpression.score, testData["score"])
        self.assertEqual(gaExpression.units, testData["units"])

    def testGetFeatureGroupById(self):
        # TODO: fix me
        rnaQuantification = self._gaObject
        groupId = _expressionTestData["feature_compound_id"]
        gaFeatureGroup = rnaQuantification.getFeatureGroup(groupId)
        self.assertFeatureGroupEqual(gaFeatureGroup, _featureGroupTestData)

    def testGetFeatureGroupByBadIdFails(self):
        rnaQuantification = self._gaObject
        rawId = datamodel.FeatureGroupCompoundId.deobfuscate(
            _expressionTestData["feature_compound_id"])
        badId = rawId + "not_in_database"
        featureGroupId = datamodel.FeatureGroupCompoundId.obfuscate(badId)
        with self.assertRaises(exceptions.FeatureGroupNotFoundException):
            rnaQuantification.getFeatureGroup(featureGroupId)

    def getRandomFeatureGroupId(self):
        idList = self._expressionLevelInfo.expressionLevel.keys()
        expressionId = random.choice(idList)
        expression = self._expressionLevelInfo.expressionLevel[expressionId]
        return expression["featureGroupId"]

    def assertFeatureGroupEqual(self, gaFeatureGroupObj, testData):
        gaFeatureGroup = gaFeatureGroupObj.toProtocolElement()
        self.assertEqual(gaFeatureGroup.id, testData["id"])
        self.assertEqual(gaFeatureGroup.analysisId, testData["analysisId"])
        self.assertEqual(gaFeatureGroup.name, testData["name"])
