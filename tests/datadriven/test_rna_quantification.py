"""
Data-driven tests for rna quantification.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.protocol as protocol
import ga4gh.datamodel.rna_quantification as rna_quantification
import tests.datadriven as datadriven
import ga4gh.exceptions as exceptions


_datasetName = "ds"


_rnaQuantTestData = {
    "annotation_ids": ["Gencodev16"],
    "description": "RNAseq data from ENCODE evaluation",
    "name": "ENCFF305LZB",
    "read_group_id": ""
}


_expressionTestData = {
    "bad_id": "MWtnLXAzLXN1YnNldDpybmFfZXhhbXBsZV8yOm1tOV9leGFtcGxlXzI=",
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
    "analysisId": "ENCFF305LZB",
    "name": "ENSG00000076984.13"
}


def _getRnaQuantCompoundId(dataSetName, rnaQuant):
    return datamodel.CompoundId.obfuscate(":".join(
        [dataSetName, rnaQuant]))


def _getExpressionCompoundId(dataSetName, rnaQuant, expressionId):
    return datamodel.CompoundId.obfuscate(":".join(
        [dataSetName, rnaQuant, expressionId]))


def testRnaQuantification():
    testDataDir = "tests/data/datasets/dataset1/rnaQuant"
    for test in datadriven.makeTests(
            testDataDir, RnaQuantificationTest, '*.db'):
        yield test


class RnaQuantificationTest(datadriven.DataDrivenTest):
    """
    Data driven test class for rna quantification. Builds an alternative model
    of a rna quantification, and verifies that it is consistent with the model
    built by the rna_quantification.RNASeqResult object.
    """
    def __init__(self, rnaQuantificationLocalId, baseDir):
        self._dataset = datasets.AbstractDataset("ds")
        rnaQuantificationId = rnaQuantificationLocalId[:-3]  # remove '.db'
        super(RnaQuantificationTest, self).__init__(
            rnaQuantificationId, baseDir)

    def getDataModelInstance(self, localId, dataPath):
        return rna_quantification.RNASeqResult(
            self._dataset, localId, dataPath, None)

    def getProtocolClass(self):
        return protocol.RnaQuantification

    def testValidateObjects(self):
        rnaQuantification = self._gaObject
        rnaQuantificationPe = rnaQuantification.toProtocolElement()
        self.assertValid(
            protocol.RnaQuantification, rnaQuantificationPe.toJsonDict())

    def testRnaQuantificationObject(self):
        gaRnaQuant = self._gaObject.toProtocolElement()
        idString = _getRnaQuantCompoundId(
            _datasetName,
            _rnaQuantTestData["name"])
        compoundId = datamodel.RnaQuantificationCompoundId.parse(idString)
        self.assertEqual(gaRnaQuant.id, str(compoundId))
        self.assertEqual(
            gaRnaQuant.annotationIds, _rnaQuantTestData["annotation_ids"])
        self.assertEqual(
            gaRnaQuant.description, _rnaQuantTestData["description"])
        self.assertEqual(gaRnaQuant.name, _rnaQuantTestData["name"])
        self.assertEqual(
            gaRnaQuant.readGroupId, _rnaQuantTestData["read_group_id"])

    def testGetExpressionLevelById(self):
        rnaQuantification = self._gaObject
        idString = _getExpressionCompoundId(
            _datasetName,
            _rnaQuantTestData["name"],
            _expressionTestData["name"])
        compoundId = datamodel.ExpressionLevelCompoundId.parse(idString)
        gaExpression = rnaQuantification.getExpressionLevel(compoundId)
        self.assertExpressionEqual(gaExpression, _expressionTestData)

    def testGetExpressionLevelByBadIdFails(self):
        rnaQuantification = self._gaObject
        badId = datamodel.ExpressionLevelCompoundId(
            rnaQuantification.getCompoundId(), "bad_id")
        with self.assertRaises(exceptions.ExpressionLevelNotFoundException):
            rnaQuantification.getExpressionLevel(badId)

    def assertExpressionEqual(self, gaExpressionObj, testData):
        gaExpression = gaExpressionObj.toProtocolElement()
        idString = _getExpressionCompoundId(
            _datasetName,
            _rnaQuantTestData["name"],
            _expressionTestData["name"])
        compoundId = datamodel.ExpressionLevelCompoundId.parse(idString)
        self.assertEqual(gaExpression.id, str(compoundId))
        self.assertEqual(
            gaExpression.annotationId, testData["annotation_id"])
        self.assertEqual(
            gaExpression.expression, testData["expression"])
        self.assertEqual(
            gaExpression.featureGroupId, testData["feature_group_id"])
        self.assertEqual(
            gaExpression.isNormalized, testData["is_normalized"])
        self.assertEqual(
            gaExpression.rawReadCount, testData["raw_read_count"])
        self.assertEqual(gaExpression.score, testData["score"])
        self.assertEqual(gaExpression.units, testData["units"])

    def testGetFeatureGroupById(self):
        rnaQuantification = self._gaObject
        idString = _getExpressionCompoundId(
            _datasetName,
            _rnaQuantTestData["name"],
            _featureGroupTestData["name"])
        compoundId = datamodel.FeatureGroupCompoundId.parse(idString)
        gaFeatureGroup = rnaQuantification.getFeatureGroup(
            compoundId)
        self.assertFeatureGroupEqual(
            gaFeatureGroup, _featureGroupTestData)

    def testGetFeatureGroupByBadIdFails(self):
        rnaQuantification = self._gaObject
        badId = datamodel.FeatureGroupCompoundId.parse(
            _expressionTestData["bad_id"])
        with self.assertRaises(
                exceptions.FeatureGroupNotFoundException):
            rnaQuantification.getFeatureGroup(badId)

    def assertFeatureGroupEqual(
            self, gaFeatureGroupObj, testData):
        gaFeatureGroup = gaFeatureGroupObj.toProtocolElement()
        idString = _getExpressionCompoundId(
            _datasetName,
            _rnaQuantTestData["name"],
            _featureGroupTestData["name"])
        compoundId = datamodel.FeatureGroupCompoundId.parse(idString)
        self.assertEqual(gaFeatureGroup.id, str(compoundId))
        self.assertEqual(
            gaFeatureGroup.analysisId, testData["analysisId"])
        self.assertEqual(gaFeatureGroup.name, testData["name"])
