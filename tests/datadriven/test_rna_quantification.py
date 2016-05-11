"""
Data-driven tests for rna quantification.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import ga4gh.datarepo as datarepo
import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.rna_quantification as rna_quantification
import ga4gh.protocol as protocol
import tests.datadriven as datadriven
import tests.paths as paths


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
    "quantification_group_id": "ENSG00000076984.13",
    "is_normalized": True,
    "raw_read_count": 4317.0,
    "score": 23.34315,
    "units": "TPM",
    "conf_low": 24.1,
    "conf_hi": 24.6,
    "num_expression_entries": 2,
    "num_entries_over_threshold": 1
}


_quantificationGroupTestData = {
    "analysisId": "ENCFF305LZB",
    "name": "ENSG00000076984.13",
    "num_quantification_group_entries": 2
}


def _getRnaQuantCompoundId(dataSetName, rnaQuant):
    splits = [dataSetName, rnaQuant]
    joined = datamodel.CompoundId.join(splits)
    obfuscated = datamodel.CompoundId.obfuscate(joined)
    return obfuscated


def _getExpressionCompoundId(dataSetName, rnaQuant, expressionId):
    splits = [dataSetName, rnaQuant, expressionId]
    joined = datamodel.CompoundId.join(splits)
    obfuscated = datamodel.CompoundId.obfuscate(joined)
    return obfuscated


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
        self._dataset = datasets.Dataset(_datasetName)
        self._repo = datarepo.SqlDataRepository(paths.testDataRepo)
        self._repo.open(datarepo.MODE_READ)
        rnaQuantificationId = rnaQuantificationLocalId[:-3]  # remove '.db'
        super(RnaQuantificationTest, self).__init__(
            rnaQuantificationId, baseDir)

    def getDataModelInstance(self, localId, dataPath):
        return rna_quantification.RNASeqResult(
            self._dataset, localId, dataPath)

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
            gaExpression.quantificationGroupId,
            testData["quantification_group_id"])
        self.assertEqual(
            gaExpression.isNormalized, testData["is_normalized"])
        self.assertEqual(
            gaExpression.rawReadCount, testData["raw_read_count"])
        self.assertEqual(gaExpression.score, testData["score"])
        self.assertEqual(gaExpression.units, testData["units"])
        self.assertEqual(
            gaExpression.confInterval,
            [testData["conf_low"], testData["conf_hi"]])

    def testSearchExpressionLevels(self):
        rnaQuantification = self._gaObject
        rnaQuantID = rnaQuantification.getLocalId()
        expressionLevels = rnaQuantification.getExpressionLevels(rnaQuantID)
        self.assertEqual(
            _expressionTestData["num_expression_entries"],
            len(expressionLevels))
        overThreshold = rnaQuantification.getExpressionLevels(
            rnaQuantID, threshold=100.0)
        self.assertEqual(
            _expressionTestData["num_entries_over_threshold"],
            len(overThreshold))

    def testGetQuantificationGroupById(self):
        rnaQuantification = self._gaObject
        idString = _getExpressionCompoundId(
            _datasetName,
            _rnaQuantTestData["name"],
            _quantificationGroupTestData["name"])
        compoundId = datamodel.QuantificationGroupCompoundId.parse(idString)
        gaQuantificationGroup = rnaQuantification.getQuantificationGroup(
            compoundId)
        self.assertQuantificationGroupEqual(
            gaQuantificationGroup, _quantificationGroupTestData)

    def assertQuantificationGroupEqual(
            self, gaQuantificationGroupObj, testData):
        gaQuantificationGroup = gaQuantificationGroupObj.toProtocolElement()
        idString = _getExpressionCompoundId(
            _datasetName,
            _rnaQuantTestData["name"],
            _quantificationGroupTestData["name"])
        compoundId = datamodel.QuantificationGroupCompoundId.parse(idString)
        self.assertEqual(gaQuantificationGroup.id, str(compoundId))
        self.assertEqual(
            gaQuantificationGroup.analysisId, testData["analysisId"])
        self.assertEqual(gaQuantificationGroup.name, testData["name"])

    def testSearchQuantificationGroups(self):
        rnaQuantification = self._gaObject
        rnaQuantID = rnaQuantification.getLocalId()
        quantificationGroups = rnaQuantification.getQuantificationGroups(
            rnaQuantID)
        self.assertEqual(
            _quantificationGroupTestData["num_quantification_group_entries"],
            len(quantificationGroups))
