"""
Data-driven tests for sequence annotation Features.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.server.datarepo as datarepo
import ga4gh.server.datamodel as datamodel
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.datamodel.references as references
import ga4gh.server.datamodel.sequence_annotations as sequence_annotations
import tests.datadriven as datadriven
import tests.paths as paths

import ga4gh.schemas.protocol as protocol

_datasetName = "ds"

_discontinuousTestData = {
    "featureSetName": "discontinuous",
    "referenceName": "apidb|Pf3D7_13",
    "totalFeatures": 30,
    "sampleFeatureId": 4409955920,
    "sampleParentId": 4409956112,
    "sampleStart": 820942,
    "sampleEnd": 821379,
    "sampleStrand": protocol.POS_STRAND,
    "sampleSiblings": 2,
    "region": [0, 2**32],
    "ontologyRestriction": ["gene", ],
    "featuresWithOntology": 3
}

_gencodeV21Set1TestData = {
    "featureSetName": "gencodeV21Set1",
    "referenceName": "chr1",
    "totalFeatures": 543,
    "sampleFeatureId": 4397111632,
    "sampleParentId": 4397111312,
    "sampleStart": 804776,
    "sampleEnd": 804832,
    "sampleStrand": protocol.POS_STRAND,
    "sampleSiblings": 5,
    "region": [0, 2**32],
    "ontologyRestriction": ["gene", ],
    "featuresWithOntology": 58
}

_sacCerTestTestData = {
    "featureSetName": "sacCerTest",
    "referenceName": "chrI",
    "totalFeatures": 33,
    "sampleFeatureId": 4324354512,
    "sampleParentId": None,
    "sampleStart": 337,
    "sampleEnd": 801,
    "sampleStrand": protocol.NEG_STRAND,
    "sampleSiblings": 33,
    "region": [0, 2**32],
    "ontologyRestriction": ["gene", ],
    "featuresWithOntology": 11
}

_specialCasesTestTestData = {
    "featureSetName": "specialCasesTest",
    "referenceName": "2L",
    "totalFeatures": 4,
    "sampleFeatureId": 4554844304,
    "sampleParentId": None,
    "sampleStart": 22229583,
    "sampleEnd": 22229699,
    "sampleStrand": protocol.NEG_STRAND,
    "sampleSiblings": 4,
    "region": [0, 2**32],
    "ontologyRestriction": ["gene", ],
    "featuresWithOntology": 0
}

_testDataForFeatureSetName = {
    "discontinuous": _discontinuousTestData,
    "gencodeV21Set1": _gencodeV21Set1TestData,
    "sacCerTest": _sacCerTestTestData,
    "specialCasesTest": _specialCasesTestTestData
}


def _getFeatureCompoundId(dataSetName, featureSetName, featureId):
    splits = [dataSetName, featureSetName, str(featureId)]
    joined = datamodel.CompoundId.join(splits)
    obfuscated = datamodel.CompoundId.obfuscate(joined)
    return obfuscated


def testFeatureSets():
    testDataDir = "tests/data/datasets/dataset1/sequenceAnnotations"
    for test in datadriven.makeTests(testDataDir, FeatureSetTests, "*.db"):
        yield test


class FeatureSetTests(datadriven.DataDrivenTest):
    """
    Re-parses source GFF3 files, compares the results with the contents
    of corresponding sequence_annotations.Feature objects.
    """
    def __init__(self, featureSetLocalName, dataPath):
        """
        :param localId: Name of the GFF3 resource corresponding to a pair
        of files, .db and .gff3
        :param dataPath: string representing full path to the .db file
        :return:
        """
        self._dataset = datasets.Dataset(_datasetName)
        self._repo = datarepo.SqlDataRepository(paths.testDataRepo)
        self._repo.open(datarepo.MODE_READ)
        self._ontology = self._repo.getOntologyByName(paths.ontologyName)
        self._referenceSet = references.AbstractReferenceSet("test_rs")
        featureSetLocalName = featureSetLocalName[:-3]  # remove '.db'
        self._testData = _testDataForFeatureSetName[featureSetLocalName]
        super(FeatureSetTests, self).__init__(featureSetLocalName, dataPath)

    def getProtocolClass(self):
        return protocol.FeatureSet

    def getDataModelInstance(self, localId, dataPath):
        featureSet = sequence_annotations.Gff3DbFeatureSet(
            self._dataset, localId)
        featureSet.setOntology(self._ontology)
        featureSet.setReferenceSet(self._referenceSet)
        featureSet.populateFromFile(dataPath)
        return featureSet

    def testGetFeatureById(self):
        idString = _getFeatureCompoundId(
            _datasetName,
            self._testData["featureSetName"],
            self._testData["sampleFeatureId"])
        compoundId = datamodel.FeatureCompoundId.parse(idString)
        feature = self._gaObject.getFeature(compoundId)
        self.assertIsNotNone(feature)
        self.assertEqual(feature.id, idString)
        if self._testData["sampleParentId"] is not None:
            self.assertEqual(
                datamodel.FeatureCompoundId.parse(
                    feature.parent_id).featureId,
                str(self._testData["sampleParentId"]))
        else:
            self.assertEqual(feature.parent_id, '')
        self.assertEqual(
            feature.reference_name,
            self._testData["referenceName"])
        self.assertEqual(
            feature.start,
            self._testData["sampleStart"])
        self.assertEqual(
            feature.end,
            self._testData["sampleEnd"])
        self.assertEqual(
            feature.strand,
            self._testData["sampleStrand"])

    def testFetchAllFeaturesInRegion(self):
        features = []
        for feature in self._gaObject.getFeatures(
                self._testData["referenceName"],
                self._testData["region"][0],
                self._testData["region"][1],
                None, 1000):
            features.append(feature)
        self.assertEqual(len(features), self._testData["totalFeatures"])

    def testFetchFeaturesRestrictedByOntology(self):
        features = []
        for feature in self._gaObject.getFeatures(
                self._testData["referenceName"],
                self._testData["region"][0],
                self._testData["region"][1],
                None, 1000,
                featureTypes=self._testData["ontologyRestriction"]):
            features.append(feature)
        self.assertEqual(len(features),
                         self._testData["featuresWithOntology"])

    def testFetchFeaturesRestrictedByParent(self):
        parentId = ""
        if self._testData["sampleParentId"] is not None:
            parentIdString = _getFeatureCompoundId(
                _datasetName,
                self._testData["featureSetName"],
                self._testData["sampleParentId"])
            parentId = datamodel.FeatureCompoundId.parse(
                parentIdString).featureId
        features = []
        for feature in self._gaObject.getFeatures(
                self._testData["referenceName"],
                self._testData["region"][0],
                self._testData["region"][1],
                None, 1000,
                parentId=parentId):
            features.append(feature)
        self.assertEqual(len(features),
                         self._testData["sampleSiblings"])
