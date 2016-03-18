"""
Data-driven tests for sequence annotation Features.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.backend as backend
import ga4gh.datarepo as datarepo
import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.sequenceAnnotations as sequenceAnnotations
import ga4gh.protocol as protocol
import tests.datadriven as datadriven
import ga4gh.exceptions as exceptions

_discontinuousTestData = {
    "featureSetId": "ZHM6ZGlzY29udGludW91cw",
    "referenceName": "apidb|Pf3D7_13",
    "totalFeatures": 21,
    "sampleFeatureId": ("ZHM6ZGlzY29udGludW91czphcGlk"
                        "YnxleG9uX01BTDEzUDEuMTAzLTEw"),
    "sampleParentId": ("ZHM6ZGlzY29udGludW91czphcGlk"
                       "YnxybmFfTUFMMTNQMS4xMDMtMQ"),
    "sampleStart": 805985,
    "sampleEnd": 806230,
    "sampleStrand": protocol.Strand.POS_STRAND,
    "sampleSiblings": 11,
    "region": [0, 2**32],
    "ontologyRestriction": ["gene", ],
    "featuresWithOntology": 3
}

_gencodeV21Set1TestData = {
    "featureSetId": "ZHM6Z2VuY29kZVYyMVNldDE",
    "referenceName": "chr1",
    "totalFeatures": 542,
    "sampleFeatureId": ("ZHM6Z2VuY29kZVYyMVNldDE6ZXhvb"
                        "jtFTlNUMDAwMDA1OTA4NDguMzs1"),
    "sampleParentId": ("ZHM6Z2VuY29kZVYyMVNldDE6RU5TVD"
                       "AwMDAwNTkwODQ4LjM"),
    "sampleStart": 804776,
    "sampleEnd": 804832,
    "sampleStrand": protocol.Strand.POS_STRAND,
    "sampleSiblings": 5,
    "region": [0, 2**32],
    "ontologyRestriction": ["gene", ],
    "featuresWithOntology": 58
}

_sacCerTestTestData = {
    "featureSetId": "ZHM6c2FjQ2VyVGVzdA",
    "referenceName": "chrI",
    "totalFeatures": 21,
    "sampleFeatureId": ("ZHM6c2FjQ2VyVGVzdDpURUwwMUwtWEM"),
    "sampleParentId": (""),
    "sampleStart": 337,
    "sampleEnd": 801,
    "sampleStrand": protocol.Strand.NEG_STRAND,
    "sampleSiblings": 20,
    "region": [0, 2**32],
    "ontologyRestriction": ["gene", ],
    "featuresWithOntology": 11
}

_specialCasesTestTestData = {
    "featureSetId": "ZHM6c3BlY2lhbENhc2VzVGVzdA",
    "referenceName": "2L",
    "totalFeatures": 3,
    "sampleFeatureId": ("ZHM6c3BlY2lhbENhc2VzVGVzdDo7MTU3Nzg3N"
                        "19ibGFzdHhfbWFza2VkX2FhX1NQVFIuZG1lbA"),
    "sampleParentId": (""),
    "sampleStart": 22229583,
    "sampleEnd": 22229699,
    "sampleStrand": protocol.Strand.NEG_STRAND,
    "sampleSiblings": 2,
    "region": [0, 2**32],
    "ontologyRestriction": ["gene", ],
    "featuresWithOntology": 0
}

_testDataForFeatureSetId = {
    "discontinuous": _discontinuousTestData,
    "gencodeV21Set1": _gencodeV21Set1TestData,
    "sacCerTest": _sacCerTestTestData,
    "specialCasesTest": _specialCasesTestTestData
}


def testFeatureSets():
    testDataDir = "tests/data/datasets/dataset1/sequenceAnnotations"
    for test in datadriven.makeTests(testDataDir, FeatureSetTests, "*.db"):
        yield test


class FeatureSetTests(datadriven.DataDrivenTest):
    """
    Re-parses source GFF3 files, compares the results with the contents
    of corresponding sequenceAnnotations.Feature objects.
    """
    def __init__(self, featureSetLocalId, dataPath):
        """
        :param localId: Name of the GFF3 resource corresponding to a pair
        of files, .db and .gff3
        :param dataPath: string representing full path to the .db file
        :return:
        """
        self._dataset = datasets.AbstractDataset("ds")
        self._datarepo = datarepo.FileSystemDataRepository("tests/data")
        self._backend = backend.Backend(self._datarepo)
        featureSetLocalId = featureSetLocalId[:-3]  # remove '.db'
        self._testData = _testDataForFeatureSetId[featureSetLocalId]
        super(FeatureSetTests, self).__init__(featureSetLocalId, dataPath)

    def getProtocolClass(self):
        return protocol.FeatureSet

    def getDataModelInstance(self, localId, dataPath):
        return sequenceAnnotations.Gff3DbFeatureSet(
            self._dataset, localId, dataPath, self._datarepo)

    def testGetFeatureById(self):
        idString = self._testData["sampleFeatureId"]
        compoundId = datamodel.FeatureCompoundId.parse(idString)
        feature = self._gaObject.getFeature(compoundId)
        self.assertIsNotNone(feature)
        self.assertEqual(feature.id, idString)
        self.assertEqual(feature.parentId,
                         self._testData["sampleParentId"])
        self.assertEqual(feature.referenceName,
                         self._testData["referenceName"])
        self.assertEqual(feature.start,
                         self._testData["sampleStart"])
        self.assertEqual(feature.end,
                         self._testData["sampleEnd"])
        self.assertEqual(feature.strand,
                         self._testData["sampleStrand"])

    def testGetFeatureFailsWithBadId(self):
        idString = self._testData["sampleFeatureId"] + "W00t"
        try:
            compoundId = datamodel.FeatureCompoundId.parse(idString)
            self._gaObject.getFeature(compoundId)
            self.assertFalse("Exception should be thrown by this point")
        except exceptions.ObjectWithIdNotFoundException:
            pass

    def testFetchAllFeaturesInRegion(self):
        features = []
        nextPageTokens = []
        for (feature, nextPageToken) in self._gaObject.getFeatures(
                self._testData["referenceName"],
                self._testData["region"][0],
                self._testData["region"][1],
                None, 1000):
            features.append(feature)
            nextPageTokens.append(nextPageToken)
        self.assertEqual(len(features), self._testData["totalFeatures"])
        self.assertIsNone(nextPageTokens[-1])

    def testFetchFeaturesRestrictedByOntology(self):
        features = []
        for (feature, _) in self._gaObject.getFeatures(
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
        if self._testData["sampleParentId"] != "":
            parentId = datamodel.FeatureCompoundId.parse(
                self._testData["sampleParentId"]).featureId
        features = []
        for (feature, _) in self._gaObject.getFeatures(
                self._testData["referenceName"],
                self._testData["region"][0],
                self._testData["region"][1],
                None, 1000,
                parentId=parentId):
            features.append(feature)
        self.assertEqual(len(features),
                         self._testData["sampleSiblings"])
