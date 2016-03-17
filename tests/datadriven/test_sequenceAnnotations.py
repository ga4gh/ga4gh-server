"""
Data-driven tests for sequence annotation Features.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.backend as backend
import ga4gh.datarepo as datarepo
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.sequenceAnnotations as sequenceAnnotations
import ga4gh.protocol as protocol
import tests.datadriven as datadriven

_discontinuousTestData = {
    "featureSetId": "YnJjYTE6ZGlzY29udGludW91cw",
    "referenceName": "apidb|Pf3D7_13",
    "totalFeatures": 21,
    "sampleFeatureId": ("YnJjYTE6ZGlzY29udGludW91czphcGlk"
                        "YnxleG9uX01BTDEzUDEuMTAzLTEw"),
    "sampleParentId": ("YnJjYTE6ZGlzY29udGludW91czphcGlk"
                       "YnxybmFfTUFMMTNQMS4xMDMtMQ"),
    "sampleStart": 805985,
    "sampleEnd": 806230,
    "sampleStrand": protocol.Strand.POS_STRAND,
    "sampleSiblings": 15,
    "regionRestriction": [0, 0],
    "featuresInRegion": 0,
    "ontologyRestriction": "gene",
    "featuresWithOntology": 0
}

_gencodeV21Set1TestData = {
    "featureSetId": "YnJjYTE6Z2VuY29kZVYyMVNldDE",
    "referenceName": "chr1",
    "totalFeatures": 542,
    "sampleFeatureId": ("YnJjYTE6Z2VuY29kZVYyMVNldDE6ZXhvb"
                        "jtFTlNUMDAwMDA1OTA4NDguMzs1"),
    "sampleParentId": ("YnJjYTE6Z2VuY29kZVYyMVNldDE6RU5TVD"
                       "AwMDAwNTkwODQ4LjM"),
    "sampleStart": 804776,
    "sampleEnd": 804832,
    "sampleStrand": protocol.Strand.POS_STRAND,
    "sampleSiblings": 5,
    "regionRestriction": [0, 0],
    "featuresInRegion": 0,
    "ontologyRestriction": "gene",
    "featuresWithOntology": 0

}

_sacCerTestTestData = {
    "featureSetId": "YnJjYTE6c2FjQ2VyVGVzdA",
    "referenceName": "chrI",
    "totalFeatures": 21,
    "sampleFeatureId": ("YnJjYTE6c2FjQ2VyVGVzdDpURUwwMUwtWEM"),
    "sampleParentId": (""),
    "sampleStart": 337,
    "sampleEnd": 801,
    "sampleStrand": protocol.Strand.NEG_STRAND,
    "sampleSiblings": 0,
    "regionRestriction": [0, 0],
    "featuresInRegion": 0,
    "ontologyRestriction": "gene",
    "featuresWithOntology": 0

}

_specialCasesTestTestData = {
    "featureSetId": "",
    "referenceName": "",
    "totalFeatures": 0,
    "sampleFeatureId": (""),
    "sampleParentId": (""),
    "sampleStart": 0,
    "sampleEnd": 0,
    "sampleStrand": protocol.Strand.POS_STRAND,
    "sampleSiblings": 0,
    "regionRestriction": [0, 0],
    "featuresInRegion": 0,
    "ontologyRestriction": "gene",
    "featuresWithOntology": 0

}

_testDataForFeatureSetId = {
    "discontinuous.db": _discontinuousTestData,
    "gencodeV21Set1.db": _gencodeV21Set1TestData,
    "sacCerTest.db": _sacCerTestTestData,
    "specialCasesTest.db": _specialCasesTestTestData
}


def testFeatureSets():
    testDataDir = "tests/data/datasets/dataset1/sequenceAnnotations"
    for test in datadriven.makeTests(testDataDir, FeatureSetTests):
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
        self._backend = backend.Backend(datarepo.AbstractDataRepository())
        self._dataset = datasets.AbstractDataset("ds")
        self._datarepo = datarepo.AbstractDataRepository()
        self._testData = _testDataForFeatureSetId[featureSetLocalId]
        super(FeatureSetTests, self).__init__(featureSetLocalId, dataPath)

    def getProtocolClass(self):
        return protocol.FeatureSet

    def getDataModelInstance(self, localId, dataPath):
        return sequenceAnnotations.Gff3DbFeatureSet(
            self._dataset, localId, dataPath, self._datarepo)

    def testGetFeatureById(self):
        pass

    def testFetchAllFeatures(self):
        pass

    def testFetchFeaturesInRegion(self):
        pass

    def testFetchFeaturesRestrictedByOntology(self):
        pass

    def testFetchFeaturesRestrictedByParent(self):
        pass
