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


def testFeatureSets():
    testDataDir = "tests/data/datasets/dataset1/sequenceAnnotations"
    for test in datadriven.makeTests(testDataDir, FeatureSetDataValidityTests):
        yield test


class FeatureSetDataValidityTests(datadriven.DataDrivenTest):
    """
    Re-parses source GFF3 files, compares the results with the contents
    of corresponding sequenceAnnotations.Feature objects.
    """
    def __init__(self, localId, dataPath):
        """
        :param localId: Name of the GFF3 resource corresponding to a pair
        of files, .db and .gff3
        :param dataPath: string representing full path to the .db file
        :return:
        """
        self._dataset = datasets.AbstractDataset("ds")
        self._backend = backend.Backend(datarepo.AbstractDataRepository())
        self._datarepo = datarepo.FileSystemDataRepository("tests/data")
        self._referenceSet = None
        self._gff3Records = []
        super(FeatureSetDataValidityTests, self).__init__(localId, dataPath)

    def getProtocolClass(self):
        return protocol.FeatureSet

    def getDataModelInstance(self, localId, dataPath):
        return sequenceAnnotations.Gff3DbFeatureSet(
            self._dataset, localId, dataPath, self._datarepo)
