"""
Data-driven tests for g2p.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.datamodel.genotype_phenotype as g2p
import ga4gh.datamodel.datasets as datasets
import ga4gh.protocol as protocol
import tests.datadriven as datadriven
import tests.paths as paths


def testG2P():
    testDataDir = os.path.join(
        paths.testDataDir, "datasets/dataset1/phenotypes")
    for test in datadriven.makeTests(testDataDir, PhenotypeAssociationSetTest):
        yield test


class PhenotypeAssociationSetTest(datadriven.DataDrivenTest):
    def __init__(self, localId, baseDir):
        self._dataset = datasets.AbstractDataset("ds")
        super(PhenotypeAssociationSetTest, self).__init__(localId, baseDir)
        # TODO compare
        pass

    def getDataModelInstance(self, localId, dataPath):
        return g2p.PhenotypeAssociationSet(self._dataset, localId, dataPath)

    def getProtocolClass(self):
        return protocol.PhenotypeAssociationSet

    def testQuery(self):
        # TODO the output of G2Pdataset should be compared
        # with the results of some RDF query/munging here
        pass
