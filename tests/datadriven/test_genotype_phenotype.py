"""
Data-driven tests for g2p.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.datamodel.genotype_phenotype as g2p
import ga4gh.protocol as protocol
import tests.datadriven as datadriven
import tests.paths as paths


def testG2P():
    testDataDir = os.path.join(
        paths.testDataDir, "g2pdatasets")
    for test in datadriven.makeTests(testDataDir, G2PDatasetTest):
        yield test


class G2PDatasetTest(datadriven.DataDrivenTest):
    def __init__(self, localId, baseDir):
        super(G2PDatasetTest, self).__init__(localId, baseDir)
        # TODO compare
        pass

    def getDataModelInstance(self, localId, dataPath):
        print(localId, dataPath)
        return g2p.G2PDataset(setName=localId, relativePath=dataPath)

    def getProtocolClass(self):
        # there is no G2P set protocol type so we return an
        # empty annotation
        return protocol.FeaturePhenotypeAssociation

    def testQuery(self):
        # TODO the output of G2Pdataset should be compared
        # with the results of some RDF query/munging here
        pass