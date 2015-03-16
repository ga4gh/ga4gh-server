"""
Data-driven tests for reads
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.protocol as protocol
import ga4gh.datamodel.reads as reads
import tests.datadriven as datadriven


def testVariantSets():
    testDataDir = "tests/data/reads"
    for test in datadriven.makeTests(testDataDir, ReadGroupSetTest):
        yield test


class ReadGroupSetTest(datadriven.DataDrivenTest):

    def getDataModelClass(self):
        return reads.ReadGroupSet

    def getProtocolClass(self):
        return protocol.GAReadGroupSet

    def testFixMe(self):
        readGroupSet = self._gaObject
        # TODO
        self.assertIsNotNone(readGroupSet)
