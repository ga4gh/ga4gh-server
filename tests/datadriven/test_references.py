"""
Data-driven tests for references.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob

# TODO it may be a bit circular to use pysam as our interface for
# accessing reference information, since this is the method we use
# in the main code. We should examine the possibility of using
# a different library here --- perhaps BioPython?
import pysam

import ga4gh.protocol as protocol
import ga4gh.datamodel.references as references
import tests.datadriven as datadriven


def testReferenceSets():
    testDataDir = "tests/data/references"
    for test in datadriven.makeTests(testDataDir, ReferenceSetTest):
        yield test


class ReferenceSetTest(datadriven.DataDrivenTest):
    """
    Data drive test class for reference sets. Builds an alternative model of
    a reference set, and verifies that it is consistent with the model
    built by the references.ReferenceSet object.
    """
    def __init__(self, referenceSetId, baseDir):
        super(ReferenceSetTest, self).__init__(referenceSetId, baseDir)
        # Read in all the FASTA files in dataDir. Each reference within
        # the reference set maps to a single FASTA file.
        self._idFastaFileMap = {}
        for path in glob.glob(os.path.join(self._dataDir, "*.fa.gz")):
            filename = os.path.split(path)[1]
            referenceId = "{}:{}".format(
                referenceSetId, filename.split(".")[0])
            self._idFastaFileMap[referenceId] = pysam.FastaFile(path)

    def getDataModelClass(self):
        return references.LinearReferenceSet

    def getProtocolClass(self):
        return protocol.ReferenceSet

    def testNumReferences(self):
        references = list(self._gaObject.getReferences())
        self.assertEqual(len(self._idFastaFileMap), len(references))
