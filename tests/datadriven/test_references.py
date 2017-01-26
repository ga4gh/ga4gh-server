"""
Data-driven tests for references.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import hashlib
import os
import unittest

# TODO it may be a bit circular to use pysam as our interface for
# accessing reference information, since this is the method we use
# in the main code. We should examine the possibility of using
# a different library here --- perhaps BioPython?
import pysam

import ga4gh.server.datamodel.references as references
import ga4gh.server.exceptions as exceptions
import tests.datadriven as datadriven
import tests.paths as paths

import ga4gh.schemas.protocol as protocol


def testReferenceSets():
    testDataDir = os.path.join(paths.testDataDir, "referenceSets")
    pattern = "*.fa.gz"
    for test in datadriven.makeTests(testDataDir, ReferenceSetTest, pattern):
        yield test


class ReferenceSetTest(datadriven.DataDrivenTest):
    """
    Data drive test class for reference sets. Builds an alternative model of
    a reference set, and verifies that it is consistent with the model
    built by the references.ReferenceSet object.
    """
    def __init__(self, referenceSetId, fastaFile):
        super(ReferenceSetTest, self).__init__(referenceSetId, fastaFile)
        self._fastaFile = pysam.FastaFile(fastaFile)

    def getDataModelInstance(self, localId, dataPath):
        referenceSet = references.HtslibReferenceSet(localId)
        referenceSet.populateFromFile(dataPath)
        return referenceSet

    def getProtocolClass(self):
        return protocol.ReferenceSet

    def testValidateObjects(self):
        # test that validation works on reference sets and references
        referenceSet = self._gaObject
        referenceSetPe = referenceSet.toProtocolElement()
        self.assertValid(
            protocol.ReferenceSet, protocol.toJson(referenceSetPe))
        for gaReference in referenceSet.getReferences():
            reference = protocol.toJson(gaReference.toProtocolElement())
            self.assertValid(protocol.Reference, reference)

    def testGetBases(self):
        # test searching with no arguments succeeds
        referenceSet = self._gaObject
        for gaReference in referenceSet.getReferences():
            self.assertReferencesEqual(gaReference)

    def testGetBasesStart(self):
        # test searching with start only succeeds
        self.doRangeTest(5, None)

    def testGetBasesEnd(self):
        # test searching with end only succeeds
        self.doRangeTest(None, 5)

    def testGetBasesRanges(self):
        # test searching with start and end succeeds
        self.doRangeTest(2, 5)

    @unittest.skip("We assume that the 0 is unset, as protobuf3 has 0 == None")
    def testGetBasesEmpty(self):
        self.doRangeTest(0, 0)

    def testOutOfBounds(self):
        referenceSet = self._gaObject
        for reference in referenceSet.getReferences():
            length = reference.getLength()
            badRanges = [
                (-1, 1), (0, length + 1), (0, 2**34), (-2**32, 1),
                (1, 0), (length, length - 1),
            ]
            for start, end in badRanges:
                self.assertRaises(
                    exceptions.ReferenceRangeErrorException,
                    reference.getBases, start, end)

    def testMd5checksums(self):
        referenceSet = self._gaObject
        referenceMd5s = []
        for gaReference in referenceSet.getReferences():
            bases = self._fastaFile.fetch(gaReference.getLocalId())
            basesChecksum = hashlib.md5(bases).hexdigest()
            self.assertEqual(basesChecksum, gaReference.getMd5Checksum())
            referenceMd5s.append(gaReference.getMd5Checksum())
        referenceMd5s.sort()
        checksumsString = ''.join(referenceMd5s)
        md5checksum = hashlib.md5(checksumsString).hexdigest()
        referenceSetMd5 = referenceSet.getMd5Checksum()
        self.assertEqual(md5checksum, referenceSetMd5)

    def doRangeTest(self, start=None, end=None):
        referenceSet = self._gaObject
        for gaReference in referenceSet.getReferences():
            self.assertReferencesEqual(gaReference, start, end)

    def assertReferencesEqual(
            self, gaReference, start=None, end=None):
        theStart = 0 if start is None else start
        bases = self._fastaFile.fetch(gaReference.getLocalId())
        theEnd = len(bases) if end is None else end
        gaBases = gaReference.getBases(theStart, theEnd)
        pysamBases = bases[theStart:theEnd]
        self.assertEqual(gaBases, pysamBases)
