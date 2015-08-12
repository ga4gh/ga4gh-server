"""
Data-driven tests for references.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import hashlib

# TODO it may be a bit circular to use pysam as our interface for
# accessing reference information, since this is the method we use
# in the main code. We should examine the possibility of using
# a different library here --- perhaps BioPython?
import pysam

import ga4gh.datamodel.references as references
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
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
    class ReferenceInfo(object):
        """
        Container class for information about a reference
        """
        def __init__(self, referenceSetId, fastaFileName):
            self.fastaFile = pysam.FastaFile(fastaFileName)
            self.bases = self.fastaFile.fetch(self.fastaFile.references[0])
            self.length = len(self.bases)

    def __init__(self, referenceSetId, baseDir):
        super(ReferenceSetTest, self).__init__(None, referenceSetId, baseDir)
        self._referenceInfos = {}
        for fastaFileName in glob.glob(
                os.path.join(self._dataDir, "*.fa.gz")):
            self._readFasta(referenceSetId, fastaFileName)

    def _readFasta(self, referenceSetId, fastaFileName):
        referenceInfo = self.ReferenceInfo(referenceSetId, fastaFileName)
        self._referenceInfos[fastaFileName] = referenceInfo

    def getDataModelClass(self):
        class HtslibReferenceSetWrapper(references.HtslibReferenceSet):
            def __init__(self, parentContainer, localId, dataDir):
                super(HtslibReferenceSetWrapper, self).__init__(
                    localId, dataDir)
        return HtslibReferenceSetWrapper

    def getProtocolClass(self):
        return protocol.ReferenceSet

    def testValidateObjects(self):
        # test that validation works on reference sets and references
        referenceSet = self._gaObject
        referenceSetPe = referenceSet.toProtocolElement()
        self.assertValid(
            protocol.ReferenceSet, referenceSetPe.toJsonDict())
        self.assertGreater(len(referenceSetPe.referenceIds), 0)
        for gaReference in referenceSet.getReferences():
            reference = gaReference.toProtocolElement().toJsonDict()
            self.assertValid(protocol.Reference, reference)

    def testGetBases(self):
        # test searching with no arguments succeeds
        referenceSet = self._gaObject
        for gaReference in referenceSet.getReferences():
            pysamReference = self._referenceInfos[
                gaReference.getFastaFilePath()]
            self.assertReferencesEqual(gaReference, pysamReference)

    def testGetBasesStart(self):
        # test searching with start only succeeds
        self.doRangeTest(5, None)

    def testGetBasesEnd(self):
        # test searching with end only succeeds
        self.doRangeTest(None, 5)

    def testGetBasesRanges(self):
        # test searching with start and end succeeds
        self.doRangeTest(5, 10)

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
            # pysamReference = self._referenceInfos[
            #     gaReference.getFastaFilePath()]
            # basesChecksum = hashlib.md5(pysamReference.bases).hexdigest()
            self.assertEqual("TODO", gaReference.getMd5Checksum())
            referenceMd5s.append(gaReference.getMd5Checksum())
        checksumsString = ''.join(referenceMd5s)
        md5checksum = hashlib.md5(checksumsString).hexdigest()
        referenceSetMd5 = referenceSet.getMd5Checksum()
        self.assertEqual(md5checksum, referenceSetMd5)

    def doRangeTest(self, start=None, end=None):
        referenceSet = self._gaObject
        for gaReference in referenceSet.getReferences():
            pysamReference = self._referenceInfos[
                gaReference.getFastaFilePath()]
            self.assertReferencesEqual(
                gaReference, pysamReference, start, end)

    def assertReferencesEqual(
            self, gaReference, pysamReference, start=None, end=None):
        theStart = 0 if start is None else start
        theEnd = pysamReference.length if end is None else end
        gaBases = gaReference.getBases(theStart, theEnd)
        pysamBases = pysamReference.bases[theStart:theEnd]
        self.assertEqual(gaBases, pysamBases)
