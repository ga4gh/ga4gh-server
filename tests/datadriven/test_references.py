"""
Data-driven tests for references.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import glob
import hashlib
import json
import os

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
    testDataDir = "tests/data/referenceSets"
    for test in datadriven.makeTests(testDataDir, ReferenceSetTest):
        yield test


class ReferenceSetInfo(object):
    """
    Container class for information about a reference set
    """
    def __init__(self, dataDir):
        metadataFilePath = '{}.json'.format(dataDir)
        with open(metadataFilePath) as metadataFile:
            self.metadata = json.load(metadataFile)


class ReferenceInfo(object):
    """
    Container class for information about a reference
    """
    def __init__(self, referenceName, dataDir):
        fastaFileName = os.path.join(
            dataDir, "{}.fa.gz".format(referenceName))
        self.fastaFile = pysam.FastaFile(fastaFileName)
        assert len(self.fastaFile.references) == 1
        assert self.fastaFile.references[0] == referenceName
        self.bases = self.fastaFile.fetch(self.fastaFile.references[0])
        self.length = len(self.bases)
        metadataFileName = os.path.join(
            dataDir, "{}.json".format(referenceName))
        with open(metadataFileName) as metadataFile:
            self.metadata = json.load(metadataFile)


class ReferenceSetTest(datadriven.DataDrivenTest):
    """
    Data drive test class for reference sets. Builds an alternative model of
    a reference set, and verifies that it is consistent with the model
    built by the references.ReferenceSet object.
    """
    def __init__(self, referenceSetId, baseDir):
        super(ReferenceSetTest, self).__init__(referenceSetId, baseDir)
        self._referenceSetInfo = ReferenceSetInfo(baseDir)
        self._referenceInfos = {}
        for fastaFilePath in glob.glob(
                os.path.join(self._dataPath, "*.fa.gz")):
            fastaFileName = os.path.split(fastaFilePath)[1]
            referenceName = fastaFileName.split(".")[0]
            referenceInfo = ReferenceInfo(referenceName, self._dataPath)
            self._referenceInfos[referenceName] = referenceInfo

    def getDataModelInstance(self, localId, dataPath):
        return references.HtslibReferenceSet(localId, dataPath, None)

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

    def testMetadata(self):
        # Check that the metadata loaded from the JSON file is
        # consistent with returned objects.
        referenceSet = self._gaObject
        metadata = self._referenceSetInfo.metadata
        self.assertEqual(
            referenceSet.getAssemblyId(), metadata['assemblyId'])
        self.assertEqual(
            referenceSet.getDescription(), metadata['description'])
        self.assertEqual(
            referenceSet.getIsDerived(), metadata['isDerived'])
        self.assertEqual(
            referenceSet.getNcbiTaxonId(), metadata['ncbiTaxonId'])
        self.assertEqual(
            referenceSet.getSourceAccessions(), metadata['sourceAccessions'])
        self.assertEqual(
            referenceSet.getSourceUri(), metadata['sourceUri'])
        for reference in referenceSet.getReferences():
            referenceInfo = self._referenceInfos[reference.getLocalId()]
            metadata = referenceInfo.metadata
            self.assertEqual(
                reference.getNcbiTaxonId(), metadata["ncbiTaxonId"])
            self.assertEqual(
                reference.getMd5Checksum(), metadata["md5checksum"])
            self.assertEqual(
                reference.getSourceUri(), metadata["sourceUri"])
            self.assertEqual(
                reference.getIsDerived(), metadata["isDerived"])
            self.assertEqual(
                reference.getSourceDivergence(), metadata["sourceDivergence"])
            self.assertEqual(
                reference.getSourceAccessions(), metadata["sourceAccessions"])

    def testGetBases(self):
        # test searching with no arguments succeeds
        referenceSet = self._gaObject
        for gaReference in referenceSet.getReferences():
            pysamReference = self._referenceInfos[gaReference.getLocalId()]
            self.assertReferencesEqual(gaReference, pysamReference)

    def testGetBasesStart(self):
        # test searching with start only succeeds
        self.doRangeTest(5, None)

    def testGetBasesEnd(self):
        # test searching with end only succeeds
        self.doRangeTest(None, 5)

    def testGetBasesRanges(self):
        # test searching with start and end succeeds
        self.doRangeTest(2, 5)

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
            pysamReference = self._referenceInfos[gaReference.getLocalId()]
            basesChecksum = hashlib.md5(pysamReference.bases).hexdigest()
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
            pysamReference = self._referenceInfos[gaReference.getLocalId()]
            self.assertReferencesEqual(
                gaReference, pysamReference, start, end)

    def assertReferencesEqual(
            self, gaReference, pysamReference, start=None, end=None):
        theStart = 0 if start is None else start
        theEnd = pysamReference.length if end is None else end
        gaBases = gaReference.getBases(theStart, theEnd)
        pysamBases = pysamReference.bases[theStart:theEnd]
        self.assertEqual(gaBases, pysamBases)
