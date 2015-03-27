"""
Data-driven tests for reads
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import glob
import os

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
import ga4gh.datamodel.reads as reads
import tests.datadriven as datadriven

import pysam


def testReads():
    testDataDir = "tests/data/reads"
    for test in datadriven.makeTests(testDataDir, ReadGroupSetTest):
        yield test


class ReadGroupSetTest(datadriven.DataDrivenTest):
    """
    Data driven test for read group sets
    """
    class ReadGroupInfo(object):
        """
        Container class for information about a read group
        """
        def __init__(self, readGroupSetId, samFileName):
            filename = os.path.split(samFileName)[1]
            localId = os.path.splitext(filename)[0]
            readGroupId = "{}:{}".format(readGroupSetId, localId)
            self.id = readGroupId
            self.samFile = pysam.AlignmentFile(samFileName)
            self.reads = []
            self.refNames = collections.defaultdict(list)
            self.refIds = collections.defaultdict(list)
            for read in self.samFile:
                refName = self.samFile.getrname(read.reference_id)
                self.refNames[refName].append(read)
                refId = read.reference_id
                self.refIds[refId].append(read)
                self.reads.append(read)

    def __init__(self, readGroupSetId, baseDir):
        super(ReadGroupSetTest, self).__init__(readGroupSetId, baseDir)
        self._readGroupInfos = {}
        for samFileName in glob.glob(
                os.path.join(self._dataDir, "*.bam")):
            self._readSam(readGroupSetId, samFileName)

    def _readSam(self, readGroupSetId, samFileName):
        readGroupInfo = self.ReadGroupInfo(readGroupSetId, samFileName)
        self._readGroupInfos[samFileName] = readGroupInfo

    def getDataModelClass(self):
        return reads.HtslibReadGroupSet

    def getProtocolClass(self):
        return protocol.GAReadGroupSet

    def testGetReadAlignmentsBothRefs(self):
        # test that querying by both referenceName and referenceId fails
        with self.assertRaises(exceptions.BadReadsSearchRequestBothRefs):
            readGroupSet = self._gaObject
            for readGroup in readGroupSet.getReadGroups():
                list(readGroup.getReadAlignments("a", 5, 0))

    def testGetReadAlignments(self):
        # test that searching with no arguments succeeds
        for gaAlignment, pysamAlignment, readGroup, readGroupInfo in \
                self.getReadAlignmentsGenerator():
            self.assertAlignmentsEqual(
                gaAlignment, pysamAlignment, readGroupInfo)

    def testValidateObjects(self):
        # test that validation works on read groups and reads
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            protocol.GAReadGroup.validate(
                readGroup.toProtocolElement().toJsonDict())
            alignments = list(readGroup.getReadAlignments())
            for gaAlignment in alignments:
                protocol.GAReadAlignment.validate(gaAlignment.toJsonDict())

    def testGetReadAlignmentsRefName(self):
        # test that searching with a reference name succeeds
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getSamFilePath()]
            for refName, refNameReads in readGroupInfo.refNames.items():
                alignments = list(readGroup.getReadAlignments(refName))
                for gaAlignment, pysamAlignment in zip(
                        alignments, refNameReads):
                    self.assertAlignmentsEqual(
                        gaAlignment, pysamAlignment, readGroupInfo)

    def testGetReadAlignmentsRefId(self):
        # test that searching with a reference id succeeds
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getSamFilePath()]
            for refId, refIdReads in readGroupInfo.refIds.items():
                alignments = list(
                    readGroup.getReadAlignments(referenceId=refId))
                for gaAlignment, pysamAlignment in zip(
                        alignments, refIdReads):
                    self.assertAlignmentsEqual(
                        gaAlignment, pysamAlignment, readGroupInfo)

    def testGetReadAlignmentsStartEnd(self):
        # test that searching with start and end coords succeeds
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getSamFilePath()]
            for refName, refNameReads in readGroupInfo.refNames.items():
                bigNumThatPysamWontChokeOn = 2**30
                alignments = list(readGroup.getReadAlignments(
                    refName, None, 0, bigNumThatPysamWontChokeOn))
                for gaAlignment, pysamAlignment in zip(
                        alignments, refNameReads):
                    self.assertAlignmentsEqual(
                        gaAlignment, pysamAlignment, readGroupInfo)

    def testGetReadAlignmentSearchRanges(self):
        # test that various range searches work
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getSamFilePath()]
            for refName, refNameReads in readGroupInfo.refNames.items():
                alignments = list(readGroup.getReadAlignments(refName))
                length = len(alignments)
                if length < 2:
                    continue
                positions = [read.alignment.position.position
                             for read in alignments]
                if length != len(set(positions)):
                    continue
                begin = positions[0]
                beginLength = len(alignments[0].alignedSequence)
                end = positions[-1]
                self.assertGetReadAlignmentsRangeResult(
                    readGroup, refName, begin, end + 1, length)
                self.assertGetReadAlignmentsRangeResult(
                    readGroup, refName, begin, end, length - 1)
                self.assertGetReadAlignmentsRangeResult(
                    readGroup, refName, begin, begin, 0)
                self.assertGetReadAlignmentsRangeResult(
                    readGroup, refName, begin + beginLength,
                    end + 1, length - 1)
                self.assertGetReadAlignmentsRangeResult(
                    readGroup, refName, begin + beginLength,
                    end, length - 2)

    def assertGetReadAlignmentsRangeResult(self, readGroup, refName,
                                           start, end, result):
        alignments = list(readGroup.getReadAlignments(
            refName, None, start, end))
        self.assertEqual(len(alignments), result)

    def assertAlignmentsEqual(self, gaAlignment, pysamAlignment,
                              readGroupInfo):
        self.assertEqual(
            gaAlignment.alignedQuality,
            list(pysamAlignment.query_qualities))
        self.assertEqual(
            gaAlignment.alignedSequence,
            pysamAlignment.query_sequence)
        self.assertEqual(
            gaAlignment.alignment.mappingQuality,
            pysamAlignment.mapping_quality)
        self.assertEqual(
            gaAlignment.alignment.position.referenceName,
            readGroupInfo.samFile.getrname(pysamAlignment.reference_id))
        self.assertEqual(
            gaAlignment.alignment.position.position,
            pysamAlignment.reference_start)
        self.assertCigarEqual(
            gaAlignment.alignment.cigar,
            pysamAlignment.cigar)
        self.assertFlag(
            gaAlignment.duplicateFragment,
            pysamAlignment, reads.SamFlags.DUPLICATE_FRAGMENT)
        self.assertFlag(
            gaAlignment.failedVendorQualityChecks,
            pysamAlignment, reads.SamFlags.FAILED_VENDOR_QUALITY_CHECKS)
        self.assertEqual(
            gaAlignment.fragmentLength,
            pysamAlignment.template_length)
        self.assertEqual(
            gaAlignment.fragmentName,
            pysamAlignment.query_name)
        self.assertEqual(
            gaAlignment.id,
            "{}:{}".format(readGroupInfo.id, pysamAlignment.query_name))
        self.assertEqual(
            gaAlignment.info,
            dict(pysamAlignment.tags))
        self.assertEqual(
            gaAlignment.nextMatePosition.position,
            pysamAlignment.next_reference_start)
        if pysamAlignment.next_reference_id != -1:
            self.assertEqual(
                gaAlignment.nextMatePosition.referenceName,
                readGroupInfo.samFile.getrname(
                    pysamAlignment.next_reference_id))
        else:
            self.assertIsNone(gaAlignment.nextMatePosition.referenceName)
        self.assertFlag(
            gaAlignment.numberReads,
            pysamAlignment, reads.SamFlags.NUMBER_READS)
        self.assertFlag(
            gaAlignment.properPlacement,
            pysamAlignment, reads.SamFlags.PROPER_PLACEMENT)
        self.assertEqual(
            gaAlignment.readGroupId,
            readGroupInfo.id)
        self.assertFlag(
            gaAlignment.readNumber,
            pysamAlignment, reads.SamFlags.READ_NUMBER_ONE)
        self.assertFlag(
            gaAlignment.readNumber,
            pysamAlignment, reads.SamFlags.READ_NUMBER_TWO)
        self.assertFlag(
            gaAlignment.secondaryAlignment,
            pysamAlignment, reads.SamFlags.SECONDARY_ALIGNMENT)
        self.assertFlag(
            gaAlignment.supplementaryAlignment,
            pysamAlignment, reads.SamFlags.SUPPLEMENTARY_ALIGNMENT)

    def assertFlag(self, gaAlignmentAttr, pysamAlignment, mask):
        flagSet = reads.SamFlags.isFlagSet(pysamAlignment.flag, mask)
        self.assertEqual(gaAlignmentAttr, flagSet)

    def assertCigarEqual(self, gaCigar, pysamCigar):
        self.assertEqual(len(gaCigar), len(pysamCigar))
        for i, gaCigarUnit in enumerate(gaCigar):
            operation, length = pysamCigar[i]
            gaCigarUnitOperation = reads.SamCigar.ga2int(
                gaCigarUnit.operation)
            self.assertEqual(
                gaCigarUnitOperation, operation)
            self.assertEqual(
                gaCigarUnit.operationLength, length)

    def getReadAlignmentsGenerator(self):
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getSamFilePath()]
            alignments = list(readGroup.getReadAlignments())
            for gaAlignment, pysamAlignment in zip(
                    alignments, readGroupInfo.reads):
                yield gaAlignment, pysamAlignment, readGroup, readGroupInfo
