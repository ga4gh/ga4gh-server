"""
Data-driven tests for reads
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections

import ga4gh.backend as backend
import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.references as references
import ga4gh.protocol as protocol
import tests.datadriven as datadriven
import tests.utils as utils

import pysam


def testReads():
    testDataDir = "tests/data/datasets/dataset1/reads"
    for test in datadriven.makeTests(
            testDataDir, ReadGroupSetTest, '*.bam'):
        yield test


class ReadGroupSetInfo(object):
    """
    Container class for information about a read group set
    """
    def __init__(self, samFile):
        self.numAlignedReads = samFile.mapped
        self.numUnalignedReads = samFile.unmapped


class ReadGroupInfo(object):
    """
    Container class for information about a read group
    """
    def __init__(self, gaReadGroupSet, samFile, readGroupName):
        self.gaReadGroup = reads.AbstractReadGroup(
            gaReadGroupSet, readGroupName)
        self.id = self.gaReadGroup.getId()
        self.samFile = samFile
        self.mappedReads = collections.defaultdict(list)
        for read in self.samFile:
            tags = dict(read.tags)
            if 'RG' not in tags or tags['RG'] != readGroupName:
                continue
            if read.reference_id != -1:
                # mapped read
                referenceName = self.samFile.getrname(read.reference_id)
                self.mappedReads[referenceName].append(read)
        self.numAlignedReads = -1
        self.numUnalignedReads = -1
        self.programs = []
        if 'PG' in self.samFile.header:
            self.programs = self.samFile.header['PG']
        self.sampleId = None
        self.description = None
        self.predictedInsertSize = None
        self.instrumentModel = None
        self.sequencingCenter = None
        self.experimentDescription = None
        self.library = None
        self.platformUnit = None
        self.runTime = None
        if 'RG' in self.samFile.header:
            readGroupHeader = [
                rgHeader for rgHeader in self.samFile.header['RG']
                if rgHeader['ID'] == readGroupName][0]
            self.sampleId = readGroupHeader.get('SM', None)
            self.description = readGroupHeader.get('DS', None)
            if 'PI' in readGroupHeader:
                self.predictedInsertSize = int(readGroupHeader['PI'])
            self.instrumentModel = readGroupHeader.get('PL', None)
            self.sequencingCenter = readGroupHeader.get('CN', None)
            self.experimentDescription = readGroupHeader.get('DS', None)
            self.library = readGroupHeader.get('LB', None)
            self.platformUnit = readGroupHeader.get('PU', None)
            self.runTime = readGroupHeader.get('DT', None)


class ReadGroupSetTest(datadriven.DataDrivenTest):
    """
    Data driven test for read group sets
    """
    def __init__(self, localId, dataPath):
        self._backend = backend.AbstractBackend()
        self._referenceSet = None
        self._dataset = datasets.AbstractDataset("ds")
        self._readGroupInfos = {}
        self._readGroupSetInfo = None
        self._samFile = pysam.AlignmentFile(dataPath)
        self._readReferences()
        super(ReadGroupSetTest, self).__init__(localId, dataPath)
        self._readAlignmentInfo()

    def _readReferences(self):
        # Read the reference information from the samfile
        referenceSetName = None
        for referenceInfo in self._samFile.header['SQ']:
            if 'AS' not in referenceInfo:
                infoDict = reads.parseMalformedBamHeader(referenceInfo)
            # If there's still no reference set name in there we use
            # a default name.
            name = infoDict.get("AS", "Default")
            if referenceSetName is None:
                referenceSetName = name
                self._addReferenceSet(referenceSetName)
            else:
                self.assertEqual(referenceSetName, name)
            self._addReference(infoDict['SN'])

    def _addReferenceSet(self, referenceSetName):
        self._referenceSet = references.AbstractReferenceSet(referenceSetName)
        self._backend.addReferenceSet(self._referenceSet)

    def _addReference(self, referenceName):
        reference = references.AbstractReference(
            self._referenceSet, referenceName)
        self._referenceSet.addReference(reference)

    def _readAlignmentInfo(self):
        self._readGroupSetInfo = ReadGroupSetInfo(self._samFile)
        if 'RG' in self._samFile.header:
            readGroupHeaders = self._samFile.header['RG']
            readGroupNames = [
                readGroupHeader['ID'] for readGroupHeader
                in readGroupHeaders]
        else:
            readGroupNames = ['default']
        for readGroupName in readGroupNames:
            readGroupInfo = ReadGroupInfo(
                self._gaObject, self._samFile, readGroupName)
            self._readGroupInfos[readGroupName] = readGroupInfo

    def getDataModelInstance(self, localId, dataPath):
        return reads.HtslibReadGroupSet(
            self._dataset, localId, dataPath, self._backend)

    def getProtocolClass(self):
        return protocol.ReadGroupSet

    def testSampleIdEtc(self):
        # test that sampleId and other misc fields are set correctly
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getLocalId()]
            gaReadGroup = readGroup.toProtocolElement()
            self.assertEqual(
                readGroupInfo.sampleId,
                gaReadGroup.sampleId)
            self.assertEqual(
                readGroupInfo.predictedInsertSize,
                gaReadGroup.predictedInsertSize)
            self.assertEqual(
                readGroupInfo.description,
                gaReadGroup.description)

    def testExperiments(self):
        # test that the experiment field is set correctly
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getLocalId()]
            gaReadGroup = readGroup.toProtocolElement()
            self.assertIn(
                "experiment",
                datamodel.CompoundId.deobfuscate(gaReadGroup.experiment.id))
            self.assertEqual(
                readGroupInfo.instrumentModel,
                gaReadGroup.experiment.instrumentModel)
            self.assertEqual(
                readGroupInfo.sequencingCenter,
                gaReadGroup.experiment.sequencingCenter)
            self.assertEqual(
                readGroupInfo.experimentDescription,
                gaReadGroup.experiment.description)
            self.assertEqual(
                readGroupInfo.library,
                gaReadGroup.experiment.library)
            self.assertEqual(
                readGroupInfo.platformUnit,
                gaReadGroup.experiment.platformUnit)
            self.assertEqual(
                readGroupInfo.runTime,
                gaReadGroup.experiment.runTime)

    def testPrograms(self):
        # test that program info is set correctly
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getLocalId()]
            gaPrograms = readGroup.getPrograms()
            htslibPrograms = readGroupInfo.programs
            for gaProgram, htslibProgram in utils.zipLists(
                    gaPrograms, htslibPrograms):
                self.assertEqual(
                    gaProgram.id, htslibProgram.get('ID'))
                self.assertEqual(
                    gaProgram.commandLine, htslibProgram.get('CL', None))
                self.assertEqual(
                    gaProgram.name, htslibProgram.get('PN', None))
                self.assertEqual(
                    gaProgram.prevProgramId, htslibProgram.get('PP', None))
                self.assertEqual(
                    gaProgram.version, htslibProgram.get('VN', None))

    def testReadGroupStats(self):
        # test that the stats attrs are populated correctly
        readGroupSet = self._gaObject
        gaReadGroupSet = readGroupSet.toProtocolElement()
        readGroupSetInfo = self._readGroupSetInfo
        self.assertEqual(
            readGroupSet.getNumAlignedReads(),
            readGroupSetInfo.numAlignedReads)
        self.assertEqual(
            readGroupSet.getNumUnalignedReads(),
            readGroupSetInfo.numUnalignedReads)
        self.assertEqual(
            gaReadGroupSet.stats.alignedReadCount,
            readGroupSetInfo.numAlignedReads)
        self.assertEqual(
            gaReadGroupSet.stats.unalignedReadCount,
            readGroupSetInfo.numUnalignedReads)
        for readGroup in readGroupSet.getReadGroups():
            gaReadGroup = readGroup.toProtocolElement()
            self.assertEqual(
                readGroup.getNumAlignedReads(), -1)
            self.assertEqual(
                readGroup.getNumUnalignedReads(), -1)
            self.assertEqual(
                gaReadGroup.stats.alignedReadCount, -1)
            self.assertEqual(
                gaReadGroup.stats.unalignedReadCount, -1)

    def testValidateObjects(self):
        # test that validation works on read groups and reads
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            self.assertValid(
                protocol.ReadGroup,
                readGroup.toProtocolElement().toJsonDict())
            for reference in self._referenceSet.getReferences():
                for gaAlignment in readGroup.getReadAlignments(reference):
                    self.assertValid(
                        protocol.ReadAlignment,
                        gaAlignment.toJsonDict())

    def testGetReadAlignmentsRefId(self):
        # test that searching with a reference id succeeds
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getLocalId()]
            for name, alignments in readGroupInfo.mappedReads.items():
                reference = self._referenceSet.getReferenceByName(name)
                self.assertAlignmentListsEqual(
                    list(readGroup.getReadAlignments(reference)), alignments,
                    readGroupInfo)

    def testGetReadAlignmentsStartEnd(self):
        # test that searching with start and end coords succeeds
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getLocalId()]
            for name, alignments, in readGroupInfo.mappedReads.items():
                bigNumThatPysamWontChokeOn = 2**30
                reference = self._referenceSet.getReferenceByName(name)
                gaAlignments = list(readGroup.getReadAlignments(
                    reference, 0, bigNumThatPysamWontChokeOn))
                self.assertAlignmentListsEqual(
                    gaAlignments, alignments, readGroupInfo)

    def testGetReadAlignmentSearchRanges(self):
        # test that various range searches work
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getLocalId()]
            for name in readGroupInfo.mappedReads.keys():
                reference = self._referenceSet.getReferenceByName(name)
                alignments = list(readGroup.getReadAlignments(reference))
                length = len(alignments)
                if length < 2:
                    continue
                positions = [read.alignment.position.position
                             for read in alignments]
                if length != len(set(positions)):
                    continue
                begin = positions[0]
                end = positions[-1]
                self.assertGetReadAlignmentsRangeResult(
                    readGroup, reference, begin, end + 1, length)
                self.assertGetReadAlignmentsRangeResult(
                    readGroup, reference, begin, end, length - 1)
                self.assertGetReadAlignmentsRangeResult(
                    readGroup, reference, begin, begin, 0)

    def assertGetReadAlignmentsRangeResult(
            self, readGroup, reference, start, end, result):
        alignments = list(readGroup.getReadAlignments(reference, start, end))
        self.assertEqual(len(alignments), result)

    def assertAlignmentListsEqual(
            self, gaAlignments, pysamAlignments, readGroupInfo):
        for gaAlignment, pysamAlignment in utils.zipLists(
                gaAlignments, pysamAlignments):
            self.assertAlignmentsEqual(
                gaAlignment, pysamAlignment, readGroupInfo)

    def assertAlignmentsEqual(self, gaAlignment, pysamAlignment,
                              readGroupInfo):
        if pysamAlignment.query_qualities is None:
            self.assertEqual(gaAlignment.alignedQuality, [])
        else:
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
        # TODO test reverseStrand on position and on nextMatePosition once
        # it has been implemented.
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
        compoundId = datamodel.ReadAlignmentCompoundId(
            readGroupInfo.gaReadGroup.getCompoundId(),
            pysamAlignment.query_name)
        self.assertEqual(gaAlignment.id, str(compoundId))
        self.assertEqual(
            gaAlignment.info,
            {key: [str(value)] for key, value in pysamAlignment.tags})
        if pysamAlignment.next_reference_id != -1:
            self.assertEqual(
                gaAlignment.nextMatePosition.position,
                pysamAlignment.next_reference_start)
            self.assertEqual(
                gaAlignment.nextMatePosition.referenceName,
                readGroupInfo.samFile.getrname(
                    pysamAlignment.next_reference_id))
        else:
            self.assertIsNone(gaAlignment.nextMatePosition)
        self.assertFlag(
            gaAlignment.properPlacement,
            pysamAlignment, reads.SamFlags.PROPER_PLACEMENT)
        self.assertEqual(
            gaAlignment.readGroupId,
            readGroupInfo.id)
        self.assertFlag(
            gaAlignment.secondaryAlignment,
            pysamAlignment, reads.SamFlags.SECONDARY_ALIGNMENT)
        self.assertFlag(
            gaAlignment.supplementaryAlignment,
            pysamAlignment, reads.SamFlags.SUPPLEMENTARY_ALIGNMENT)
        # TODO test readNumber and numberReads (nice naming guys...) once
        # we have figured out what they mean and how the map back to
        # the SAM flags 0x1, 0x40 and 0x80

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
