"""
Data-driven tests for reads
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import os

import ga4gh.server.backend as backend
import ga4gh.server.datamodel as datamodel
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.datamodel.reads as reads
import ga4gh.server.datamodel.references as references
import ga4gh.server.datarepo as datarepo
import tests.datadriven as datadriven
import tests.paths as paths

import ga4gh.common.utils as utils
import ga4gh.schemas.protocol as protocol

import pysam


def testReads():
    testDataDir = os.path.join(paths.testDataDir, "datasets/dataset1/reads")
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
        self.sampleName = None
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
            self.sampleName = readGroupHeader.get('SM', None)
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
        self._backend = backend.Backend(datarepo.AbstractDataRepository())
        self._referenceSet = None
        self._dataset = datasets.Dataset("ds")
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
        self._backend.getDataRepository().addReferenceSet(self._referenceSet)

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
        readGroupSet = reads.HtslibReadGroupSet(self._dataset, localId)
        readGroupSet.populateFromFile(dataPath)
        return readGroupSet

    def getProtocolClass(self):
        return protocol.ReadGroupSet

    def testSampleNameEtc(self):
        # test that sampleId and other misc fields are set correctly
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            readGroupInfo = self._readGroupInfos[readGroup.getLocalId()]
            gaReadGroup = readGroup.toProtocolElement()
            self.assertEqual(
                readGroupInfo.sampleName,
                gaReadGroup.sample_name)
            self.assertEqual(
                readGroupInfo.predictedInsertSize,
                gaReadGroup.predicted_insert_size)
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
                gaReadGroup.experiment.instrument_model)
            self.assertEqual(
                readGroupInfo.sequencingCenter,
                gaReadGroup.experiment.sequencing_center)
            self.assertEqual(
                readGroupInfo.experimentDescription,
                gaReadGroup.experiment.description)
            self.assertEqual(
                readGroupInfo.library,
                gaReadGroup.experiment.library)
            self.assertEqual(
                readGroupInfo.platformUnit,
                gaReadGroup.experiment.platform_unit)
            self.assertEqual(
                readGroupInfo.runTime,
                gaReadGroup.experiment.run_time)

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
                    gaProgram.command_line, htslibProgram.get('CL', None))
                self.assertEqual(
                    gaProgram.name, htslibProgram.get('PN', None))
                self.assertEqual(
                    gaProgram.prev_program_id, htslibProgram.get('PP', None))
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
            gaReadGroupSet.stats.aligned_read_count,
            readGroupSetInfo.numAlignedReads)
        self.assertEqual(
            gaReadGroupSet.stats.unaligned_read_count,
            readGroupSetInfo.numUnalignedReads)
        for readGroup in readGroupSet.getReadGroups():
            gaReadGroup = readGroup.toProtocolElement()
            self.assertEqual(
                readGroup.getNumAlignedReads(), -1)
            self.assertEqual(
                readGroup.getNumUnalignedReads(), -1)
            self.assertEqual(
                gaReadGroup.stats.aligned_read_count, -1)
            self.assertEqual(
                gaReadGroup.stats.unaligned_read_count, -1)

    def testValidateObjects(self):
        # test that validation works on read groups and reads
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.getReadGroups():
            self.assertIsInstance(
                readGroup.toProtocolElement(), protocol.ReadGroup)
            for reference in self._referenceSet.getReferences():
                for gaAlignment in readGroup.getReadAlignments(reference):
                    self.assertIsInstance(
                        gaAlignment, protocol.ReadAlignment)

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
                positions = [
                    read.alignment.position.position for read in alignments
                    if read.alignment is not None]
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

    def getDictFromMessageMap(self, messageMap):
        return dict([
            (k, [protocol.getValueFromValue(x) for x in v.values])
            for (k, v) in messageMap._values.items()])

    def assertAlignmentsEqual(self, gaAlignment, pysamAlignment,
                              readGroupInfo):
        if pysamAlignment.query_qualities is None:
            self.assertEqual(gaAlignment.aligned_quality, [])
        else:
            self.assertEqual(
                gaAlignment.aligned_quality,
                list(pysamAlignment.query_qualities))
        self.assertEqual(
            gaAlignment.aligned_sequence,
            pysamAlignment.query_sequence)
        if reads.SamFlags.isFlagSet(
                pysamAlignment.flag, reads.SamFlags.READ_UNMAPPED):
            self.assertEqual(0, gaAlignment.alignment.ByteSize())
        else:
            self.assertEqual(
                gaAlignment.alignment.mapping_quality,
                pysamAlignment.mapping_quality)
            self.assertEqual(
                gaAlignment.alignment.position.reference_name,
                readGroupInfo.samFile.getrname(pysamAlignment.reference_id))
            self.assertEqual(
                gaAlignment.alignment.position.position,
                pysamAlignment.reference_start)
            # TODO test reverseStrand on position and on
            # nextMatePosition once it has been implemented.
            self.assertCigarEqual(
                gaAlignment.alignment.cigar,
                pysamAlignment.cigar)
        self.assertFlag(
            gaAlignment.duplicate_fragment,
            pysamAlignment, reads.SamFlags.DUPLICATE_READ)
        self.assertFlag(
            gaAlignment.failed_vendor_quality_checks,
            pysamAlignment, reads.SamFlags.FAILED_QUALITY_CHECK)
        self.assertEqual(
            gaAlignment.fragment_length,
            pysamAlignment.template_length)
        self.assertEqual(
            gaAlignment.fragment_name,
            pysamAlignment.query_name)
        compoundId = datamodel.ReadAlignmentCompoundId(
            self._gaObject.getCompoundId(),
            pysamAlignment.query_name)
        self.assertEqual(gaAlignment.id, str(compoundId))
        ret = protocol.ReadAlignment()
        for key, value in pysamAlignment.tags:
            protocol.setAttribute(ret.attributes.attr[key].values, value)
        self.assertEqual(
            protocol.toJson(gaAlignment.attributes),
            protocol.toJson(ret.attributes))
        if reads.SamFlags.isFlagSet(
                pysamAlignment.flag, reads.SamFlags.MATE_UNMAPPED):
            self.assertEqual(0, gaAlignment.next_mate_position.ByteSize())
        else:
            self.assertEqual(
                gaAlignment.next_mate_position.position,
                pysamAlignment.next_reference_start)
            if pysamAlignment.next_reference_id != -1:
                self.assertEqual(
                    gaAlignment.next_mate_position.reference_name,
                    readGroupInfo.samFile.getrname(
                        pysamAlignment.next_reference_id))
            else:
                self.assertEqual(
                    gaAlignment.next_mate_position.reference_name, "")
        if gaAlignment.number_reads == 1:
            self.assertFlag(
                False, pysamAlignment, reads.SamFlags.READ_PAIRED)
        elif gaAlignment.number_reads == 2:
            self.assertFlag(
                True, pysamAlignment, reads.SamFlags.READ_PAIRED)
        else:
            # we shouldn't be setting numberReads to anything else
            self.assertTrue(False)
        if gaAlignment.read_number is -1:
            self.assertFlag(
                False, pysamAlignment, reads.SamFlags.FIRST_IN_PAIR)
            self.assertFlag(
                False, pysamAlignment, reads.SamFlags.SECOND_IN_PAIR)
        elif gaAlignment.read_number == 0:
            self.assertFlag(
                True, pysamAlignment, reads.SamFlags.FIRST_IN_PAIR)
            self.assertFlag(
                False, pysamAlignment, reads.SamFlags.SECOND_IN_PAIR)
        elif gaAlignment.read_number == 1:
            self.assertFlag(
                False, pysamAlignment, reads.SamFlags.FIRST_IN_PAIR)
            self.assertFlag(
                True, pysamAlignment, reads.SamFlags.SECOND_IN_PAIR)
        elif gaAlignment.read_number == 2:
            self.assertFlag(
                True, pysamAlignment, reads.SamFlags.FIRST_IN_PAIR)
            self.assertFlag(
                True, pysamAlignment, reads.SamFlags.SECOND_IN_PAIR)
        else:
            # we shouldn't be setting readNumber to anything else
            self.assertTrue(False)
        self.assertFlag(
            not gaAlignment.improper_placement,
            pysamAlignment, reads.SamFlags.READ_PROPER_PAIR)
        self.assertEqual(
            gaAlignment.read_group_id,
            readGroupInfo.id)
        self.assertFlag(
            gaAlignment.secondary_alignment,
            pysamAlignment, reads.SamFlags.SECONDARY_ALIGNMENT)
        self.assertFlag(
            gaAlignment.supplementary_alignment,
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
                gaCigarUnit.operation_length, length)
