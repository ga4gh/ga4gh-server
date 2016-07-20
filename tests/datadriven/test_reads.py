"""
Data-driven tests for reads
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import os

import ga4gh.registry as registry
import ga4gh.datasource.htslib as htslib
import ga4gh.protocol as protocol
import tests.datadriven as datadriven
import tests.utils as utils
import tests.paths as paths

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


def get_reference_by_name(reference_set, name):
    """
    Returns the reference within the specified reference set with the
    specified name.

    TODO once the data driven tests code paths have been refactored
    to use the external interface, remove this method.
    """
    ret = None
    for reference in reference_set.references:
        if reference.name == name:
            ret = reference
            break
    return ret


class ReadGroupInfo(object):
    """
    Container class for information about a read group
    """
    def __init__(self, gaReadGroupSet, samFile, readGroupName):
        self.id = None
        for read_group in gaReadGroupSet.read_groups:
            if read_group.name == readGroupName:
                self.id = str(read_group.id)
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
        self._registry_db = registry.RegistryDb(paths.testDataRepoUrl)
        self._registry_db.open()
        super(ReadGroupSetTest, self).__init__(localId, dataPath)
        self._readGroupInfos = {}
        self._readGroupSetInfo = None
        self._samFile = pysam.AlignmentFile(dataPath)
        self._readAlignmentInfo()

    def _readAlignmentInfo(self):
        self._readGroupSetInfo = ReadGroupSetInfo(self._samFile)
        if 'RG' in self._samFile.header:
            readGroupHeaders = self._samFile.header['RG']
            readGroupNames = [
                readGroupHeader['ID'] for readGroupHeader
                in readGroupHeaders]
        else:
            readGroupNames = ['Default_RG']
        for readGroupName in readGroupNames:
            readGroupInfo = ReadGroupInfo(
                self._gaObject, self._samFile, readGroupName)
            self._readGroupInfos[readGroupName] = readGroupInfo

    def getDataModelInstance(self, localId, dataPath):
        # We assume the files have been added with the default name rules.
        name = localId.split(".")[0]
        read_group_set = self._registry_db.get_read_group_set_by_name(
            "dataset1", name)
        self.assertTrue(os.path.samefile(read_group_set.data_url, dataPath))
        return read_group_set

    def getProtocolClass(self):
        return protocol.ReadGroupSet

    def testSampleNameEtc(self):
        # test that sampleId and other misc fields are set correctly
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.read_groups:
            readGroupInfo = self._readGroupInfos[readGroup.name]
            gaReadGroup = readGroup.get_protobuf()
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
        for readGroup in readGroupSet.read_groups:
            readGroupInfo = self._readGroupInfos[readGroup.name]
            gaReadGroup = readGroup.get_protobuf()
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
        readGroupSet = self._gaObject.get_protobuf()
        for readGroup in readGroupSet.read_groups:
            readGroupInfo = self._readGroupInfos[readGroup.name]
            gaPrograms = readGroup.programs
            htslibPrograms = readGroupInfo.programs
            for gaProgram, htslibProgram in utils.zipLists(
                    gaPrograms, htslibPrograms):
                self.assertEqual(
                    gaProgram.command_line, htslibProgram.get('CL', ""))
                self.assertEqual(
                    gaProgram.name, htslibProgram.get('PN', ""))
                # TODO fix this when we have proper ID mapping for programs
                # self.assertEqual(
                #     gaProgram.prev_program_id, htslibProgram.get('PP', ""))
                self.assertEqual(
                    gaProgram.version, htslibProgram.get('VN', ""))

    def testReadGroupStats(self):
        # test that the stats attrs are populated correctly
        readGroupSet = self._gaObject
        gaReadGroupSet = readGroupSet.get_protobuf()
        readGroupSetInfo = self._readGroupSetInfo
        self.assertEqual(
            readGroupSet.read_stats.aligned_read_count,
            readGroupSetInfo.numAlignedReads)
        self.assertEqual(
            readGroupSet.read_stats.unaligned_read_count,
            readGroupSetInfo.numUnalignedReads)
        self.assertEqual(
            readGroupSet.read_stats.base_count, -1)
        self.assertEqual(
            gaReadGroupSet.stats.aligned_read_count,
            readGroupSetInfo.numAlignedReads)
        self.assertEqual(
            gaReadGroupSet.stats.unaligned_read_count,
            readGroupSetInfo.numUnalignedReads)
        self.assertEqual(
            gaReadGroupSet.stats.base_count, -1)
        for readGroup in readGroupSet.read_groups:
            gaReadGroup = readGroup.get_protobuf()
            self.assertEqual(
                readGroup.read_stats.aligned_read_count, -1)
            self.assertEqual(
                readGroup.read_stats.unaligned_read_count, -1)
            self.assertEqual(
                readGroup.read_stats.base_count, -1)
            self.assertEqual(
                gaReadGroup.stats.aligned_read_count, -1)
            self.assertEqual(
                gaReadGroup.stats.unaligned_read_count, -1)
            self.assertEqual(
                gaReadGroup.stats.base_count, -1)

    def testGetReadAlignmentsRefId(self):
        # test that searching with a reference id succeeds
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.read_groups:
            readGroupInfo = self._readGroupInfos[readGroup.name]
            for name, alignments in readGroupInfo.mappedReads.items():
                reference = get_reference_by_name(
                    readGroupSet.reference_set, name)
                gaAlignments = list(
                    readGroupSet.get_read_alignments(
                        reference, [readGroup]))
                self.assertAlignmentListsEqual(
                    gaAlignments, alignments, readGroupInfo)

    def testGetReadAlignmentsStartEnd(self):
        # test that searching with start and end coords succeeds
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.read_groups:
            readGroupInfo = self._readGroupInfos[readGroup.name]
            for name, alignments, in readGroupInfo.mappedReads.items():
                bigNumThatPysamWontChokeOn = 2**30
                reference = get_reference_by_name(
                    readGroupSet.reference_set, name)
                gaAlignments = list(readGroupSet.get_read_alignments(
                    reference, [readGroup], 0, bigNumThatPysamWontChokeOn))
                self.assertAlignmentListsEqual(
                    gaAlignments, alignments, readGroupInfo)

    def testGetReadAlignmentSearchRanges(self):
        # test that various range searches work
        readGroupSet = self._gaObject
        for readGroup in readGroupSet.read_groups:
            readGroupInfo = self._readGroupInfos[readGroup.name]
            for name in readGroupInfo.mappedReads.keys():
                reference = get_reference_by_name(
                    readGroupSet.reference_set, name)
                alignments = list(readGroupSet.get_read_alignments(
                    reference, [readGroup]))
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
                    readGroupSet, readGroup, reference, begin, end + 1, length)
                self.assertGetReadAlignmentsRangeResult(
                    readGroupSet, readGroup, reference, begin, end, length - 1)
                self.assertGetReadAlignmentsRangeResult(
                    readGroupSet, readGroup, reference, begin, begin, 0)

    def assertGetReadAlignmentsRangeResult(
            self, readGroupSet, readGroup, reference, start, end, result):
        alignments = list(
            readGroupSet.get_read_alignments(
                reference, [readGroup], start, end))
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
        if htslib.SamFlags.isFlagSet(
                pysamAlignment.flag, htslib.SamFlags.READ_UNMAPPED):
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
            pysamAlignment, htslib.SamFlags.DUPLICATE_READ)
        self.assertFlag(
            gaAlignment.failed_vendor_quality_checks,
            pysamAlignment, htslib.SamFlags.FAILED_QUALITY_CHECK)
        self.assertEqual(
            gaAlignment.fragment_length,
            pysamAlignment.template_length)
        self.assertEqual(
            gaAlignment.fragment_name,
            pysamAlignment.query_name)
        compoundId = htslib.HtslibReadAlignmentCompoundId(
            None, gaAlignment.read_group_id, pysamAlignment.query_name)
        self.assertEqual(gaAlignment.id, str(compoundId))
        tags1 = self.getDictFromMessageMap(gaAlignment.info)
        tags2 = {key: [str(value)] for key, value in pysamAlignment.tags}
        if 'RG' in tags2:
            del tags2['RG']
        self.assertEqual(tags1, tags2)
        if htslib.SamFlags.isFlagSet(
                pysamAlignment.flag, htslib.SamFlags.MATE_UNMAPPED):
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
                False, pysamAlignment, htslib.SamFlags.READ_PAIRED)
        elif gaAlignment.number_reads == 2:
            self.assertFlag(
                True, pysamAlignment, htslib.SamFlags.READ_PAIRED)
        else:
            # we shouldn't be setting numberReads to anything else
            self.assertTrue(False)
        if gaAlignment.read_number is -1:
            self.assertFlag(
                False, pysamAlignment, htslib.SamFlags.FIRST_IN_PAIR)
            self.assertFlag(
                False, pysamAlignment, htslib.SamFlags.SECOND_IN_PAIR)
        elif gaAlignment.read_number == 0:
            self.assertFlag(
                True, pysamAlignment, htslib.SamFlags.FIRST_IN_PAIR)
            self.assertFlag(
                False, pysamAlignment, htslib.SamFlags.SECOND_IN_PAIR)
        elif gaAlignment.read_number == 1:
            self.assertFlag(
                False, pysamAlignment, htslib.SamFlags.FIRST_IN_PAIR)
            self.assertFlag(
                True, pysamAlignment, htslib.SamFlags.SECOND_IN_PAIR)
        elif gaAlignment.read_number == 2:
            self.assertFlag(
                True, pysamAlignment, htslib.SamFlags.FIRST_IN_PAIR)
            self.assertFlag(
                True, pysamAlignment, htslib.SamFlags.SECOND_IN_PAIR)
        else:
            # we shouldn't be setting readNumber to anything else
            self.assertTrue(False)
        self.assertFlag(
            not gaAlignment.improper_placement,
            pysamAlignment, htslib.SamFlags.READ_PROPER_PAIR)
        self.assertEqual(gaAlignment.read_group_id, readGroupInfo.id)
        self.assertFlag(
            gaAlignment.secondary_alignment,
            pysamAlignment, htslib.SamFlags.SECONDARY_ALIGNMENT)
        self.assertFlag(
            gaAlignment.supplementary_alignment,
            pysamAlignment, htslib.SamFlags.SUPPLEMENTARY_ALIGNMENT)

    def assertFlag(self, gaAlignmentAttr, pysamAlignment, mask):
        flagSet = htslib.SamFlags.isFlagSet(pysamAlignment.flag, mask)
        self.assertEqual(gaAlignmentAttr, flagSet)

    def assertCigarEqual(self, gaCigar, pysamCigar):
        self.assertEqual(len(gaCigar), len(pysamCigar))
        for i, gaCigarUnit in enumerate(gaCigar):
            operation, length = pysamCigar[i]
            gaCigarUnitOperation = htslib.SamCigar.ga2int(
                gaCigarUnit.operation)
            self.assertEqual(
                gaCigarUnitOperation, operation)
            self.assertEqual(
                gaCigarUnit.operation_length, length)
