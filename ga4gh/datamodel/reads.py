"""
Module responsible for translating read data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import random

import pysam

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.references as references
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol
import ga4gh.pb as pb


def parseMalformedBamHeader(headerDict):
    """
    Parses the (probably) intended values out of the specified
    BAM header dictionary, which is incompletely parsed by pysam.
    This is caused by some tools incorrectly using spaces instead
    of tabs as a seperator.
    """
    headerString = " ".join(
        "{}:{}".format(k, v) for k, v in headerDict.items())
    ret = {}
    for item in headerString.split():
        key, value = item.split(":", 1)
        ret[key] = value
    return ret


class SamCigar(object):
    """
    Utility class for working with SAM CIGAR strings
    """
    # see http://pysam.readthedocs.org/en/latest/api.html
    # #pysam.AlignedSegment.cigartuples
    cigarStrings = [
        protocol.CigarUnit.ALIGNMENT_MATCH,
        protocol.CigarUnit.INSERT,
        protocol.CigarUnit.DELETE,
        protocol.CigarUnit.SKIP,
        protocol.CigarUnit.CLIP_SOFT,
        protocol.CigarUnit.CLIP_HARD,
        protocol.CigarUnit.PAD,
        protocol.CigarUnit.SEQUENCE_MATCH,
        protocol.CigarUnit.SEQUENCE_MISMATCH,
    ]

    @classmethod
    def ga2int(cls, value):
        for i, cigarString in enumerate(cls.cigarStrings):
            if value == cigarString:
                return i

    @classmethod
    def int2ga(cls, value):
        return cls.cigarStrings[value]


class SamFlags(object):
    """
    Utility class for working with SAM flags
    """
    READ_PAIRED = 0x1
    READ_PROPER_PAIR = 0x2
    READ_UNMAPPED = 0x4
    MATE_UNMAPPED = 0x8
    READ_REVERSE_STRAND = 0x10
    MATE_REVERSE_STRAND = 0x20
    FIRST_IN_PAIR = 0x40
    SECOND_IN_PAIR = 0x80
    SECONDARY_ALIGNMENT = 0x100
    FAILED_QUALITY_CHECK = 0x200
    DUPLICATE_READ = 0x400
    SUPPLEMENTARY_ALIGNMENT = 0x800

    @staticmethod
    def isFlagSet(flagAttr, flag):
        return flagAttr & flag == flag

    @staticmethod
    def setFlag(flagAttr, flag):
        return flagAttr | flag


class AbstractReadGroupSet(datamodel.DatamodelObject):
    """
    The base class of a read group set
    """
    compoundIdClass = datamodel.ReadGroupSetCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractReadGroupSet, self).__init__(parentContainer, localId)
        self._readGroupIdMap = {}
        self._readGroupIds = []
        self._referenceSet = None

    def addReadGroup(self, readGroup):
        """
        Adds the specified ReadGroup to this ReadGroupSet.
        """
        id_ = readGroup.getId()
        self._readGroupIdMap[id_] = readGroup
        self._readGroupIds.append(id_)

    def getReadGroups(self):
        """
        Returns the list of ReadGroups in this ReadGroupSet.
        """
        return [self._readGroupIdMap[id_] for id_ in self._readGroupIds]

    def getReadGroup(self, id_):
        """
        Returns the ReadGroup with the specified id if it exists in this
        ReadGroupSet, or raises a ReadGroupNotFoundException otherwise.
        """
        if id_ not in self._readGroupIdMap:
            raise exceptions.ReadGroupNotFoundException(id_)
        return self._readGroupIdMap[id_]

    def getReferenceSet(self):
        """
        Returns the ReferenceSet that this ReadGroupSet is aligned to.
        """
        return self._referenceSet

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReadGroupSet.
        """
        readGroupSet = protocol.ReadGroupSet()
        readGroupSet.id = self.getId()
        readGroupSet.read_groups.extend(
            [readGroup.toProtocolElement()
                for readGroup in self.getReadGroups()]
        )
        readGroupSet.name = self.getLocalId()
        readGroupSet.dataset_id = self.getParentContainer().getId()
        stats = readGroupSet.stats
        stats.aligned_read_count = pb.int(self.getNumAlignedReads())
        stats.unaligned_read_count = pb.int(self.getNumUnalignedReads())
        return readGroupSet

    def getNumAlignedReads(self):
        """
        Return the number of aligned reads in this read group set
        """
        raise NotImplementedError()

    def getNumUnalignedReads(self):
        """
        Return the number of unaligned reads in this read group set
        """
        raise NotImplementedError()

    def getPrograms(self):
        """
        Returns an array of Programs used to generate this read group set
        """
        raise NotImplementedError()


class SimulatedReadGroupSet(AbstractReadGroupSet):
    """
    A simulated read group set
    """
    def __init__(
            self, parentContainer, localId, referenceSet, randomSeed=1,
            numReadGroups=1, numAlignments=2):
        super(SimulatedReadGroupSet, self).__init__(
            parentContainer, localId)
        self._referenceSet = referenceSet
        self._numAlignments = numAlignments
        for i in range(numReadGroups):
            localId = "rg{}".format(i)
            readGroup = SimulatedReadGroup(
                self, localId, randomSeed + i, numAlignments)
            self.addReadGroup(readGroup)

    def getNumAlignedReads(self):
        return self._numAlignments

    def getNumUnalignedReads(self):
        return 0

    def getPrograms(self):
        return []


class HtslibReadGroupSet(datamodel.PysamDatamodelMixin, AbstractReadGroupSet):
    """
    Class representing a logical collection ReadGroups.
    """
    def __init__(
            self, parentContainer, localId, samFilePath, dataRepository):
        super(HtslibReadGroupSet, self).__init__(parentContainer, localId)
        self._samFilePath = samFilePath
        samFile = self.getFileHandle(self._samFilePath)
        self._setHeaderFields(samFile)
        if 'RG' not in samFile.header or len(samFile.header['RG']) == 0:
            self._defaultReadGroup = True
            readGroup = HtslibReadGroup(self, 'default')
            self.addReadGroup(readGroup)
        else:
            self._defaultReadGroup = False
            for readGroupHeader in samFile.header['RG']:
                readGroup = HtslibReadGroup(
                    self, readGroupHeader['ID'], readGroupHeader)
                self.addReadGroup(readGroup)
        # Find the reference set name (if there is one) by looking at
        # the BAM headers.
        referenceSetName = None
        for referenceInfo in samFile.header['SQ']:
            if 'AS' not in referenceInfo:
                infoDict = parseMalformedBamHeader(referenceInfo)
            else:
                infoDict = referenceInfo
            name = infoDict.get('AS', references.DEFAULT_REFERENCESET_NAME)
            if referenceSetName is None:
                referenceSetName = name
            elif referenceSetName != name:
                raise exceptions.MultipleReferenceSetsInReadGroupSet(
                    samFilePath, name, referenceSetName)
        self._referenceSet = None
        if referenceSetName is not None:
            self._referenceSet = dataRepository.getReferenceSetByName(
                referenceSetName)
            # TODO verify that the references in the BAM file exist
            # in the reference set. Otherwise, we won't be able to
            # query for them.

    def _setHeaderFields(self, samFile):
        programs = []
        if 'PG' in samFile.header:
            htslibPrograms = samFile.header['PG']
            for htslibProgram in htslibPrograms:
                program = protocol.ReadGroup.Program()
                program.id = htslibProgram['ID']
                program.command_line = htslibProgram.get(
                    'CL', pb.DEFAULT_STRING)
                program.name = htslibProgram.get('PN', pb.DEFAULT_STRING)
                program.prev_program_id = htslibProgram.get(
                    'PP', pb.DEFAULT_STRING)
                program.version = htslibProgram.get('VN', pb.DEFAULT_STRING)
                programs.append(program)
        self._programs = programs

    def openFile(self, dataFile):
        return pysam.AlignmentFile(dataFile)

    def getSamFilePath(self):
        """
        Returns the file path of the sam file
        """
        return self._samFilePath

    def isUsingDefaultReadGroup(self):
        """
        Returns whether the readGroupSet is using a default read group
        """
        return self._defaultReadGroup

    def getNumAlignedReads(self):
        samFile = self.getFileHandle(self._samFilePath)
        return samFile.mapped

    def getNumUnalignedReads(self):
        samFile = self.getFileHandle(self._samFilePath)
        return samFile.unmapped

    def getPrograms(self):
        return self._programs


class AbstractReadGroup(datamodel.DatamodelObject):
    """
    Class representing a ReadGroup. A ReadGroup is all the data that's
    processed the same way by the sequencer.  There are typically 1-10
    ReadGroups in a ReadGroupSet.
    """
    compoundIdClass = datamodel.ReadGroupCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractReadGroup, self).__init__(parentContainer, localId)
        datetimeNow = datetime.datetime.now()
        now = protocol.convertDatetime(datetimeNow)
        self._iso8601 = datetimeNow.strftime("%Y-%m-%dT%H:%M:%SZ")
        self._creationTime = now
        self._updateTime = now

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReadGroup.
        """
        # TODO this is very incomplete, but we don't have the
        # implementation to fill out the rest of the fields currently
        readGroup = protocol.ReadGroup()
        readGroup.id = self.getId()
        readGroup.created = self._creationTime
        readGroup.updated = self._updateTime
        dataset = self.getParentContainer().getParentContainer()
        readGroup.dataset_id = dataset.getId()
        readGroup.name = self.getLocalId()
        readGroup.predicted_insert_size = pb.int(self.getPredictedInsertSize())
        referenceSet = self._parentContainer.getReferenceSet()
        readGroup.sample_id = pb.string(self.getSampleId())
        if referenceSet is not None:
            readGroup.reference_set_id = referenceSet.getId()
        stats = readGroup.stats
        stats.aligned_read_count = self.getNumAlignedReads()
        stats.unaligned_read_count = self.getNumUnalignedReads()
        readGroup.programs.extend(self.getPrograms())
        readGroup.description = pb.string(self.getDescription())
        experiment = readGroup.experiment
        experiment.id = self.getExperimentId()
        experiment.instrument_model = pb.string(self.getInstrumentModel())
        experiment.sequencing_center = pb.string(self.getSequencingCenter())
        experiment.description = pb.string(self.getExperimentDescription())
        experiment.library = pb.string(self.getLibrary())
        experiment.platform_unit = pb.string(self.getPlatformUnit())
        experiment.message_create_time = self._iso8601
        experiment.message_update_time = self._iso8601
        experiment.run_time = pb.string(self.getRunTime())
        return readGroup

    def getReadAlignmentId(self, gaAlignment):
        """
        Returns a string ID suitable for use in the specified GA
        ReadAlignment object in this ReadGroup.
        """
        compoundId = datamodel.ReadAlignmentCompoundId(
            self.getCompoundId(), gaAlignment.fragment_name)
        return str(compoundId)

    def getNumAlignedReads(self):
        """
        Return the number of aligned reads in the read group
        """
        raise NotImplementedError()

    def getNumUnalignedReads(self):
        """
        Return the number of unaligned reads in the read group
        """
        raise NotImplementedError()

    def getPrograms(self):
        """
        Returns an array of Programs used to generate this read group
        """
        raise NotImplementedError()

    def getDescription(self):
        """
        Returns a description of this read group
        """
        raise NotImplementedError()

    def getSampleId(self):
        """
        Returns the sample id of the read group
        """
        raise NotImplementedError()

    def getPredictedInsertSize(self):
        """
        Returns the predicted insert size of the read group
        """
        raise NotImplementedError()

    def getInstrumentModel(self):
        """
        Returns the instrument model used for this experiment
        """
        raise NotImplementedError()

    def getSequencingCenter(self):
        """
        Returns the sequencing center used for this experiment
        """
        raise NotImplementedError()

    def getExperimentDescription(self):
        """
        Returns the description of this read group
        """
        raise NotImplementedError()

    def getLibrary(self):
        """
        Returns the name of the library used in this experiment
        """
        raise NotImplementedError()

    def getPlatformUnit(self):
        """
        Returns the platform unit used in this experiment
        """
        raise NotImplementedError()

    def getRunTime(self):
        """
        Returns the time at which the experiment was performed
        """
        raise NotImplementedError()

    def getExperimentId(self):
        """
        Returns the id of the experiment used for this read group
        """
        return str(datamodel.ExperimentCompoundId(
            self.getCompoundId(), 'experiment'))


class SimulatedReadGroup(AbstractReadGroup):
    """
    A simulated readgroup
    """
    def __init__(self, parentContainer, localId, randomSeed, numAlignments=2):
        super(SimulatedReadGroup, self).__init__(parentContainer, localId)
        self._randomSeed = randomSeed

    def getReadAlignments(self, referenceId=None, start=None, end=None):
        rng = random.Random(self._randomSeed)

        # We seed reads with sequential seeds starting from here. We hope no
        # ranges for two different simulated read groups ever overlap (because
        # then we'd start seeing identical reads in the two groups.)
        read_seed_start = rng.getrandbits(64)

        for i in range(self.getNumAlignedReads()):
            seed = read_seed_start + i
            yield self._createReadAlignment(i, seed)

    def _createReadAlignment(self, i, seed):
        # TODO fill out a bit more
        rng = random.Random(seed)
        alignment = protocol.ReadAlignment()
        alignment.fragment_length = rng.randint(10, 100)
        alignment.aligned_sequence = ""
        for i in range(alignment.fragment_length):
            # TODO: are these reasonable quality values?
            alignment.aligned_quality.append(rng.randint(1, 20))
            alignment.aligned_sequence += rng.choice("ACGT")

        gaLinearAlignment = alignment.alignment
        gaPosition = gaLinearAlignment.position
        gaPosition.position = 0
        gaPosition.reference_name = "NotImplemented"
        gaPosition.strand = protocol.POS_STRAND
        alignment.duplicate_fragment = False
        alignment.failed_vendor_quality_checks = False

        alignment.fragment_name = "simulated{}".format(i)
        alignment.info.clear()
        alignment.next_mate_position.Clear()
        alignment.number_reads = 0
        alignment.proper_placement = False
        alignment.read_group_id = self.getId()
        alignment.read_number = -1
        alignment.secondary_alignment = False
        alignment.supplementary_alignment = False
        alignment.id = self.getReadAlignmentId(alignment)
        return alignment

    def getNumAlignedReads(self):
        return self._parentContainer.getNumAlignedReads()

    def getNumUnalignedReads(self):
        return 0

    def getPrograms(self):
        return []

    def getDescription(self):
        return None

    def getSampleId(self):
        return 'sampleId'

    def getPredictedInsertSize(self):
        return 0

    def getInstrumentModel(self):
        return None

    def getSequencingCenter(self):
        return None

    def getExperimentDescription(self):
        return None

    def getLibrary(self):
        return None

    def getPlatformUnit(self):
        return None

    def getRunTime(self):
        return None


class HtslibReadGroup(datamodel.PysamDatamodelMixin, AbstractReadGroup):
    """
    A readgroup based on htslib's reading of a given file
    """
    def __init__(self, parentContainer, localId, readGroupHeader=None):
        super(HtslibReadGroup, self).__init__(parentContainer, localId)
        self._parentSamFilePath = parentContainer.getSamFilePath()
        self._filterReads = not parentContainer.isUsingDefaultReadGroup()
        self._sampleId = None
        self._description = None
        self._predictedInsertSize = None
        self._instrumentModel = None
        self._sequencingCenter = None
        self._experimentDescription = None
        self._library = None
        self._platformUnit = None
        self._runTime = None
        if readGroupHeader is not None:
            self._sampleId = readGroupHeader.get('SM', None)
            self._description = readGroupHeader.get('DS', None)
            if 'PI' in readGroupHeader:
                self._predictedInsertSize = int(readGroupHeader['PI'])
            self._instrumentModel = readGroupHeader.get('PL', None)
            self._sequencingCenter = readGroupHeader.get('CN', None)
            self._experimentDescription = readGroupHeader.get('DS', None)
            self._library = readGroupHeader.get('LB', None)
            self._platformUnit = readGroupHeader.get('PU', None)
            self._runTime = readGroupHeader.get('DT', None)

    def getSamFilePath(self):
        return self._parentSamFilePath

    def getReadAlignments(self, reference, start=None, end=None):
        """
        Returns an iterator over the specified reads
        """
        # TODO If reference is None, return against all references,
        # including unmapped reads.
        samFile = self._parentContainer.getFileHandle(self._parentSamFilePath)
        referenceName = reference.getLocalId().encode()
        # TODO deal with errors from htslib
        start, end = self.sanitizeAlignmentFileFetch(start, end)
        readAlignments = samFile.fetch(referenceName, start, end)
        if self._filterReads:
            for readAlignment in readAlignments:
                tags = dict(readAlignment.tags)
                if 'RG' in tags and tags['RG'] == self._localId:
                    yield self.convertReadAlignment(readAlignment)
        else:
            for readAlignment in readAlignments:
                yield self.convertReadAlignment(readAlignment)

    def convertReadAlignment(self, read):
        """
        Convert a pysam ReadAlignment to a GA4GH ReadAlignment
        """
        samFile = self._parentContainer.getFileHandle(
            self._parentSamFilePath)
        # TODO fill out remaining fields
        # TODO refine in tandem with code in converters module
        ret = protocol.ReadAlignment()
        if read.query_qualities is None:
            pass  # ret.aligned_quality.clear()
        else:
            ret.aligned_quality.extend(read.query_qualities)
        ret.aligned_sequence = read.query_sequence
        if SamFlags.isFlagSet(read.flag, SamFlags.READ_UNMAPPED):
            ret.alignment.Clear()
        else:
            ret.alignment.mapping_quality = read.mapping_quality
            ret.alignment.position.reference_name = samFile.getrname(
                read.reference_id)
            ret.alignment.position.position = read.reference_start
            ret.alignment.position.strand = protocol.POS_STRAND
            if SamFlags.isFlagSet(read.flag, SamFlags.READ_REVERSE_STRAND):
                ret.alignment.position.strand = protocol.NEG_STRAND
            for operation, length in read.cigar:
                gaCigarUnit = protocol.CigarUnit()
                gaCigarUnit.operation = SamCigar.int2ga(operation)
                gaCigarUnit.operation_length = length
                gaCigarUnit.reference_sequence = ""  # TODO fix this!
                ret.alignment.cigar.extend([gaCigarUnit])
        ret.duplicate_fragment = SamFlags.isFlagSet(
            read.flag, SamFlags.DUPLICATE_READ)
        ret.failed_vendor_quality_checks = SamFlags.isFlagSet(
            read.flag, SamFlags.FAILED_QUALITY_CHECK)
        ret.fragment_length = read.template_length
        ret.fragment_name = read.query_name
        for key, value in read.tags:
            ret.info[key].values.add(string_value=str(value))
        if SamFlags.isFlagSet(read.flag, SamFlags.MATE_UNMAPPED):
            pass
        else:
            if read.next_reference_id != -1:
                ret.next_mate_position.reference_name = samFile.getrname(
                    read.next_reference_id)
            else:
                ret.next_mate_position.reference_name = ""
            ret.next_mate_position.position = read.next_reference_start
            ret.next_mate_position.strand = protocol.POS_STRAND
            if SamFlags.isFlagSet(read.flag, SamFlags.MATE_REVERSE_STRAND):
                ret.next_mate_position.strand = protocol.NEG_STRAND
        if SamFlags.isFlagSet(read.flag, SamFlags.READ_PAIRED):
            ret.number_reads = 2
        else:
            ret.number_reads = 1
        ret.read_number = -1
        if SamFlags.isFlagSet(read.flag, SamFlags.FIRST_IN_PAIR):
            if SamFlags.isFlagSet(read.flag, SamFlags.SECOND_IN_PAIR):
                ret.read_number = 2
            else:
                ret.read_number = 0
        elif SamFlags.isFlagSet(read.flag, SamFlags.SECOND_IN_PAIR):
            ret.read_number = 1
        ret.proper_placement = SamFlags.isFlagSet(
            read.flag, SamFlags.READ_PROPER_PAIR)
        ret.read_group_id = self.getId()
        ret.secondary_alignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SECONDARY_ALIGNMENT)
        ret.supplementary_alignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SUPPLEMENTARY_ALIGNMENT)
        ret.id = self.getReadAlignmentId(ret)
        return ret

    def getNumAlignedReads(self):
        return -1  # TODO populate with metadata

    def getNumUnalignedReads(self):
        return -1  # TODO populate with metadata

    def getPrograms(self):
        return self._parentContainer.getPrograms()

    def getDescription(self):
        return self._description

    def getSampleId(self):
        return self._sampleId

    def getPredictedInsertSize(self):
        return self._predictedInsertSize

    def getInstrumentModel(self):
        return self._instrumentModel

    def getSequencingCenter(self):
        return self._sequencingCenter

    def getExperimentDescription(self):
        return self._experimentDescription

    def getLibrary(self):
        return self._library

    def getPlatformUnit(self):
        return self._platformUnit

    def getRunTime(self):
        return self._runTime
