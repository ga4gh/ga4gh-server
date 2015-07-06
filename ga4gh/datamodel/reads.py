"""
Module responsible for translating read data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import os

import pysam

import ga4gh.protocol as protocol
import ga4gh.datamodel as datamodel
import ga4gh.exceptions as exceptions


class SamCigar(object):
    """
    Utility class for working with SAM CIGAR strings
    """
    # see http://pysam.readthedocs.org/en/latest/api.html
    # #pysam.AlignedSegment.cigartuples
    cigarStrings = [
        protocol.CigarOperation.ALIGNMENT_MATCH,
        protocol.CigarOperation.INSERT,
        protocol.CigarOperation.DELETE,
        protocol.CigarOperation.SKIP,
        protocol.CigarOperation.CLIP_SOFT,
        protocol.CigarOperation.CLIP_HARD,
        protocol.CigarOperation.PAD,
        protocol.CigarOperation.SEQUENCE_MATCH,
        protocol.CigarOperation.SEQUENCE_MISMATCH,
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
    NUMBER_READS = 0x1
    PROPER_PLACEMENT = 0x2
    READ_NUMBER_ONE = 0x40
    READ_NUMBER_TWO = 0x80
    SECONDARY_ALIGNMENT = 0x100
    FAILED_VENDOR_QUALITY_CHECKS = 0x200
    DUPLICATE_FRAGMENT = 0x400
    SUPPLEMENTARY_ALIGNMENT = 0x800

    @staticmethod
    def isFlagSet(flagAttr, flag):
        return flagAttr & flag == flag

    @staticmethod
    def setFlag(flagAttr, flag):
        flagAttr |= flag


class AbstractReadGroupSet(datamodel.DatamodelObject):
    """
    The base class of a read group set
    """
    def __init__(self, id_):
        self._id = id_
        self._readGroups = []

    def getId(self):
        # TODO move into the superclass
        return self._id

    def getReadGroups(self):
        """
        Returns the read groups in this read group set
        """
        return self._readGroups

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReadGroupSet.
        """
        readGroupSet = protocol.ReadGroupSet()
        readGroupSet.id = self._id
        readGroupSet.readGroups = [
            readGroup.toProtocolElement() for readGroup in self._readGroups]
        readGroupSet.name = None
        readGroupSet.dataset_id = None
        return readGroupSet


class SimulatedReadGroupSet(AbstractReadGroupSet):
    """
    A simulated read group set
    """
    def __init__(self, id_):
        super(SimulatedReadGroupSet, self).__init__(id_)
        read_group_id = "{}:one".format(id_)
        readGroup = SimulatedReadGroup(read_group_id)
        self._readGroups.append(readGroup)


class HtslibReadGroupSet(datamodel.PysamDatamodelMixin, AbstractReadGroupSet):
    """
    Class representing a logical collection ReadGroups.
    """
    def __init__(self, id_, dataDir):
        super(HtslibReadGroupSet, self).__init__(id_)
        self._dataDir = dataDir
        self._readGroups = []
        self._setAccessTimes(dataDir)
        self._scanDataFiles(dataDir, ["*.bam"])

    def _addDataFile(self, path):
        filename = os.path.split(path)[1]
        localId = os.path.splitext(filename)[0]
        read_group_id = "{}:{}".format(self._id, localId)
        readGroup = HtslibReadGroup(read_group_id, path)
        self._readGroups.append(readGroup)


class AbstractReadGroup(object):
    """
    Class representing a ReadGroup. A ReadGroup is all the data that's
    processed the same way by the sequencer.  There are typically 1-10
    ReadGroups in a ReadGroupSet.
    """
    def __init__(self, id_):
        self._id = id_
        now = protocol.convertDatetime(datetime.datetime.now())
        self._creationTime = now
        self._updateTime = now

    def getId(self):
        """
        Returns the id of the read group
        """
        return self._id

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReadGroup.
        """
        # TODO this is very incomplete, but we don't have the
        # implementation to fill out the rest of the fields currently
        readGroup = protocol.ReadGroup()
        readGroup.id = self._id
        readGroup.created = self._creationTime
        readGroup.updated = self._updateTime
        readGroup.dataset_id = None
        readGroup.description = None
        readGroup.experiment = None
        readGroup.info = {}
        readGroup.name = readGroup.id
        readGroup.predictedInsertSize = None
        readGroup.programs = []
        readGroup.reference_set_id = None
        readGroup.sampleId = None
        return readGroup


class SimulatedReadGroup(AbstractReadGroup):
    """
    A simulated readgroup
    """
    def __init__(self, id_):
        super(SimulatedReadGroup, self).__init__(id_)

    def getReadAlignments(self, reference_id=None, start=None, end=None):
        for i in range(2):
            alignment = self._createReadAlignment(i)
            yield alignment

    def _createReadAlignment(self, i):
        # TODO fill out a bit more
        id_ = "{}:simulated{}".format(self._id, i)
        alignment = protocol.ReadAlignment()
        alignment.aligned_quality = [1, 2, 3]
        alignment.aligned_sequence = "ACT"
        gaPosition = protocol.Position()
        gaPosition.position = 0
        gaPosition.reference_name = "whatevs"
        gaPosition.strand = protocol.Strand.POS_STRAND
        gaLinearAlignment = protocol.LinearAlignment()
        gaLinearAlignment.position = gaPosition
        alignment.alignment = gaLinearAlignment
        alignment.duplicate_fragment = False
        alignment.failed_vendor_quality_checks = False
        alignment.fragment_length = 3
        alignment.fragment_name = id_
        alignment.id = id_
        alignment.info = {}
        alignment.next_mate_position = None
        alignment.number_reads = None
        alignment.proper_placement = False
        alignment.read_group_id = self._id
        alignment.read_number = None
        alignment.secondary_alignment = False
        alignment.supplementary_alignment = False
        return alignment


class HtslibReadGroup(datamodel.PysamDatamodelMixin, AbstractReadGroup):
    """
    A readgroup based on htslib's reading of a given file
    """
    def __init__(self, id_, dataFile):
        super(HtslibReadGroup, self).__init__(id_)
        self._samFilePath = dataFile
        try:
            self._samFile = pysam.AlignmentFile(dataFile)
        except ValueError:
            raise exceptions.FileOpenFailedException(dataFile)

    def getSamFilePath(self):
        """
        Returns the file path of the sam file
        """
        return self._samFilePath

    def getReadAlignments(self, reference_id=None, start=None, end=None):
        """
        Returns an iterator over the specified reads
        """
        # TODO If reference_id is None, return against all references,
        # including unmapped reads.
        reference_name = ""
        if reference_id is not None:
            self.sanitizeGetRName(reference_id)
            reference_name = self._samFile.getrname(reference_id)
        reference_name, start, end = self.sanitizeAlignmentFileFetch(
            reference_name, start, end)
        # TODO deal with errors from htslib
        readAlignments = self._samFile.fetch(reference_name, start, end)
        for readAlignment in readAlignments:
            yield self.convertReadAlignment(readAlignment)

    def convertReadAlignment(self, read):
        """
        Convert a pysam ReadAlignment to a GA4GH ReadAlignment
        """
        # TODO fill out remaining fields
        # TODO refine in tandem with code in converters module
        ret = protocol.ReadAlignment()
        ret.aligned_quality = list(read.query_qualities)
        ret.aligned_sequence = read.query_sequence
        ret.alignment = protocol.LinearAlignment()
        ret.alignment.mapping_quality = read.mapping_quality
        ret.alignment.position = protocol.Position()
        self.sanitizeGetRName(read.reference_id)
        ret.alignment.position.reference_name = self._samFile.getrname(
            read.reference_id)
        ret.alignment.position.position = read.reference_start
        ret.alignment.position.strand = \
            protocol.Strand.POS_STRAND  # TODO fix this!
        ret.alignment.cigar = []
        for operation, length in read.cigar:
            gaCigarUnit = protocol.CigarUnit()
            gaCigarUnit.operation = SamCigar.int2ga(operation)
            gaCigarUnit.operation_length = length
            gaCigarUnit.referenceSequence = None  # TODO fix this!
            ret.alignment.cigar.append(gaCigarUnit)
        ret.duplicate_fragment = SamFlags.isFlagSet(
            read.flag, SamFlags.DUPLICATE_FRAGMENT)
        ret.failed_vendor_quality_checks = SamFlags.isFlagSet(
            read.flag, SamFlags.FAILED_VENDOR_QUALITY_CHECKS)
        ret.fragment_length = read.template_length
        ret.fragment_name = read.query_name
        ret.id = "{}:{}".format(self._id, read.query_name)
        ret.info = {key: [str(value)] for key, value in read.tags}
        ret.next_mate_position = None
        if read.next_reference_id != -1:
            ret.next_mate_position = protocol.Position()
            self.sanitizeGetRName(read.next_reference_id)
            ret.next_mate_position.reference_name = self._samFile.getrname(
                read.next_reference_id)
            ret.next_mate_position.position = read.next_reference_start
            ret.next_mate_position.strand = \
                protocol.Strand.POS_STRAND  # TODO fix this!
        # TODO Is this the correct mapping between number_reads and
        # sam flag 0x1? What about the mapping between number_reads
        # and 0x40 and 0x80?
        ret.number_reads = None
        ret.read_number = None
        if SamFlags.isFlagSet(read.flag, SamFlags.NUMBER_READS):
            ret.number_reads = 2
            if SamFlags.isFlagSet(read.flag, SamFlags.READ_NUMBER_ONE):
                ret.read_number = 0
            elif SamFlags.isFlagSet(read.flag, SamFlags.READ_NUMBER_TWO):
                ret.read_number = 1
        ret.proper_placement = SamFlags.isFlagSet(
            read.flag, SamFlags.PROPER_PLACEMENT)
        ret.read_group_id = self._id
        ret.secondary_alignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SECONDARY_ALIGNMENT)
        ret.supplementary_alignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SUPPLEMENTARY_ALIGNMENT)
        return ret
