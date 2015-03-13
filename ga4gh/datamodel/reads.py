"""
Module responsible for translating read data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import os
import glob

import pysam

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


class SamCigar(object):
    """
    Utility class for working with SAM CIGAR strings
    """
    # see http://pysam.readthedocs.org/en/latest/api.html
    # #pysam.AlignedSegment.cigartuples
    cigarStrings = [
        protocol.GACigarOperation.ALIGNMENT_MATCH,
        protocol.GACigarOperation.INSERT,
        protocol.GACigarOperation.DELETE,
        protocol.GACigarOperation.SKIP,
        protocol.GACigarOperation.CLIP_SOFT,
        protocol.GACigarOperation.CLIP_HARD,
        protocol.GACigarOperation.PAD,
        protocol.GACigarOperation.SEQUENCE_MATCH,
        protocol.GACigarOperation.SEQUENCE_MISMATCH,
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


class ReadGroupSet(object):
    """
    Class representing a logical collection ReadGroups.
    """
    def __init__(self, id_, dataDir):
        self._id = id_
        self._dataDir = dataDir
        self._readGroups = []
        # we only support BAM files right now;
        # SAM files (which can't have index files) would be too slow
        pattern = "*.bam"
        for path in glob.glob(os.path.join(self._dataDir, pattern)):
            filename = os.path.basename(path)
            localId = os.path.splitext(filename)[0]
            readGroupId = "{}:{}".format(self._id, localId)
            readGroup = ReadGroup(readGroupId, path)
            self._readGroups.append(readGroup)

    def getReadGroups(self):
        """
        Returns the read groups in this read group set
        """
        return self._readGroups

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReadGroupSet.
        """
        readGroupSet = protocol.GAReadGroupSet()
        readGroupSet.id = self._id
        readGroupSet.readGroups = [
            readGroup.toProtocolElement() for readGroup in self._readGroups]
        readGroupSet.name = None
        readGroupSet.datasetId = None
        return readGroupSet


class ReadGroup(object):
    """
    Class representing a ReadGroup. A ReadGroup is all the data that's
    processed the same way by the sequencer.  There are typically 1-10
    ReadGroups in a ReadGroupSet.
    """
    def __init__(self, id_, dataFile):
        self._id = id_
        self._samFilePath = dataFile
        self._samFile = pysam.AlignmentFile(dataFile)

    def getId(self):
        """
        Returns the id of the read group
        """
        return self._id

    def getSamFilePath(self):
        """
        Returns the file path of the sam file
        """
        return self._samFilePath

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReadGroup.
        """
        # TODO this is very incomplete, but we don't have the
        # implementation to fill out the rest of the fields currently
        now = protocol.convertDatetime(datetime.datetime.now())
        readGroup = protocol.GAReadGroup()
        readGroup.id = self._id
        readGroup.created = now
        readGroup.updated = now
        readGroup.datasetId = None
        readGroup.description = None
        readGroup.experiment = None
        readGroup.info = {}
        readGroup.name = readGroup.id
        readGroup.predictedInsertSize = None
        readGroup.programs = []
        readGroup.referenceSetId = None
        readGroup.sampleId = None
        return readGroup

    def getReadAlignments(self, referenceName=None, referenceId=None,
                          start=None, end=None):
        """
        Returns an iterator over the specified reads
        """
        if referenceName is not None and referenceId is not None:
            raise exceptions.BadReadsSearchRequestBothRefs()
        if referenceId is not None:
            referenceName = self._samFile.getrname(referenceId)
        # TODO deal with errors from htslib
        readAlignments = self._samFile.fetch(referenceName, start, end)
        for readAlignment in readAlignments:
            yield self.convertReadAlignment(readAlignment)

    def convertReadAlignment(self, read):
        """
        Convert a pysam ReadAlignment to a GAReadAlignment
        """
        # TODO fill out remaining fields
        # TODO refine in tandem with code in converters module
        ret = protocol.GAReadAlignment()
        ret.alignedQuality = read.query_qualities
        ret.alignedSequence = read.query_sequence
        ret.alignment = protocol.GALinearAlignment()
        ret.alignment.mappingQuality = read.mapping_quality
        ret.alignment.position = protocol.GAPosition()
        ret.alignment.position.referenceName = self._samFile.getrname(
            read.reference_id)
        ret.alignment.position.position = read.reference_start
        ret.alignment.cigar = []
        for operation, length in read.cigar:
            gaCigarUnit = protocol.GACigarUnit()
            gaCigarUnit.operation = SamCigar.int2ga(operation)
            gaCigarUnit.operationLength = length
            ret.alignment.cigar.append(gaCigarUnit)
        ret.duplicateFragment = SamFlags.isFlagSet(
            read.flag, SamFlags.DUPLICATE_FRAGMENT)
        ret.failedVendorQualityChecks = SamFlags.isFlagSet(
            read.flag, SamFlags.FAILED_VENDOR_QUALITY_CHECKS)
        ret.fragmentLength = read.template_length
        ret.fragmentName = read.query_name
        ret.id = None  # TODO
        ret.info = dict(read.tags)
        ret.nextMatePosition = protocol.GAPosition()
        ret.nextMatePosition.position = read.next_reference_start
        if read.next_reference_id != -1:
            ret.nextMatePosition.referenceName = \
                self._samFile.getrname(read.next_reference_id)
        ret.numberReads = SamFlags.isFlagSet(
            read.flag, SamFlags.NUMBER_READS)
        ret.properPlacement = SamFlags.isFlagSet(
            read.flag, SamFlags.PROPER_PLACEMENT)
        ret.readGroupId = self._id
        ret.readNumber = (SamFlags.isFlagSet(
            read.flag, SamFlags.READ_NUMBER_ONE) and
            SamFlags.isFlagSet(
                read.flag, SamFlags.READ_NUMBER_TWO))
        ret.secondaryAlignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SECONDARY_ALIGNMENT)
        ret.supplementaryAlignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SUPPLEMENTARY_ALIGNMENT)
        return ret
