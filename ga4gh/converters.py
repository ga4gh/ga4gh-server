"""
Provides classes that take protocol requests, send that request to
the server, and write a particular genomics file type with the results.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections

import pysam

import ga4gh.protocol as protocol


class AbstractConverter(object):
    """
    Abstract base class for converter classes
    """
    def __init__(self, httpClient):
        self._httpClient = httpClient


##############################################################################
# SAM
##############################################################################


class SamException(Exception):
    """
    Something that went wrong during converting a SAM file
    """


class SamConverter(AbstractConverter):
    """
    Converts a request to a SAM file
    """
    def __init__(self, httpClient, searchReadsRequest, outputFile,
                 binaryOutput):
        super(SamConverter, self).__init__(httpClient)
        self._searchReadsRequest = searchReadsRequest
        self._outputFile = outputFile
        self._binaryOutput = binaryOutput

    def convert(self):
        header = self._getHeader()
        targetIds = self._getTargetIds(header)
        # pysam can't write to file streams (except for stdout)
        # http://pysam.readthedocs.org/en/latest/usage.html#using-streams
        if self._binaryOutput:
            flags = "wb"
        else:
            flags = "wh"  # h for header
        fileString = "-"
        if self._outputFile is not None:
            fileString = self._outputFile
        alignmentFile = pysam.AlignmentFile(
            fileString, flags, header=header)
        iterator = self._httpClient.searchReads(self._searchReadsRequest)
        for read in iterator:
            alignedSegment = SamLine.toAlignedSegment(read, targetIds)
            alignmentFile.write(alignedSegment)
        alignmentFile.close()

    def _getHeader(self):
        # TODO where to get actual values for header?
        # need some kind of getReadGroup(readGroupId) method in protocol
        # just add these dummy lines for now
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [
                {'LN': 1575, 'SN': 'chr1'},
                {'LN': 1584, 'SN': 'chr2'},
            ],
        }
        return header

    def _getTargetIds(self, header):
        # this seems to be how pysam sets the target ids
        targetIds = collections.defaultdict(int)
        targetId = 0
        if 'SQ' in header:
            headerLines = header['SQ']
            for headerLine in headerLines:
                refName = headerLine['SN']
                targetIds[refName] = targetId
                targetId += 1
        return targetIds


class SamLine(object):
    """
    Methods for processing a line in a SAM file
    """
    _encoding = 'utf8'

    # see http://pysam.readthedocs.org/en/latest/api.html
    # #pysam.AlignedSegment.cigartuples
    _cigarMap = {
        protocol.GACigarOperation.ALIGNMENT_MATCH: 0,
        protocol.GACigarOperation.INSERT: 1,
        protocol.GACigarOperation.DELETE: 2,
        protocol.GACigarOperation.SKIP: 3,
        protocol.GACigarOperation.CLIP_SOFT: 4,
        protocol.GACigarOperation.CLIP_HARD: 5,
        protocol.GACigarOperation.PAD: 6,
        protocol.GACigarOperation.SEQUENCE_MATCH: 7,
        protocol.GACigarOperation.SEQUENCE_MISMATCH: 8,
    }

    # see tables in SAM spec, section 1.5
    _tagReservedFieldPrefixes = set(["X", "Y", "Z", ])
    _tagIntegerFields = set([
        "AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", "HI", "IH", "MQ",
        "NH", "NM", "OP", "PQ", "SM", "TC", "UQ", ])
    _tagStringFields = set([
        "BC", "BQ", "CC", "CO", "CQ", "CS", "CT", "E2", "FS", "LB", "MC",
        "MD", "OQ", "OC", "PG", "PT", "PU", "QT", "Q2", "R2", "RG", "RT",
        "SA", "U2", ])
    _tagIntegerArrayFields = set(["FZ", ])

    def __init__(self):
        raise SamException("SamLine can't be instantiated")

    @classmethod
    def toAlignedSegment(cls, read, targetIds):
        ret = pysam.AlignedSegment()
        # QNAME
        ret.query_name = read.fragmentName.encode(cls._encoding)
        # SEQ
        ret.query_sequence = read.alignedSequence.encode(cls._encoding)
        # FLAG
        ret.flag = cls.toSamFlag(read)
        # RNAME
        refName = read.alignment.position.referenceName
        ret.reference_id = targetIds[refName]
        # POS
        ret.reference_start = int(read.alignment.position.position)
        # MAPQ
        ret.mapping_quality = read.alignment.mappingQuality
        # CIGAR
        ret.cigar = cls.toCigar(read)
        # RNEXT
        nextRefName = read.nextMatePosition.referenceName
        ret.next_reference_id = targetIds[nextRefName]
        # PNEXT
        ret.next_reference_start = int(read.nextMatePosition.position)
        # TLEN
        ret.template_length = read.fragmentLength
        # QUAL
        ret.query_qualities = read.alignedQuality
        ret.tags = cls.toTags(read)
        return ret

    @classmethod
    def toSamFlag(cls, read):
        flag = 0
        if read.numberReads:
            flag += 0x1
        if read.properPlacement:
            flag += 0x2
        if read.readNumber:
            flag += 0x40
            flag += 0x80
        if read.secondaryAlignment:
            flag += 0x100
        if read.failedVendorQualityChecks:
            flag += 0x200
        if read.duplicateFragment:
            flag += 0x400
        if read.supplementaryAlignment:
            flag += 0x800
        return flag

    @classmethod
    def toCigar(cls, read):
        cigarTuples = []
        for gaCigarUnit in read.alignment.cigar:
            operation = cls._cigarMap[gaCigarUnit.operation]
            length = int(gaCigarUnit.operationLength)
            cigarTuple = (operation, length)
            cigarTuples.append(cigarTuple)
        return tuple(cigarTuples)

    @classmethod
    def _parseTagValue(cls, tag, value):
        if tag[0] in cls._tagReservedFieldPrefixes:
            # user reserved fields... not really sure what to do here
            return value[0].encode(cls._encoding)
        elif tag in cls._tagIntegerFields:
            return int(value[0])
        elif tag in cls._tagStringFields:
            return value[0].encode(cls._encoding)
        elif tag in cls._tagIntegerArrayFields:
            return [int(integerString) for integerString in value]
        else:
            raise SamException("unrecognized tag '{}'".format(tag))

    @classmethod
    def toTags(cls, read):
        tags = []
        for tag, value in read.info.items():
            val = cls._parseTagValue(tag, value)
            tagTuple = (tag, val)
            tags.append(tagTuple)
        retval = tuple(tags)
        return retval


##############################################################################
# VCF
##############################################################################


class VcfException(Exception):
    pass


class VcfConverter(AbstractConverter):
    """
    Converts a request to a VCF file
    """
    def __init__(self, httpClient, outputStream, searchVariantSetsRequest,
                 searchVariantsRequest):
        raise NotImplementedError
