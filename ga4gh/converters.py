"""
Provides classes that take protocol requests, send that request to
the server, and write a particular genomics file type with the results.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections

import pysam

import ga4gh.datamodel.reads as reads
import ga4gh.protocol as protocol


class AbstractConverter(object):
    """
    Abstract base class for converter classes
    """
    def __init__(
            self, container, objectIterator, outputFile, binaryOutput):
        self._container = container
        self._objectIterator = objectIterator
        self._outputFile = outputFile
        self._binaryOutput = binaryOutput


##############################################################################
# SAM
##############################################################################


class SamException(Exception):
    """
    Something that went wrong during converting a SAM file
    """


class SamConverter(object):
    """
    Converts a requested range from a GA4GH server into a SAM file.
    """
    def __init__(
            self, client, readGroupId=None, referenceId=None,
            start=None, end=None, outputFileName=None, binaryOutput=False):
        self._client = client
        self._readGroup = self._client.getReadGroup(readGroupId)
        self._reference = self._client.getReference(referenceId)
        self._start = start
        self._end = end
        self._outputFileName = outputFileName
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
        if self._outputFileName is not None:
            fileString = self._outputFileName
        alignmentFile = pysam.AlignmentFile(fileString, flags, header=header)
        iterator = self._client.searchReads(
            [self._readGroup.id], self._reference.id, self._start, self._end)
        for read in iterator:
            alignedSegment = SamLine.toAlignedSegment(read, targetIds)
            alignmentFile.write(alignedSegment)
        alignmentFile.close()

    def _getHeader(self):
        # Create header information using self._reference
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{
                'LN': self._reference.length,
                'SN': self._reference.name
            }]
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
        if read.nextMatePosition is None:
            ret.next_reference_id = -1
        else:
            nextRefName = read.nextMatePosition.referenceName
            ret.next_reference_id = targetIds[nextRefName]
        # PNEXT
        if read.nextMatePosition is None:
            ret.next_reference_start = -1
        else:
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
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.NUMBER_READS)
        if read.properPlacement:
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.PROPER_PLACEMENT)
        if read.alignment.position.strand == protocol.Strand.NEG_STRAND:
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.REVERSED)
        if (read.nextMatePosition is not None and
                read.nextMatePosition.strand == protocol.Strand.NEG_STRAND):
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.NEXT_MATE_REVERSED)
        if read.readNumber == 0:
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.READ_NUMBER_ONE)
        elif read.readNumber == 1:
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.READ_NUMBER_TWO)
        elif read.readNumber == 2:
            flag = reads.SamFlags.setFlag(
                flag,
                reads.SamFlags.READ_NUMBER_ONE |
                reads.SamFlags.READ_NUMBER_TWO)
        if read.secondaryAlignment:
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.SECONDARY_ALIGNMENT)
        if read.failedVendorQualityChecks:
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.FAILED_VENDOR_QUALITY_CHECKS)
        if read.duplicateFragment:
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.DUPLICATE_FRAGMENT)
        if read.supplementaryAlignment:
            flag = reads.SamFlags.setFlag(
                flag, reads.SamFlags.SUPPLEMENTARY_ALIGNMENT)
        return flag

    @classmethod
    def toCigar(cls, read):
        cigarTuples = []
        for gaCigarUnit in read.alignment.cigar:
            operation = reads.SamCigar.ga2int(gaCigarUnit.operation)
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
    Converts the Variants represented by a SearchVariantsRequest into
    VCF format using pysam.
    """
    def _writeHeader(self):
        variantSet = self._container
        # TODO convert this into pysam types and write to the output file.
        # For now, just print out some stuff to demonstrate how to get the
        # attributes we have.
        print("ID = ", variantSet.id)
        print("Dataset ID = ", variantSet.datasetId)
        print("Metadata = ")
        for metadata in variantSet.metadata:
            print("\t", metadata)

    def _writeBody(self):
        for variant in self._objectIterator:
            # TODO convert each variant object into pysam objects and write to
            # the output file. For now, just print the first variant and break.
            print(variant)
            break

    def convert(self):
        """
        Run the conversion process.
        """
        # TODO allocate the pysam VCF object which can be used for the
        # conversion process. See the convert method for ga2sam above.
        self._writeHeader()
        self._writeBody()
