"""
Provides classes that take protocol requests, send that request to
the server, and write a particular genomics file type with the results.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import re

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
        super(VcfConverter, self).__init__(httpClient)
        self.vcf = None
        self.variantSet = None
        self.outputStream = outputStream
        self.searchVariantsRequest = searchVariantsRequest
        self.searchVariantSetsRequest = searchVariantSetsRequest

    def convert(self):
        self.getVariantSet()
        self.getCallSets()
        self.getVariants()
        self.createVcf()

    def getVariantSet(self):
        response = self._httpClient.searchVariantSets(
            self.searchVariantSetsRequest)
        variantSets = [variantSet for variantSet in response]
        if len(variantSets) == 0:
            raise VcfException("VariantSet not found")
        if len(variantSets) > 1:
            raise VcfException("More than one VariantSet returned")
        self.variantSet = variantSets[0]
        return variantSet

    def getCallSets(self):
        request = protocol.GASearchCallSetsRequest()
        request.variantSetIds = [self.variantSet.id]
        response = self._httpClient.searchCallSets(request)
        callSets = [callSet for callSet in response]
        self.callSets = callSets
        return callSets

    def getVariants(self):
        response = self._httpClient.searchVariants(
            self.searchVariantsRequest)
        self.variantGenerator = response
        return response

    def createVcf(self):
        self.vcf = Vcf(self.outputStream, self.variantSet, self.callSets)
        self.vcf.writeHeader()
        self.vcf.writeData(self.variantGenerator)


class VcfMetadataInfoLine(object):
    """
    A line in the metadata of a VCF file
    """
    # see pysam source file cvcf.pyx
    numberTypes = {
        'NT_UNKNOWN': 0,
        'NT_NUMBER': 1,
        'NT_ALLELES': 2,
        'NT_NR_ALLELES': 3,
        'NT_GENOTYPES': 4,
        'NT_PHASED_GENOTYPES': 5,
    }

    def __init__(self, infoDict):
        self.infoDict = infoDict
        for key, value in self.infoDict.items():
            setattr(self, key, value)

    def getNumbertype(self):
        # The google api seems to return one of the following for number:
        # - an integer
        # - a period
        # - None (only on custom header lines)
        if re.match('[0-9]+', self.number) is not None:
            return self.numberTypes['NT_NUMBER']
        elif self.number == 'R':
            return self.numberTypes['NT_ALLELES']
        elif self.number == 'A':
            return self.numberTypes['NT_NR_ALLELES']
        elif self.number == 'G':
            return self.numberTypes['NT_GENOTYPES']
        else:
            return self.numberTypes['NT_UNKNOWN']

    def __getattr__(self, name):
        # the pysam implementation looks up the following attributes
        # on this object
        if name == 'numbertype':
            return self.getNumbertype()
        elif name == "missingvalue":
            return None
        raise VcfException(
            "VcfMetadataInfoLine has no attribute {0}".format(name))


class VcfMetadata(object):
    """
    The metadata of a VCF file
    """
    def __init__(self):
        self.info = {}
        self.filter_ = {}
        self.format_ = {}
        self.header = []
        self.switch = {
            'INFO': self.info,
            'FILTER': self.filter_,
            'FORMAT': self.format_,
        }

    def processVariantSet(self, variantSet):
        for metadata in variantSet.metadata:
            infoLine = self._processMetadata(metadata)
            if infoLine.key in self.switch:
                self.switch[infoLine.key][infoLine.id] = infoLine
            elif infoLine.key == 'fileformat':
                # fileformat already written by pysam
                continue
            else:
                keyVal = self._processHeader(infoLine)
                self.header.append(keyVal)

    def _processHeader(self, infoLine):
        if infoLine.key == 'contig':
            valueStr = []
            for ilKey, ilVal in infoLine.info.items():
                ilVal = ilVal[0]
                if ilVal.find(' ') != -1:
                    ilVal = '"{0}"'.format(ilVal)
                fragment = '{0}={1}'.format(ilKey, ilVal)
                valueStr.append(fragment)
            joinedValueStr = ','.join(valueStr)
            value = "<{0}>".format(joinedValueStr)
        elif infoLine.key == 'ALT':
            value = "<ID={id}, Description={description}>".format(
                **infoLine.infoDict)
        # TODO special cases for other fields: assembly, etc...
        # see VCF doc section 1.2
        else:
            value = infoLine.value
        return (infoLine.key, value)

    def _processMetadata(self, metadata):
        infoDict = metadata.toJsonDict()
        infoLine = VcfMetadataInfoLine(infoDict)
        return infoLine


class VcfRow(object):
    """
    Represents a row in a VCF file
    """
    defaultChrom = 'defaultChrom'
    defaultQuality = '50'
    defaultFilter = []
    defaultFormat = []
    defaultNoSample = []

    def __init__(self, variant, samples):
        self.variantDict = {
            'chrom': self.defaultChrom,
            # pysam increments pos, so it needs to be an int
            'pos': int(variant.start),
            'id': variant.id,
            'ref': variant.referenceBases,
            'alt': variant.alternateBases,
            'qual': self.defaultQuality,
            'filter': self.defaultFilter,
            'info': variant.info,
            'format': self.defaultFormat,
        }
        for call in variant.calls:
            self.variantDict[call.callSetName] = call.info
        for sample in samples:
            if sample not in self.variantDict:
                self.variantDict[sample] = self.defaultNoSample

    def toDict(self):
        return self.variantDict


class Vcf(object):
    """
    A convenience wrapper around pysam.VCF functionality
    providing methods to write a VCF file.
    """
    def __init__(self, stream, variantSet, callSets):
        self.stream = stream
        self.vcf = pysam.VCF()
        self.variantSet = variantSet
        self.callSets = callSets
        self.data = []
        self.headerWritten = False

    def writeHeader(self):
        self.vcfMetadata = VcfMetadata()
        self.vcfMetadata.processVariantSet(self.variantSet)
        self.vcf.setinfo(self.vcfMetadata.info)
        self.vcf.setfilter(self.vcfMetadata.filter_)
        self.vcf.setformat(self.vcfMetadata.format_)
        self.vcf.setheader(self.vcfMetadata.header)
        self.samples = [callSet.sampleId for callSet in self.callSets]
        self.vcf.setsamples(self.samples)
        self.vcf.writeheader(self.stream)
        self.headerWritten = True

    def writeData(self, variantGenerator):
        if self.headerWritten is False:
            raise VcfException("Need to write header before writing data")
        requestLimit = 1  # TODO find beter way to limit data
        requestNumber = 0
        for variant in variantGenerator:
            if requestNumber >= requestLimit:
                break
            row = VcfRow(variant, self.samples)
            # for debugging pysam, uncomment these lines
            # and comment out the write_data call
            # fakePysam = FakePysam(
            #    self.samples, self.vcfMetadata.format_, self.vcfMetadata.info)
            # fakePysam.write_data(self.stream, row.toDict())
            self.vcf.write_data(self.stream, row.toDict())
            requestNumber += 1


#############################
# TODO take out below
# this is copy-and-pasted from pysam's implementation to help with debugging
#############################


class FakePysam(object):

    NT_UNKNOWN = 0
    NT_NUMBER = 1
    NT_ALLELES = 2
    NT_NR_ALLELES = 3
    NT_GENOTYPES = 4
    NT_PHASED_GENOTYPES = 5

    def __init__(self, samples, format_, info):
        self._version = 40
        self._samples = samples
        self._format = format_
        self._info = info

    def write_data(self, stream, data):
        required = ['chrom','pos','id','ref','alt','qual','filter','info','format'] + self._samples
        for k in required:
            if k not in data: raise ValueError("Required key %s not found in data" % str(k))
        if data['alt'] == []: alt = "."
        else: alt = ",".join(data['alt'])
        if data['filter'] == None: filter = "."
        elif data['filter'] == []: 
            if self._version == 33: filter = "0"
            else: filter = "PASS"
        else: filter = ';'.join(data['filter'])
        if data['qual'] == -1: qual = "."
        else: qual = str(data['qual'])

        output = [data['chrom'], 
                  str(data['pos']+1),   # change to 1-based position
                  data['id'],
                  data['ref'],
                  alt,
                  qual,
                  filter,
                  self.format_formatdata(
                      data['info'], self._info, separator=";"),
                  self.format_formatdata(
                      data['format'], self._format, value=False)]
        
        for s in self._samples:
            output.append(self.format_formatdata(
                data[s], self._format, key=False))
        
        stream.write( "\t".join(output) + "\n" )

    def format_formatdata( self, data, format, key=True, value=True, separator=":" ):
        output, sdata = [], []
        if type(data) == type([]): # for FORMAT field, make data with dummy values
            d = {}
            for k in data: d[k] = []
            data = d
        # convert missing values; and silently add definitions if required
        for k in data:
            self._add_definition( format, k, data[k], "(output)" )
            for idx,v in enumerate(data[k]):
                if v == format[k].missingvalue: data[k][idx] = "."
        # make sure GT comes first; and ensure fixed ordering; also convert GT data back to string
        for k in data: 
            if k != 'GT': sdata.append( (k,data[k]) )
        sdata.sort()
        if 'GT' in data:
            sdata = [('GT',map(self.convertGTback,data['GT']))] + sdata
        for k,v in sdata:
            if v == []: v = None
            if key and value:
                if v != None: output.append( k+"="+','.join(map(str,v)) )
                else: output.append( k )
            elif key: output.append(k)
            elif value:
                if v != None: output.append( ','.join(map(str,v)) )
                else: output.append( "." )                    # should not happen
        # snip off trailing missing data
        while len(output) > 1:
            last = output[-1].replace(',','').replace('.','')
            if len(last)>0: break
            output = output[:-1]
        return separator.join(output)

    def _add_definition(self, formatdict, key, data, line ):
        if key in formatdict: return
        self.error(line,self.ERROR_UNKNOWN_KEY,key)
        if data == None:
            formatdict[key] = FORMAT(key,self.NT_NUMBER,0,"Flag","(Undefined tag)",".")
            return
        if data == []: data = [""]             # unsure what type -- say string
        if type(data[0]) == type(0.0):
            formatdict[key] = FORMAT(key,self.NT_UNKNOWN,-1,"Float","(Undefined tag)",None)
            return
        if type(data[0]) == type(0):
            formatdict[key] = FORMAT(key,self.NT_UNKNOWN,-1,"Integer","(Undefined tag)",None)
            return
        formatdict[key] = FORMAT(key,self.NT_UNKNOWN,-1,"String","(Undefined tag)",".")
