"""
Module responsible for translating read data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import json
import os.path
import random

import pysam

import ga4gh.server.datamodel as datamodel
import ga4gh.server.datamodel.references as references
import ga4gh.server.exceptions as exceptions

import ga4gh.schemas.pb as pb
import ga4gh.schemas.protocol as protocol


def parseMalformedBamHeader(headerDict):
    """
    Parses the (probably) intended values out of the specified
    BAM header dictionary, which is incompletely parsed by pysam.
    This is caused by some tools incorrectly using spaces instead
    of tabs as a seperator.
    """
    headerString = " ".join(
        "{}:{}".format(k, v) for k, v in headerDict.items() if k != 'CL')
    ret = {}
    for item in headerString.split():
        key, value = item.split(":", 1)
        # build up dict, casting everything back to original type
        ret[key] = type(headerDict.get(key, ""))(value)
    if 'CL' in headerDict:
        ret['CL'] = headerDict['CL']
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


class AlignmentDataMixin(datamodel.PysamDatamodelMixin):
    """
    Mixin class that provides methods for getting read alignments
    from bam files
    """
    def _getReadAlignments(
            self, reference, start, end, readGroupSet, readGroup):
        """
        Returns an iterator over the specified reads
        """
        # TODO If reference is None, return against all references,
        # including unmapped reads.
        samFile = self.getFileHandle(self._dataUrl)
        referenceName = reference.getLocalId().encode()
        # TODO deal with errors from htslib
        start, end = self.sanitizeAlignmentFileFetch(start, end)
        readAlignments = samFile.fetch(referenceName, start, end)
        for readAlignment in readAlignments:
            tags = dict(readAlignment.tags)
            if readGroup is None:
                if 'RG' in tags:
                    alignmentReadGroupLocalId = tags['RG']
                    readGroupCompoundId = datamodel.ReadGroupCompoundId(
                        readGroupSet.getCompoundId(),
                        str(alignmentReadGroupLocalId))
                yield self.convertReadAlignment(
                    readAlignment, readGroupSet, str(readGroupCompoundId))
            else:
                if self._filterReads:
                    if 'RG' in tags and tags['RG'] == self._localId:
                        yield self.convertReadAlignment(
                            readAlignment, readGroupSet,
                            str(readGroup.getCompoundId()))
                else:
                    yield self.convertReadAlignment(
                        readAlignment, readGroupSet,
                        str(readGroup.getCompoundId()))

    def convertReadAlignment(self, read, readGroupSet, readGroupId):
        """
        Convert a pysam ReadAlignment to a GA4GH ReadAlignment
        """
        samFile = self.getFileHandle(self._dataUrl)
        # TODO fill out remaining fields
        # TODO refine in tandem with code in converters module
        ret = protocol.ReadAlignment()
        # ret.fragmentId = 'TODO'
        ret.aligned_quality.extend(read.query_qualities)
        ret.aligned_sequence = read.query_sequence
        if SamFlags.isFlagSet(read.flag, SamFlags.READ_UNMAPPED):
            ret.ClearField("alignment")
        else:
            ret.alignment.CopyFrom(protocol.LinearAlignment())
            ret.alignment.mapping_quality = read.mapping_quality
            ret.alignment.position.CopyFrom(protocol.Position())
            ret.alignment.position.reference_name = samFile.getrname(
                read.reference_id)
            ret.alignment.position.position = read.reference_start
            ret.alignment.position.strand = protocol.POS_STRAND
            if SamFlags.isFlagSet(read.flag, SamFlags.READ_REVERSE_STRAND):
                ret.alignment.position.strand = protocol.NEG_STRAND
            for operation, length in read.cigar:
                gaCigarUnit = ret.alignment.cigar.add()
                gaCigarUnit.operation = SamCigar.int2ga(operation)
                gaCigarUnit.operation_length = length
                gaCigarUnit.reference_sequence = ""  # TODO fix this!
        ret.duplicate_fragment = SamFlags.isFlagSet(
            read.flag, SamFlags.DUPLICATE_READ)
        ret.failed_vendor_quality_checks = SamFlags.isFlagSet(
            read.flag, SamFlags.FAILED_QUALITY_CHECK)
        ret.fragment_length = read.template_length
        ret.fragment_name = read.query_name
        for key, value in read.tags:
            # Useful for inspecting the structure of read tags
            # print("{key} {ktype}: {value}, {vtype}".format(
            #     key=key, ktype=type(key), value=value, vtype=type(value)))
            protocol.setAttribute(ret.attributes.attr[key].values, value)

        if SamFlags.isFlagSet(read.flag, SamFlags.MATE_UNMAPPED):
            ret.next_mate_position.Clear()
        else:
            ret.next_mate_position.Clear()
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
        ret.improper_placement = not SamFlags.isFlagSet(
            read.flag, SamFlags.READ_PROPER_PAIR)
        ret.read_group_id = readGroupId
        ret.secondary_alignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SECONDARY_ALIGNMENT)
        ret.supplementary_alignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SUPPLEMENTARY_ALIGNMENT)
        ret.id = readGroupSet.getReadAlignmentId(ret)
        return ret

    def openFile(self, dataFile):
        # We need to check to see if the path exists here as pysam does
        # not throw an error if the index is missing.
        if not os.path.exists(self._indexFile):
            raise exceptions.FileOpenFailedException(self._indexFile)
        try:
            return pysam.AlignmentFile(
                self._dataUrl, filepath_index=self._indexFile)
        except IOError as exception:
            # IOError thrown when the index file passed in is not actually
            # an index file... may also happen in other cases?
            raise exceptions.DataException(exception.message)


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
        self._numAlignedReads = -1
        self._numUnalignedReads = -1

    def setReferenceSet(self, referenceSet):
        """
        Sets the reference set for this ReadGroupSet to the specified value.
        """
        self._referenceSet = referenceSet

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

    def getReadGroupIds(self):
        """
        Returns the list of readGroupIds in this ReadGroupSet.
        """
        return self._readGroupIds

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
        readGroupSet.stats.CopyFrom(self.getStats())
        self.serializeAttributes(readGroupSet)
        return readGroupSet

    def getNumAlignedReads(self):
        """
        Return the number of aligned reads in this read group set
        """
        return self._numAlignedReads

    def getNumUnalignedReads(self):
        """
        Return the number of unaligned reads in this read group set
        """
        return self._numUnalignedReads

    def getPrograms(self):
        """
        Returns an array of Programs used to generate this read group set
        """
        raise NotImplementedError()

    def getReadAlignmentId(self, gaAlignment):
        """
        Returns a string ID suitable for use in the specified GA
        ReadAlignment object in this ReadGroupSet.
        """
        compoundId = datamodel.ReadAlignmentCompoundId(
            self.getCompoundId(), gaAlignment.fragment_name)
        return str(compoundId)

    def getStats(self):
        """
        Returns the GA4GH protocol representation of this read group set's
        ReadStats.
        """
        stats = protocol.ReadStats()
        stats.aligned_read_count = self._numAlignedReads
        stats.unaligned_read_count = self._numUnalignedReads
        return stats


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
        self._numAlignedReads = self._numAlignments
        self._numUnalignedReads = 0
        for i in range(numReadGroups):
            localId = "rg{}".format(i)
            readGroup = SimulatedReadGroup(
                self, localId, randomSeed + i, numAlignments)
            self.addReadGroup(readGroup)

    def getPrograms(self):
        return []

    def getReadAlignments(self, referenceId=None, start=None, end=None):
        for readGroup in self.getReadGroups():
            iterator = readGroup.getReadAlignments(referenceId, start, end)
            for alignment in iterator:
                yield alignment


class HtslibReadGroupSet(AlignmentDataMixin, AbstractReadGroupSet):
    """
    Class representing a logical collection ReadGroups.
    """
    defaultReadGroupName = "default"

    def __init__(self, parentContainer, localId):
        super(HtslibReadGroupSet, self).__init__(parentContainer, localId)
        self._programs = []
        self._dataUrl = None
        self._indexFile = None
        # Used when we populate from a file. Not defined when we populate
        # from the DB.
        self._bamHeaderReferenceSetName = None

    def getReadAlignments(self, reference, start=None, end=None):
        """
        Returns an iterator over the specified reads
        """
        return self._getReadAlignments(reference, start, end, self, None)

    def getBamHeaderReferenceSetName(self):
        """
        Returns the ReferenceSet name using in the BAM header.
        """
        return self._bamHeaderReferenceSetName

    def populateFromRow(self, readGroupSetRecord):
        """
        Populates the instance variables of this ReadGroupSet from the
        specified database row.
        """
        self._dataUrl = readGroupSetRecord.dataurl
        self._indexFile = readGroupSetRecord.indexfile
        self._programs = []
        for jsonDict in json.loads(readGroupSetRecord.programs):
            program = protocol.fromJson(json.dumps(jsonDict),
                                        protocol.Program)
            self._programs.append(program)
        stats = protocol.fromJson(readGroupSetRecord.stats, protocol.ReadStats)
        self._numAlignedReads = stats.aligned_read_count
        self._numUnalignedReads = stats.unaligned_read_count

    def populateFromFile(self, dataUrl, indexFile=None):
        """
        Populates the instance variables of this ReadGroupSet from the
        specified dataUrl and indexFile. If indexFile is not specified
        guess usual form.
        """
        self._dataUrl = dataUrl
        self._indexFile = indexFile
        if indexFile is None:
            self._indexFile = dataUrl + ".bai"
        samFile = self.getFileHandle(self._dataUrl)
        self._setHeaderFields(samFile)
        if 'RG' not in samFile.header or len(samFile.header['RG']) == 0:
            readGroup = HtslibReadGroup(self, self.defaultReadGroupName)
            self.addReadGroup(readGroup)
        else:
            for readGroupHeader in samFile.header['RG']:
                readGroup = HtslibReadGroup(self, readGroupHeader['ID'])
                readGroup.populateFromHeader(readGroupHeader)
                self.addReadGroup(readGroup)
        self._bamHeaderReferenceSetName = None
        for referenceInfo in samFile.header['SQ']:
            if 'AS' not in referenceInfo:
                infoDict = parseMalformedBamHeader(referenceInfo)
            else:
                infoDict = referenceInfo
            name = infoDict.get('AS', references.DEFAULT_REFERENCESET_NAME)
            if self._bamHeaderReferenceSetName is None:
                self._bamHeaderReferenceSetName = name
            elif self._bamHeaderReferenceSetName != name:
                raise exceptions.MultipleReferenceSetsInReadGroupSet(
                    self._dataUrl, name, self._bamFileReferenceName)
        self._numAlignedReads = samFile.mapped
        self._numUnalignedReads = samFile.unmapped

    def checkConsistency(self, dataRepository):
        pass
        # TODO verify that the references in the BAM file exist
        # in the reference set. Otherwise, we won't be able to
        # query for them.

    def _setHeaderFields(self, samFile):
        programs = []
        if 'PG' in samFile.header:
            htslibPrograms = samFile.header['PG']
            for htslibProgram in htslibPrograms:
                program = protocol.Program()
                program.id = htslibProgram['ID']
                program.command_line = htslibProgram.get(
                    'CL', pb.DEFAULT_STRING)
                program.name = htslibProgram.get('PN', pb.DEFAULT_STRING)
                program.prev_program_id = htslibProgram.get(
                    'PP', pb.DEFAULT_STRING)
                program.version = htslibProgram.get('VN', pb.DEFAULT_STRING)
                programs.append(program)
        self._programs = programs

    def getPrograms(self):
        return self._programs

    def getDataUrl(self):
        """
        Returns the data URL for this ReadGroupSet.
        """
        return self._dataUrl

    def getIndexFile(self):
        """
        Returns the index file for this ReadGroupSet.
        """
        return self._indexFile


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
        self._biosampleId = None

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
        readGroup.sample_name = pb.string(self.getSampleName())
        readGroup.biosample_id = pb.string(self.getBiosampleId())
        if referenceSet is not None:
            readGroup.reference_set_id = referenceSet.getId()
        readGroup.stats.CopyFrom(self.getStats())
        readGroup.programs.extend(self.getPrograms())
        readGroup.description = pb.string(self.getDescription())
        readGroup.experiment.CopyFrom(self.getExperiment())
        self.serializeAttributes(readGroup)
        return readGroup

    def getStats(self):
        """
        Returns the GA4GH protocol representation of this read group's
        ReadStats.
        """
        stats = protocol.ReadStats()
        stats.aligned_read_count = self.getNumAlignedReads()
        stats.unaligned_read_count = self.getNumUnalignedReads()
        # TODO base_count requires iterating through all reads
        return stats

    def getExperiment(self):
        """
        Returns the GA4GH protocol representation of this read group's
        Experiment.
        """
        experiment = protocol.Experiment()
        experiment.id = self.getExperimentId()
        experiment.instrument_model = pb.string(self.getInstrumentModel())
        experiment.sequencing_center = pb.string(self.getSequencingCenter())
        experiment.description = pb.string(self.getExperimentDescription())
        experiment.library = pb.string(self.getLibrary())
        experiment.platform_unit = pb.string(self.getPlatformUnit())
        experiment.message_create_time = self._iso8601
        experiment.message_update_time = self._iso8601
        experiment.run_time = pb.string(self.getRunTime())
        return experiment

    def getNumAlignedReads(self):
        """
        Return the number of aligned reads in the read group
        """
        return self._numAlignedReads

    def getNumUnalignedReads(self):
        """
        Return the number of unaligned reads in the read group
        """
        return self._numUnalignedReads

    def getPrograms(self):
        """
        Returns an array of Programs used to generate this read group
        """
        raise NotImplementedError()

    def getBiosampleId(self):
        return self._biosampleId

    def setBiosampleId(self, biosampleId):
        self._biosampleId = biosampleId

    def getDescription(self):
        """
        Returns a description of this read group
        """
        raise NotImplementedError()

    def getSampleName(self):
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
        self._numAlignedReads = self._parentContainer.getNumAlignedReads()
        self._numUnalignedReads = 0

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

        alignment.alignment.position.position = 0
        alignment.alignment.position.reference_name = "NotImplemented"
        alignment.alignment.position.strand = protocol.POS_STRAND
        alignment.duplicate_fragment = False
        alignment.failed_vendor_quality_checks = False

        alignment.fragment_name = "{}$simulated{}".format(
            self.getLocalId(), i)
        alignment.number_reads = 0
        alignment.improper_placement = False
        alignment.read_group_id = self.getId()
        alignment.read_number = -1
        alignment.secondary_alignment = False
        alignment.supplementary_alignment = False
        alignment.id = self._parentContainer.getReadAlignmentId(alignment)
        return alignment

    def getPrograms(self):
        return []

    def getDescription(self):
        return None

    def getSampleName(self):
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


class HtslibReadGroup(AlignmentDataMixin, AbstractReadGroup):
    """
    A readgroup based on htslib's reading of a given file
    """
    def __init__(self, parentContainer, localId):
        super(HtslibReadGroup, self).__init__(parentContainer, localId)
        # These attributes are used in AlignmentDataMixin.openFile
        self._dataUrl = parentContainer.getDataUrl()
        self._indexFile = parentContainer.getIndexFile()
        self._filterReads = localId != HtslibReadGroupSet.defaultReadGroupName
        self._biosampleId = None
        self._sampleName = None
        self._description = None
        self._predictedInsertSize = None
        self._instrumentModel = None
        self._sequencingCenter = None
        self._experimentDescription = None
        self._library = None
        self._platformUnit = None
        self._runTime = None
        self._numAlignedReads = -1  # TODO populate with metadata
        self._numUnalignedReads = -1  # TODO populate with metadata

    def populateFromHeader(self, readGroupHeader):
        """
        Populate the instance variables using the specified SAM header.
        """
        self._sampleName = readGroupHeader.get('SM', None)
        self._description = readGroupHeader.get('DS', None)
        if 'PI' in readGroupHeader:
            self._predictedInsertSize = int(readGroupHeader['PI'])
        self._instrumentModel = readGroupHeader.get('PL', None)
        self._sequencingCenter = readGroupHeader.get('CN', None)
        self._experimentDescription = readGroupHeader.get('DS', None)
        self._library = readGroupHeader.get('LB', None)
        self._platformUnit = readGroupHeader.get('PU', None)
        self._runTime = readGroupHeader.get('DT', None)

    def populateFromRow(self, readGroupRecord):
        """
        Populate the instance variables using the specified DB row.
        """
        self._sampleName = readGroupRecord.samplename
        self._biosampleId = readGroupRecord.biosampleid
        self._description = readGroupRecord.description
        self._predictedInsertSize = readGroupRecord.predictedinsertsize
        stats = protocol.fromJson(readGroupRecord.stats, protocol.ReadStats)
        self._numAlignedReads = stats.aligned_read_count
        self._numUnalignedReads = stats.unaligned_read_count
        experiment = protocol.fromJson(
            readGroupRecord.experiment, protocol.Experiment)
        self._instrumentModel = experiment.instrument_model
        self._sequencingCenter = experiment.sequencing_center
        self._experimentDescription = experiment.description
        self._library = experiment.library
        self._platformUnit = experiment.platform_unit
        self._runTime = experiment.run_time

    def getReadAlignments(self, reference, start=None, end=None):
        """
        Returns an iterator over the specified reads
        """
        return self._getReadAlignments(
            reference, start, end, self._parentContainer, self)

    def getPrograms(self):
        return self._parentContainer.getPrograms()

    def getDescription(self):
        return self._description

    def getSampleName(self):
        return self._sampleName

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
