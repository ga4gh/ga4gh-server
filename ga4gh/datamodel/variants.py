"""
Module responsible for translating variant data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import glob
import os
import random

import pysam
import wormtable as wt

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


def convertVCFPhaseset(vcfPhaseset):
    """
    Parses the VCF phaseset string
    """
    if vcfPhaseset is not None and vcfPhaseset != ".":
        phaseset = vcfPhaseset
    else:
        phaseset = "*"
    return phaseset


def convertVCFGenotype(vcfGenotype, vcfPhaseset):
    """
    Parses the VCF genotype and VCF phaseset strings
    """
    phaseset = None
    if vcfGenotype is not None:
        delim = "/"
        if "|" in vcfGenotype:
            delim = "|"
            phaseset = convertVCFPhaseset(vcfPhaseset)
        if "." in vcfGenotype:
            genotype = [-1]
        else:
            genotype = map(int, vcfGenotype.split(delim))
    else:
        genotype = [-1]

    return genotype, phaseset


def variantSetFactory(variantSetId, relativePath):
    if variantSetId.endswith(".wt"):
        return WormtableVariantSet(variantSetId, relativePath)
    else:
        return HtslibVariantSet(variantSetId, relativePath)


class CallSet(object):
    """
    Class representing a CallSet. A CallSet basically represents the
    metadata associated with a single VCF sample column.
    """
    def __init__(self, variantSet, callSetId, sampleName):
        self._variantSet = variantSet
        self._id = callSetId
        self._sampleName = sampleName

    def getId(self):
        # TODO this should be in the superclass, DatamodelObject.
        return self._id

    def toProtocolElement(self):
        """
        Returns the representation of this CallSet as the corresponding
        ProtocolElement.
        """
        gaCallSet = protocol.GACallSet()
        gaCallSet.created = self._variantSet.getCreationTime()
        gaCallSet.updated = self._variantSet.getUpdatedTime()
        gaCallSet.id = self._id
        gaCallSet.name = self._sampleName
        gaCallSet.sampleId = self._sampleName
        return gaCallSet


class AbstractVariantSet(object):
    """
    An abstract base class of a variant set
    """
    def __init__(self, id_):
        self._id = id_
        self._callSetIdMap = {}
        self._callSetIds = []
        self._creationTime = None
        self._updatedTime = None

    def getCreationTime(self):
        """
        Returns the creation time for this variant set.
        """
        return self._creationTime

    def getUpdatedTime(self):
        """
        Returns the time this variant set was last updated.
        """
        return self._updatedTime

    def getCallSetId(self, sampleName):
        """
        Returns the callSetId for the specified sampleName in this
        VariantSet.
        """
        return "{0}.{1}".format(self.getId(), sampleName)

    def addCallSet(self, sampleName):
        """
        Adds a CallSet for the specified sample name.
        """
        callSetId = self.getCallSetId(sampleName)
        callSet = CallSet(self, callSetId, sampleName)
        self._callSetIdMap[callSetId] = callSet
        self._callSetIds.append(callSetId)

    def getCallSetIdMap(self):
        """
        Returns the map of callSetIds to CallSet objects in this
        VariantSet.
        """
        return self._callSetIdMap

    def getCallSetIds(self):
        """
        Returns the list of callSetIds in this VariantSet.
        """
        return self._callSetIds

    def getCallSets(self):
        """
        Returns an iterator over the CallSets for this VariantSet.
        """
        return self._callSetIdMap.values()

    def toProtocolElement(self):
        """
        Converts this VariantSet into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.GAVariantSet()
        protocolElement.id = self._id
        protocolElement.datasetId = "NotImplemented"
        protocolElement.metadata = self.getMetadata()
        return protocolElement

    def getId(self):
        """
        Returns the ID of this VariantSet.

        TODO: this should be pushed into a superclass, and use an
        instance variant self._id.
        """
        return self._id

    def getNumVariants(self):
        """
        Returns the number of variants contained in this VariantSet.
        """
        raise NotImplementedError()

    def _createGaVariant(self):
        """
        Convenience method to set the common fields in a GA Variant
        object from this variant set.
        """
        ret = protocol.GAVariant()
        ret.created = self._creationTime
        ret.updated = self._updatedTime
        ret.variantSetId = self.getId()
        return ret


class SimulatedVariantSet(AbstractVariantSet):
    """
    A variant set that doesn't derive from a data store.
    Used mostly for testing.
    """
    def __init__(self, randomSeed, numCalls, variantDensity, variantSetId):
        super(SimulatedVariantSet, self).__init__(variantSetId)
        self._randomSeed = randomSeed
        self._numCalls = numCalls
        for j in range(numCalls):
            self.addCallSet("simCallSet_{}".format(j))
        self._variantDensity = variantDensity
        now = protocol.convertDatetime(datetime.datetime.now())
        self._creationTime = now
        self._updatedTime = now

    def getNumVariants(self):
        return 0

    def getMetadata(self):
        ret = []
        return ret

    def getVariants(self, referenceName, startPosition, endPosition,
                    variantName, callSetIds):
        randomNumberGenerator = random.Random()
        i = startPosition
        while i < endPosition:
            randomNumberGenerator.seed(self._randomSeed + i)
            if randomNumberGenerator.random() < self._variantDensity:
                variant = self.generateVariant(
                    self._id, referenceName, i, randomNumberGenerator)
                yield variant
            i += 1

    def generateVariant(self, variantSetId, referenceName, position,
                        randomNumberGenerator):
        """
        Generate a random variant for the specified position using the
        specified random number generator. This generator should be seeded
        with a value that is unique to this position so that the same variant
        will always be produced regardless of the order it is generated in.
        """
        variant = self._createGaVariant()
        variant.names = []
        variant.referenceName = referenceName
        variant.id = "{0}:{1}:{2}".format(
            variant.variantSetId, referenceName, position)
        variant.start = position
        variant.end = position + 1  # SNPs only for now
        bases = ["A", "C", "G", "T"]
        ref = randomNumberGenerator.choice(bases)
        variant.referenceBases = ref
        alt = randomNumberGenerator.choice(
            [base for base in bases if base != ref])
        variant.alternateBases = [alt]
        variant.calls = []
        for callSet in self.getCallSets():
            call = protocol.GACall()
            call.callSetId = callSet.getId()
            # for now, the genotype is either [0,1], [1,1] or [1,0] with equal
            # probability; probably will want to do something more
            # sophisticated later.
            randomChoice = randomNumberGenerator.choice(
                [[0, 1], [1, 0], [1, 1]])
            call.genotype = randomChoice
            # TODO What is a reasonable model for generating these likelihoods?
            # Are these log-scaled? Spec does not say.
            call.genotypeLikelihood = [-100, -100, -100]
            variant.calls.append(call)
        return variant


def _encodeValue(value):
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    else:
        return [str(value)]


_nothing = object()


def isEmptyIter(it):
    """Return True iff the iterator is empty or exhausted"""
    return next(it, _nothing) is _nothing


class HtslibVariantSet(AbstractVariantSet):
    """
    Class representing a single variant set backed by a directory of indexed
    VCF or BCF files.
    """
    def __init__(self, variantSetId, vcfPath):
        super(HtslibVariantSet, self).__init__(variantSetId)
        # ctime is in seconds, and we want milliseconds since the epoch
        ctimeInMillis = int(os.path.getctime(vcfPath) * 1000)
        self._creationTime = ctimeInMillis
        self._updatedTime = ctimeInMillis
        self._chromFileMap = {}
        self._metadata = None
        for pattern in ['*.bcf', '*.vcf.gz']:
            for filename in glob.glob(os.path.join(vcfPath, pattern)):
                self._addFile(filename)

    def _updateMetadata(self, variantFile):
        """
        Updates the metadata for his variant set based on the specified
        variant file, and ensures that it is consistent with already
        existing metadata.
        """
        expMsg = "Metadata of {} is not consistent".format(
            variantFile.filename)
        metadata = self._getMetadataFromVcf(variantFile)
        if self._metadata is None:
            self._metadata = metadata
        else:
            if self._metadata != metadata:
                # TODO CHANGE EXCEPTION
                raise Exception(expMsg)

    def getNumVariants(self):
        """
        Returns the total number of variants in this VariantSet.
        """
        # TODO How do we get the number of records in a VariantFile?
        return 0

    def _updateCallSetIds(self, variantFile):
        """
        Updates the call set IDs based on the specified variant file.
        """
        # If this is the first file, we add in the samples. If not, we check
        # for consistency.
        if len(self._callSetIdMap) == 0:
            for sample in variantFile.header.samples:
                self.addCallSet(sample)
        else:
            callSetIds = set([
                self.getCallSetId(sample)
                for sample in variantFile.header.samples])
            if callSetIds != set(self._callSetIdMap.keys()):
                # TODO CHANGE EXCEPTION
                raise Exception("Inconsistent sample names in VCF")

    def _addFile(self, filename):
        varFile = pysam.VariantFile(filename)
        if varFile.index is None:
            # TODO CHANGE EXCEPTION
            raise Exception("VCF/BCF files must be indexed")
        for chrom in varFile.index:
            # Unlike Tabix indices, CSI indices include all contigs defined
            # in the BCF header.  Thus we must test each one to see if
            # records exist or else they are likely to trigger spurious
            # overlapping errors.
            if not isEmptyIter(varFile.fetch(chrom)):
                if chrom in self._chromFileMap:
                    # TODO CHANGE EXCEPTION
                    raise Exception("cannot have overlapping VCF/BCF files.")
                self._updateMetadata(varFile)
                self._updateCallSetIds(varFile)
                self._chromFileMap[chrom] = varFile

    def _convertGaCall(self, recordId, name, pysamCall, genotypeData):
        call = protocol.GACall()
        call.callSetId = recordId
        call.callSetName = name

        # TODO:
        # NOTE: THE FOLLOWING TWO LINES IS NOT THE INTENED IMPLEMENTATION,
        ###########################################
        call.phaseset = None
        call.genotype, call.phaseset = convertVCFGenotype(
            genotypeData, call.phaseset)
        ###########################################

        # THEY SHOULD BE REPLACED BY THE FOLLOWING, ONCE NEW PYSAM
        # RELEASE SUPPORTS phaseset. AS WELL AS REMOVING genotypeData
        # FROM THE FUNCTION CALL

        ###########################################
        # call.genotype = list(pysamCall.allele_indices)
        # call.phaseset = pysamCall.phaseset
        ###########################################

        call.genotypeLikelihood = []
        for key, value in pysamCall.iteritems():
            if key == 'GL' and value is not None:
                call.genotypeLikelihood = list(value)
            elif key != 'GT':
                call.info[key] = _encodeValue(value)
        return call

    def convertVariant(self, record, callSetIds):
        """
        Converts the specified pysam variant record into a GA4GH Variant
        object. Only calls for the specified list of callSetIds will
        be included.
        """
        variant = self._createGaVariant()
        # N.B. record.pos is 1-based
        #      also consider using record.start-record.stop
        variant.id = "{0}:{1}:{2}".format(self._id,
                                          record.contig,
                                          record.pos)
        variant.referenceName = record.contig
        if record.id is not None:
            variant.names = record.id.split(';')
        variant.start = record.start          # 0-based inclusive
        variant.end = record.stop             # 0-based exclusive
        variant.referenceBases = record.ref
        if record.alts is not None:
            variant.alternateBases = list(record.alts)
        # record.filter and record.qual are also available, when supported
        # by GAVariant.
        for key, value in record.info.iteritems():
            if value is not None:
                variant.info[key] = _encodeValue(value)

        # NOTE: THE LABELED LINES SHOULD BE REMOVED ONCE PYSAM SUPPORTS
        # phaseset

        sampleData = record.__str__().split()[9:]  # REMOVAL
        variant.calls = []
        sampleIterator = 0  # REMOVAL
        for name, call in record.samples.iteritems():
            if name in callSetIds:
                genotypeData = sampleData[sampleIterator].split(
                    ":")[0]  # REMOVAL
                variant.calls.append(self._convertGaCall(
                    record.id, name, call, genotypeData))  # REPLACE
            sampleIterator += 1  # REMOVAL
        return variant

    def getVariants(self, referenceName, startPosition, endPosition,
                    variantName, callSetIds):
        """
        Returns an iterator over the specified variants. The parameters
        correspond to the attributes of a GASearchVariantsRequest object.
        """
        if variantName is not None:
            raise exceptions.NotImplementedException(
                "Searching by variantName is not supported")
        if callSetIds is None:
            callSetIds = self._callSetIdMap.keys()

        if referenceName in self._chromFileMap:
            varFile = self._chromFileMap[referenceName]
            cursor = varFile.fetch(
                referenceName, startPosition, endPosition)
            for record in cursor:
                yield self.convertVariant(record, callSetIds)

    def getMetadata(self):
        return self._metadata

    def _getMetadataFromVcf(self, varFile):
        # All the metadata is available via each varFile.header, including:
        #    records: header records
        #    version: VCF version
        #    samples -- not immediately needed
        #    contigs -- not immediately needed
        #    filters -- not immediately needed
        #    info
        #    formats

        def buildMetadata(
                key, type="String", number="1", value="", id="",
                description=""):  # All input are strings
            metadata = protocol.GAVariantSetMetadata()
            metadata.key = key
            metadata.value = value
            metadata.id = id
            metadata.type = type
            metadata.number = number
            metadata.description = description
            return metadata

        ret = []
        header = varFile.header
        ret.append(buildMetadata(key="version", value=header.version))
        for formatKey, value in header.formats.items():
            if formatKey != "GT":
                ret.append(buildMetadata(
                    key="FORMAT.{}".format(value.name), type=value.type,
                    number="{}".format(value.number)))
                # NOTE: description is not currently implemented as a member
                # of VariantMetadata in pysam/cbcf.pyx
        for infoKey, value in header.info.items():
            ret.append(buildMetadata(
                key="INFO.{}".format(value.name), type=value.type,
                number="{}".format(value.number)))
            # NOTE: description is not currently implemented as a member
            # of VariantMetadata in pysam/cbcf.pyx
        return ret

#################################################################
#
# Everything below here should be part of the Wormtable
# based implementation and can be deleted.
#
#################################################################


class WormtableVariantSet(object):
    """
    Class representing a single variant set backed by a wormtable directory.
    We assume that VCF data has been converted to wormtable format using
    vcf2wt, which maps VCF data to wormtable columns in the following
    way:
    o The zeroth column (row_id) in a row is the primary key;
    o The next 7 columns correspond to the fixed fields in VCF, and are
      called CHROM, POS, etc.
    o There is then a single column for each INFO field defined in the
      VCF header. The column names are prefixed with "INFO.", and so
      we might have INFO.AF, INFO.THETA etc.
    o There is then a set of columns for each sample, one for each field
      defined in the FORMAT rows of the header. These are prefixed by the
      sample ID and suffixed with the field name, and so we have fields
      like HG00096.GT, NA20828.GL etc.
    The ``wtadmin show`` command is useful to see the columns in a table.

    We assume that indexes on the combination of the CHROM and POS columns
    and the CHROM and ID columns have been built. This is to support
    efficient retrieval of rows based on both coordinates and names
    within a chromosome.
    """
    # Positions of the fixed columns within a row
    CHROM_COL = 1
    POS_COL = 2
    ID_COL = 3
    REF_COL = 4
    ALT_COL = 5
    QUAL_COL = 6
    FILTER_COL = 7
    # These must be bytes literals for integration with wormtable.
    GENOTYPE_LIKELIHOOD_NAME = b"GL"
    GENOTYPE_NAME = b"GT"
    PHASESET_NAME = b"PS"

    def __init__(self, variantSetId, wtDir):
        """
        Allocates a new WormtableVariantSet with the specified variantSetId
        based on the specified wormtable directory.
        """
        self._sampleNames = []
        self._variantSetId = variantSetId
        self._wtDir = wtDir
        self._table = wt.open_table(wtDir)
        self._chromPosIndex = self._table.open_index("CHROM+POS")
        self._chromIdIndex = self._table.open_index("CHROM+ID")
        self._sampleCols = {}
        self._infoCols = []
        self._firstSamplePosition = -1
        ctimeInMillis = int(os.path.getctime(wtDir) * 1000)
        # ctime is in seconds, and we want milliseconds since the epoch
        self._creationTime = ctimeInMillis
        self._updatedTime = ctimeInMillis
        cols = self._table.columns()[self.FILTER_COL + 1:]
        # We build lookup tables for the INFO and sample columns so they can
        # be easily found during conversion. For the sample columns we make
        # a dictionary mapping the sample name to a list of (sample name, col)
        # tuples for that sample.
        for col in cols:
            colName = col.get_name()
            if colName.startswith("INFO"):
                infoField = colName.split(".")[1]
                self._infoCols.append((infoField, col))
            else:
                if self._firstSamplePosition == -1:
                    # This must be a sample specific column
                    self._firstSamplePosition = col.get_position()
                sampleName, infoName = colName.split(".")
                if sampleName not in self._sampleCols:
                    self._sampleCols[sampleName] = []
                    self._sampleNames.append(sampleName)
                self._sampleCols[sampleName].append((infoName, col))

    def getNumVariants(self):
        return len(self._table)

    def convertInfoField(self, value):
        """
        Converts the specified value into an appropriate format for a protocol
        info field. Info fields are lists of strings.
        """
        if isinstance(value, tuple):
            ret = [str(x) for x in value]
        else:
            ret = [str(value)]
        return ret

    def convertVariant(self, row, sampleRowPositions):
        """
        Converts the specified wormtable row into a GAVariant object including
        the specified set of callSetIds.
        """
        variant = protocol.GAVariant()
        variant.created = self._creationTime
        variant.updated = self._updatedTime
        variant.variantSetId = self._variantSetId
        variant.referenceName = row[self.CHROM_COL]
        variant.start = row[self.POS_COL]
        names = row[self.ID_COL]
        if names is not None:
            variant.names = names.split(";")
        variant.referenceBases = row[self.REF_COL]
        variant.end = variant.start + len(variant.referenceBases)
        variant.id = "{0}.{1}".format(variant.variantSetId, row[0])
        alt = row[self.ALT_COL]
        if alt is not None:
            variant.alternateBases = alt.split(",")
        variant.info = {}
        for infoField, col in self._infoCols:
            pos = col.get_position()
            if row[pos] is not None:
                variant.info[infoField] = self.convertInfoField(row[pos])
        # All of the remaining values in the row correspond to Calls.
        for callSetId, rowPositions in sampleRowPositions.items():
            call = protocol.GACall()
            # Why do we have both of these? Wouldn't the ID be sufficient
            # to call back to searchCallSets to get more info?
            call.callSetId = callSetId
            call.callSetName = callSetId
            phaseset = None
            for (info, col), rowPosition in zip(self._sampleCols[callSetId],
                                                rowPositions):
                if info == self.GENOTYPE_LIKELIHOOD_NAME:
                    call.genotypeLikelihood = row[rowPosition]
                elif info == self.GENOTYPE_NAME:
                    genotype = row[rowPosition]
                elif info == self.PHASESET_NAME:
                    phaseset = row[rowPosition]
                else:
                    if row[rowPosition] is not None:
                        # Missing values are not included in the info array
                        call.info[info] = self.convertInfoField(
                            row[rowPosition])
            call.genotype, call.phaseset = convertVCFGenotype(genotype,
                                                              phaseset)
            variant.calls.append(call)
        return variant

    def getVariants(self, referenceName, startPosition, endPosition,
                    variantName, callSetIds):
        """
        Returns an iterator over the specified variants. The parameters
        correspond to the attributes of a GASearchVariantsRequest object.
        """
        # TODO encoding issues here? Wormtable (and most likely pysam and other
        # low level libraries) only deal with bytes. We'll need to adopt some
        # convention here for where we move from unicode to bytes.
        chrom = referenceName.encode()
        if callSetIds is None:
            readCols = self._table.columns()
            # These must be in the correct order, so we cannot use the keys of
            # self._sampleCols
            callSetIds = self._sampleNames
        else:
            readCols = self._table.columns()[:self._firstSamplePosition]
            for callSetId in callSetIds:
                try:
                    cols = [col for name, col in self._sampleCols[callSetId]]
                except KeyError:
                    raise exceptions.CallSetNotInVariantSetException(
                        callSetId, self._variantSetId)
                readCols.extend(cols)
        # Now we get the row positions for the sample columns
        sampleRowPositions = {}
        currentRowPosition = self._firstSamplePosition
        for callSetId in callSetIds:
            sampleRowPositionsList = []
            for col in self._sampleCols[callSetId]:
                sampleRowPositionsList.append(currentRowPosition)
                currentRowPosition += 1
            sampleRowPositions[callSetId] = sampleRowPositionsList
        if variantName is None:
            cursor = self._chromPosIndex.cursor(
                readCols, (chrom, startPosition), (chrom, endPosition))
            for row in cursor:
                yield self.convertVariant(row, sampleRowPositions)
        else:
            variantName = variantName.encode()  # TODO encoding?
            cursor = self._chromIdIndex.cursor(readCols, (chrom, variantName))
            for row in cursor:
                # TODO issues with having several names for a variant?
                # The result must still be within the range and must match
                # the specified name exactly. The cursor is positioned at
                # the first row >= the specified key.
                if (startPosition <= row[self.POS_COL] < endPosition and
                        row[self.ID_COL] == variantName):
                    yield self.convertVariant(row, sampleRowPositions)
                else:
                    break

    def getCreationTime(self):
        return self._creationTime

    def getUpdatedTime(self):
        return self._updatedTime

    def getCallSetIdMap(self):
        ret = {}
        for sampleName in self._sampleNames:
            callSetId = "{0}.{1}".format(self._variantSetId, sampleName)
            ret[callSetId] = CallSet(self, callSetId, sampleName)
        return ret

    def getCallSetIds(self):
        ret = []
        for sampleName in self._sampleNames:
            callSetId = "{0}.{1}".format(self._variantSetId, sampleName)
            ret.append(callSetId)
        return ret

    def getMetadata(self):
        """
        Returns a list of GAVariantSetMetadata objects for this variant set.
        """
        def buildMetadata(infoField, col):
            metadata = protocol.GAVariantSetMetadata()
            metadata.key = infoField
            metadata.value = ""
            metadata.id = ""
            # What are the encodings here? With the lack of any precise
            # definitions I'm just using wormtable values for now.
            metadata.type = col.get_type_name()
            metadata.number = str(col.get_num_elements())
            metadata.description = col.get_description()
            return metadata
        ret = []
        for infoField, col in self._infoCols:
            ret.append(buildMetadata(infoField, col))
        if len(self._sampleCols) > 0:
            # TODO this is pretty nasty, making a list just to take the head.
            sampleName = list(self._sampleCols.keys())[0]
            for infoField, col in self._sampleCols[sampleName]:
                if infoField != self.GENOTYPE_NAME:
                    ret.append(buildMetadata(infoField, col))
        return ret

    def toProtocolElement(self):
        """
        Converts this VariantSet into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.GAVariantSet()
        protocolElement.id = self._variantSetId
        protocolElement.datasetId = "NotImplemented"
        protocolElement.metadata = self.getMetadata()
        return protocolElement
