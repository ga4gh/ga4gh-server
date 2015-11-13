"""
Module responsible for translating variant data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import random
import hashlib

import pysam

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
import ga4gh.datamodel as datamodel


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


class CallSet(datamodel.DatamodelObject):
    """
    Class representing a CallSet. A CallSet basically represents the
    metadata associated with a single VCF sample column.
    """
    compoundIdClass = datamodel.CallSetCompoundId

    def toProtocolElement(self):
        """
        Returns the representation of this CallSet as the corresponding
        ProtocolElement.
        """
        variantSet = self.getParentContainer()
        gaCallSet = protocol.CallSet()
        gaCallSet.created = variantSet.getCreationTime()
        gaCallSet.updated = variantSet.getUpdatedTime()
        gaCallSet.id = self.getId()
        gaCallSet.name = self.getLocalId()
        gaCallSet.sampleId = self.getLocalId()
        gaCallSet.variantSetIds = [variantSet.getId()]
        return gaCallSet

    def getSampleName(self):
        """
        Returns the sample name for this CallSet.
        """
        return self.getLocalId()


class AbstractVariantSet(datamodel.DatamodelObject):
    """
    An abstract base class of a variant set
    """
    compoundIdClass = datamodel.VariantSetCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractVariantSet, self).__init__(parentContainer, localId)
        self._callSetIdMap = {}
        self._callSetNameMap = {}
        self._callSetIds = []
        self._creationTime = None
        self._updatedTime = None
        self._referenceSetId = ""

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

    def addCallSet(self, sampleName):
        """
        Adds a CallSet for the specified sample name.
        """
        callSet = CallSet(self, sampleName)
        callSetId = callSet.getId()
        self._callSetIdMap[callSetId] = callSet
        self._callSetNameMap[sampleName] = callSet
        self._callSetIds.append(callSetId)

    def getCallSets(self):
        """
        Returns the list of CallSets in this VariantSet.
        """
        return [self._callSetIdMap[id_] for id_ in self._callSetIds]

    def getNumCallSets(self):
        """
        Returns the number of CallSets in this variant set.
        """
        return len(self._callSetIds)

    def getCallSetByName(self, name):
        """
        Returns a CallSet with the specified name, or raises a
        CallSetNameNotFoundException if it does not exist.
        """
        if name not in self._callSetNameMap:
            raise exceptions.CallSetNameNotFoundException(name)
        return self._callSetNameMap[name]

    def getCallSetByIndex(self, index):
        """
        Returns the CallSet at the specfied index in this VariantSet.
        """
        return self._callSetIdMap[self._callSetIds[index]]

    def getCallSet(self, id_):
        """
        Returns a CallSet with the specified id, or raises a
        CallSetNotFoundException if it does not exist.
        """
        if id_ not in self._callSetIdMap:
            raise exceptions.CallSetNotFoundException(id_)
        return self._callSetIdMap[id_]

    def toProtocolElement(self):
        """
        Converts this VariantSet into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.VariantSet()
        protocolElement.id = self.getId()
        protocolElement.datasetId = self.getParentContainer().getId()
        protocolElement.referenceSetId = self._referenceSetId
        protocolElement.metadata = self.getMetadata()
        protocolElement.name = self.getLocalId()
        return protocolElement

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
        ret = protocol.Variant()
        ret.created = self._creationTime
        ret.updated = self._updatedTime
        ret.variantSetId = self.getId()
        return ret

    def getVariantId(self, gaVariant):
        """
        Returns an ID string suitable for the specified GA Variant
        object in this variant set.
        """
        md5 = self.hashVariant(gaVariant)
        compoundId = datamodel.VariantCompoundId(
            self.getCompoundId(), gaVariant.referenceName,
            gaVariant.start, md5)
        return str(compoundId)

    def getCallSetId(self, sampleName):
        """
        Returns the callSetId for the specified sampleName in this
        VariantSet.
        """
        compoundId = datamodel.CallSetCompoundId(
            self.getCompoundId(), sampleName)
        return str(compoundId)

    @classmethod
    def hashVariant(cls, gaVariant):
        """
        Produces an MD5 hash of the ga variant object to uniquely
        identify it
        """
        return hashlib.md5(
            gaVariant.referenceBases +
            str(tuple(gaVariant.alternateBases))).hexdigest()


class SimulatedVariantSet(AbstractVariantSet):
    """
    A variant set that doesn't derive from a data store.
    Used mostly for testing.
    """
    def __init__(
            self, parentContainer, localId, randomSeed=1, numCalls=1,
            variantDensity=1):
        super(SimulatedVariantSet, self).__init__(parentContainer, localId)
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
        # TODO Add simulated metadata.
        return ret

    def getVariant(self, compoundId):
        randomNumberGenerator = random.Random()
        start = int(compoundId.start)
        randomNumberGenerator.seed(self._randomSeed + start)
        variant = self.generateVariant(
            compoundId.referenceName, start, randomNumberGenerator)
        return variant

    def getVariants(self, referenceName, startPosition, endPosition,
                    callSetIds=None):
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(self._randomSeed)
        i = startPosition
        while i < endPosition:
            if randomNumberGenerator.random() < self._variantDensity:
                randomNumberGenerator.seed(self._randomSeed + i)
                yield self.generateVariant(
                    referenceName, i, randomNumberGenerator)
            i += 1

    def generateVariant(self, referenceName, position, randomNumberGenerator):
        """
        Generate a random variant for the specified position using the
        specified random number generator. This generator should be seeded
        with a value that is unique to this position so that the same variant
        will always be produced regardless of the order it is generated in.
        """
        variant = self._createGaVariant()
        variant.names = []
        variant.referenceName = referenceName
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
            call = protocol.Call()
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
        variant.id = self.getVariantId(variant)
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


class HtslibVariantSet(datamodel.PysamDatamodelMixin, AbstractVariantSet):
    """
    Class representing a single variant set backed by a directory of indexed
    VCF or BCF files.
    """
    def __init__(self, parentContainer, localId, dataDir, backend):
        super(HtslibVariantSet, self).__init__(parentContainer, localId)
        self._dataDir = dataDir
        self._setAccessTimes(dataDir)
        self._chromFileMap = {}
        self._metadata = None
        self._scanDataFiles(dataDir, ['*.bcf', '*.vcf.gz'])

    def _updateMetadata(self, variantFile):
        """
        Updates the metadata for his variant set based on the specified
        variant file, and ensures that it is consistent with already
        existing metadata.
        """
        metadata = self._getMetadataFromVcf(variantFile)
        if self._metadata is None:
            self._metadata = metadata
        else:
            if self._metadata != metadata:
                raise exceptions.InconsistentMetaDataException(
                    variantFile.filename)

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
                raise exceptions.InconsistentCallSetIdException(
                    variantFile.filename)

    def openFile(self, filename):
        return pysam.VariantFile(filename)

    def _addDataFile(self, filename):
        varFile = self.openFile(filename)
        if varFile.index is None:
            raise exceptions.NotIndexedException(filename)
        for chrom in varFile.index:
            # Unlike Tabix indices, CSI indices include all contigs defined
            # in the BCF header.  Thus we must test each one to see if
            # records exist or else they are likely to trigger spurious
            # overlapping errors.
            chrom, _, _ = self.sanitizeVariantFileFetch(chrom)
            if not isEmptyIter(varFile.fetch(chrom)):
                if chrom in self._chromFileMap:
                    raise exceptions.OverlappingVcfException(filename, chrom)
                self._updateMetadata(varFile)
                self._updateCallSetIds(varFile)
                self._chromFileMap[chrom] = filename
        varFile.close()

    def _convertGaCall(self, recordId, name, pysamCall, genotypeData):
        compoundId = self.getCallSetId(name)
        callSet = self.getCallSet(compoundId)
        call = protocol.Call()
        call.callSetId = callSet.getId()
        call.callSetName = callSet.getSampleName()
        call.sampleId = callSet.getSampleName()
        # TODO:
        # NOTE: THE FOLLOWING TWO LINES IS NOT THE INTENDED IMPLEMENTATION,
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
            if self.getCallSetId(name) in callSetIds:
                genotypeData = sampleData[sampleIterator].split(
                    ":")[0]  # REMOVAL
                variant.calls.append(self._convertGaCall(
                    record.id, name, call, genotypeData))  # REPLACE
            sampleIterator += 1  # REMOVAL
        variant.id = self.getVariantId(variant)
        return variant

    def getVariant(self, compoundId):
        if compoundId.referenceName in self._chromFileMap:
            varFileName = self._chromFileMap[compoundId.referenceName]
        else:
            raise exceptions.ObjectNotFoundException(compoundId)
        start = int(compoundId.start)
        referenceName, startPosition, endPosition = \
            self.sanitizeVariantFileFetch(
                compoundId.referenceName, start, start + 1)
        cursor = self.getFileHandle(varFileName).fetch(
            referenceName, startPosition, endPosition)
        for record in cursor:
            variant = self.convertVariant(record, self._callSetIds)
            if (record.start == start and
                    compoundId.md5 == self.hashVariant(variant)):
                return variant
            elif record.start > start:
                raise exceptions.ObjectNotFoundException()
        raise exceptions.ObjectNotFoundException(compoundId)

    def getVariants(self, referenceName, startPosition, endPosition,
                    callSetIds=None):
        """
        Returns an iterator over the specified variants. The parameters
        correspond to the attributes of a GASearchVariantsRequest object.
        """
        if callSetIds is None:
            callSetIds = self._callSetIds
        else:
            for callSetId in callSetIds:
                if callSetId not in self._callSetIds:
                    raise exceptions.CallSetNotInVariantSetException(
                        callSetId, self.getId())
        if referenceName in self._chromFileMap:
            varFileName = self._chromFileMap[referenceName]
            referenceName, startPosition, endPosition = \
                self.sanitizeVariantFileFetch(
                    referenceName, startPosition, endPosition)
            cursor = self.getFileHandle(varFileName).fetch(
                referenceName, startPosition, endPosition)
            for record in cursor:
                yield self.convertVariant(record, callSetIds)

    def getMetadata(self):
        return self._metadata

    def getMetadataId(self, metadata):
        """
        Returns the id of a metadata
        """
        return str(datamodel.VariantSetMetadataCompoundId(
            self.getCompoundId(), 'metadata:' + metadata.key))

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
                key, type_="String", number="1", value="", id_="",
                description=""):  # All input are strings
            metadata = protocol.VariantSetMetadata()
            metadata.key = key
            metadata.value = value
            metadata.type = type_
            metadata.number = number
            metadata.description = description
            if id_ == '':
                id_ = self.getMetadataId(metadata)
            metadata.id = id_
            return metadata

        ret = []
        header = varFile.header
        ret.append(buildMetadata(key="version", value=header.version))
        formats = header.formats.items()
        infos = header.info.items()
        # TODO: currently ALT field is not implemented through pysam
        # NOTE: contigs field is different between vcf files,
        # so it's not included in metadata
        # NOTE: filters in not included in metadata unless needed
        for prefix, content in [("FORMAT", formats), ("INFO", infos)]:
            for contentKey, value in content:
                attrs = dict(value.header.attrs)
                # TODO: refactor description at next pysam release
                # since description will be implemented as a member of
                # VariantMetadata
                description = attrs.get('Description', '').strip('"')
                key = "{0}.{1}".format(prefix, value.name)
                if key != "FORMAT.GT":
                    ret.append(buildMetadata(
                        key=key, type_=value.type,
                        number="{}".format(value.number),
                        description=description))
        return ret
