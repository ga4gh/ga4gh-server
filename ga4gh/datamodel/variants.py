"""
Module responsible for translating variant data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import atexit
import datetime
import glob
import json
import os
import random
import shutil
import tempfile
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

    def getSampleName(self):
        return self._sampleName

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


class AbstractVariantSet(datamodel.DatamodelObject):
    """
    An abstract base class of a variant set
    """
    def __init__(self, id_):
        super(AbstractVariantSet, self).__init__()
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
                    variantName=None, callSetIds=None):
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


def _cleanupHtslibsMess(indexDir):
    """
    Cleanup the mess that htslib has left behind with the index files.
    This is a temporary measure until we get a good interface for
    dealing with indexes for remote files.
    """
    shutil.rmtree(indexDir)


class HtslibVariantSet(datamodel.PysamSanitizer, AbstractVariantSet):
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
        numVariantFiles = 0
        for pattern in ['*.bcf', '*.vcf.gz']:
            for filename in glob.glob(os.path.join(vcfPath, pattern)):
                self._addFile(filename)
                numVariantFiles += 1
        if numVariantFiles == 0:
            raise exceptions.EmptyDirException()
        # This is a temporary workaround to allow us to use htslib's
        # facility for working with remote files. The urls.json is
        # definitely not a good idea and will be replaced later.
        # We make a temporary file for each process so that it
        # downloads its own copy and we are sure it's not overwriting
        # the copy of another process. We then register a cleanup
        # handler to get rid of these files on exit.
        urlSource = os.path.join(vcfPath, "urls.json")
        if os.path.exists(urlSource):
            with open(urlSource) as jsonFile:
                urls = json.load(jsonFile)["urls"]
            indexDir = tempfile.mkdtemp(prefix="htslib_mess.")
            cwd = os.getcwd()
            os.chdir(indexDir)
            for url in urls:
                self._addFile(url)
            os.chdir(cwd)
            atexit.register(_cleanupHtslibsMess, indexDir)

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

    def getCallSet(self, sampleName):
        """
        Returns the CallSet object for the specified sample name.
        """
        callSetId = self.getCallSetId(sampleName)
        return self._callSetIdMap[callSetId]

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

    def _addFile(self, filename):
        varFile = pysam.VariantFile(filename)
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
                self._chromFileMap[chrom] = varFile

    def _convertGaCall(self, recordId, name, pysamCall, genotypeData):
        callSet = self.getCallSet(name)
        call = protocol.GACall()
        call.callSetId = callSet.getId()
        call.callSetName = callSet.getSampleName()
        call.sampleId = callSet.getSampleName()
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
            if self.getCallSetId(name) in callSetIds:
                genotypeData = sampleData[sampleIterator].split(
                    ":")[0]  # REMOVAL
                variant.calls.append(self._convertGaCall(
                    record.id, name, call, genotypeData))  # REPLACE
            sampleIterator += 1  # REMOVAL
        return variant

    def getVariants(self, referenceName, startPosition, endPosition,
                    variantName=None, callSetIds=None):
        """
        Returns an iterator over the specified variants. The parameters
        correspond to the attributes of a GASearchVariantsRequest object.
        """
        if variantName is not None:
            raise exceptions.NotImplementedException(
                "Searching by variantName is not supported")
        # For v0.5.1, callSetIds=[] actually means return all callSets.
        # In v0.6+, callSetIds=[] means return no call sets, and
        # callSetIds=None means return all call sets. For forward
        # compatibility, we use the 0.6 interface for this function but
        # we translate back to the 0.5 interface while we support this.
        # TODO Remove this comment and workaround once we transition to
        # protocol version 0.6
        if callSetIds is None:
            callSetIds = []
        else:
            for callSetId in callSetIds:
                if callSetId not in self._callSetIds:
                    raise exceptions.CallSetNotInVariantSetException(
                        callSetId, self.getId())
        if len(callSetIds) == 0:
            callSetIds = self._callSetIds
        if referenceName in self._chromFileMap:
            varFile = self._chromFileMap[referenceName]
            referenceName, startPosition, endPosition = \
                self.sanitizeVariantFileFetch(
                    referenceName, startPosition, endPosition)
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
