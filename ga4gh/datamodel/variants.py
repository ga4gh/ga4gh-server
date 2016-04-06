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
import re
import os
import glob

import pysam

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
import ga4gh.datamodel as datamodel


def isUnspecified(str):
    """
    Checks whether a string is None or an
    empty string. Returns a boolean.
    """
    return str == "" or str is None


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

    def _createGaVariantAnnotation(self):
        """
        Convenience method to set the common fields in a GA VariantAnnotation
        object from this variant set.
        """
        ret = protocol.VariantAnnotation()
        ret.created = self._creationTime
        ret.updated = self._updatedTime
        ret.variantAnnotationSetId = self.getId()
        return ret

    def _createGaTranscriptEffect(self):
        """
        Convenience method to set the common fields in a GA TranscriptEffect
        object.
        """
        ret = protocol.TranscriptEffect()
        ret.created = self._creationTime
        ret.updated = self._updatedTime
        return ret

    def _createGaOntologyTermSo(self):
        """
        Convenience method to set the common fields in a GA OntologyTerm
        object for Sequence Ontology.
        """
        ret = protocol.OntologyTerm()
        ret.ontologySource = "Sequence Ontology"
        return ret

    def _createGaAlleleLocation(self):
        """
        Convenience method to set the common fields in a AlleleLocation
        object.
        """
        ret = protocol.AlleleLocation()
        ret.created = self._creationTime
        ret.updated = self._updatedTime
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


class SimulatedVariantAnnotationSet(AbstractVariantSet):
    """
    A variant annotation set that doesn't derive from a data store.
    Used mostly for testing.
    """
    def __init__(
            self, parentContainer, localId, variantSet, randomSeed=1):
        super(SimulatedVariantAnnotationSet, self).__init__(
            parentContainer, localId)
        self._variantSet = variantSet
        self._variantSetId = str(variantSet.getCompoundId())
        self._randomSeed = randomSeed
        # TODO refactor all the time stuff into protocol when schemas agree
        self._compoundId = datamodel.VariantAnnotationSetCompoundId(
            self.getCompoundId(), 'variantannotations')
        self._creationTime = datetime.datetime.now().isoformat() + "Z"
        self._updatedTime = datetime.datetime.now().isoformat() + "Z"
        # TODO make some reasonable connection to seqont

    def getAnalysis(self):
        analysis = protocol.Analysis()
        analysis.createDateTime = self._creationTime
        analysis.updateDateTime = self._updatedTime
        analysis.software.append("software")
        analysis.name = "name"
        analysis.description = "description"
        analysis.id = str(datamodel.VariantAnnotationSetAnalysisCompoundId(
            self._compoundId, "analysis"))
        return analysis

    def toProtocolElement(self):
        """
        Converts this VariantAnnotationSet into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.VariantAnnotationSet()
        protocolElement.variantSetId = self._variantSet.getId()
        protocolElement.id = self.getId()
        protocolElement.name = self.getLocalId()
        protocolElement.analysis = self.getAnalysis()
        return protocolElement

    def getNumVariantAnnotations(self):
        # TODO find out where this is used and how
        return 0

    def getVariantAnnotation(self, variant, randomNumberGenerator):
        ann = self.generateVariantAnnotation(
            variant, randomNumberGenerator)
        return ann

    def getVariantAnnotations(self, referenceName, start, end):
        for variant in self._variantSet.getVariants(referenceName, start, end):
            yield self.generateVariantAnnotation(variant)

    def generateVariantAnnotation(self, variant):
        """
        Generate a random variant annotation based on a given variant.
        This generator should be seeded with a value that is unique to the
        variant so that the same annotation will always be produced regardless
        of the order it is generated in.
        """
        # To make this reproducible, make a seed based on this
        # specific variant.
        seed = self._randomSeed + variant.start + variant.end
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(seed)
        ann = protocol.VariantAnnotation()
        ann.variantAnnotationSetId = str(self.getCompoundId())
        ann.variantId = variant.id
        ann.start = variant.start
        ann.end = variant.end
        ann.createDateTime = self._creationTime
        # make a transcript effect for each alternate base element
        # multiplied by a random integer (0,5)
        ann.transcriptEffects = []
        for base in variant.alternateBases * (
                randomNumberGenerator.randint(0, 5)):
            ann.transcriptEffects.append(self.generateTranscriptEffect(
                ann, base, randomNumberGenerator))
        ann.id = self.getVariantAnnotationId(variant, ann)
        return ann

    def getTranscriptEffectId(self, gaTranscriptEffect):
        effs = [eff.term for eff in gaTranscriptEffect.effects]
        return hashlib.md5(
            "{}\t{}\t{}\t{}".format(
                gaTranscriptEffect.alternateBases,
                gaTranscriptEffect.featureId,
                effs, gaTranscriptEffect.hgvsAnnotation)
            ).hexdigest()

    def getVariantAnnotationId(self, gaVariant, gaAnnotation):
        md5 = self.hashVariantAnnotation(gaVariant, gaAnnotation)
        compoundId = datamodel.VariantAnnotationCompoundId(
            self.getCompoundId(), gaVariant.referenceName,
            gaVariant.start, md5)
        return str(compoundId)

    def _addTranscriptEffectLocations(self, effect, ann):
        # TODO Make these valid HGVS values
        effect.hgvsAnnotation = protocol.HGVSAnnotation()
        effect.hgvsAnnotation.genomic = str(ann.start)
        effect.hgvsAnnotation.transcript = str(ann.start)
        effect.hgvsAnnotation.protein = str(ann.start)
        effect.proteinLocation = self._createGaAlleleLocation()
        effect.proteinLocation.start = ann.start
        effect.CDSLocation = self._createGaAlleleLocation()
        effect.CDSLocation.start = ann.start
        effect.cDNALocation = self._createGaAlleleLocation()
        effect.cDNALocation.start = ann.start
        return effect

    def _addTranscriptEffectId(self, effect):
        effect.id = str(self.getTranscriptEffectId(effect))
        return effect

    def _getRandomOntologyTerm(self, randomNumberGenerator):
        # TODO more mock options from simulated seqOnt?
        ontologyTuples = [("intron_variant", "SO:0001627"),
                          ("exon_variant", "SO:0001791")]
        term = protocol.OntologyTerm()
        ontologyTuple = randomNumberGenerator.choice(ontologyTuples)
        term.term, term.id = ontologyTuple[0], ontologyTuple[1]
        term.sourceName = "sequenceOntology"
        term.sourceVersion = "0"
        return term

    def _addTranscriptEffectOntologyTerm(self, effect, randomNumberGenerator):
        effect.effects.append(
            self._getRandomOntologyTerm(randomNumberGenerator))
        return effect

    def _generateAnalysisResult(self, effect, ann, randomNumberGenerator):
        # TODO make these sensible
        analysisResult = protocol.AnalysisResult()
        analysisResult.analysisId = "analysisId"
        analysisResult.result = "result string"
        analysisResult.score = randomNumberGenerator.randint(0, 100)
        return analysisResult

    def _addAnalysisResult(self, effect, ann, randomNumberGenerator):
        effect.analysisResults.append(
            self._generateAnalysisResult(
                effect, ann, randomNumberGenerator))
        return effect

    def generateTranscriptEffect(self, ann, alts, randomNumberGenerator):
        effect = self._createGaTranscriptEffect()
        effect.alternateBases = alts
        effect.effects = []
        effect.analysisResults = []
        # TODO how to make these featureIds sensical?
        effect.featureId = "E4TB33F"
        effect = self._addTranscriptEffectLocations(effect, ann)
        effect = self._addTranscriptEffectOntologyTerm(
            effect, randomNumberGenerator)
        effect = self._addTranscriptEffectOntologyTerm(
            effect, randomNumberGenerator)
        effect = self._addTranscriptEffectId(effect)
        effect = self._addAnalysisResult(effect, ann, randomNumberGenerator)
        return effect

    def hashVariantAnnotation(cls, gaVariant, gaVariantAnnotation):
        """
        Produces an MD5 hash of the gaVariant and gaVariantAnnotation objects
        """
        treffs = [treff.id for treff in gaVariantAnnotation.transcriptEffects]
        return hashlib.md5(
            "{}\t{}\t{}\t".format(
                gaVariant.referenceBases, tuple(gaVariant.alternateBases),
                treffs)
            ).hexdigest()


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
    def __init__(self, parentContainer, localId, dataDir, dataRepository):
        super(HtslibVariantSet, self).__init__(parentContainer, localId)
        self._dataDir = dataDir
        self._setAccessTimes(dataDir)
        self._chromFileMap = {}
        self._metadata = None
        self._patterns = ['*.bcf', '*.vcf.gz']
        self._scanDataFiles(self._dataDir, self._patterns)

    def _updateMetadata(self, variantFile):
        """
        Updates the metadata for his variant set based on the specified
        variant file
        """
        metadata = self._getMetadataFromVcf(variantFile)
        if self._metadata is None:
            self._metadata = metadata

    def _checkMetadata(self, variantFile):
        """
        Checks that metadata is consistent
        """
        metadata = self._getMetadataFromVcf(variantFile)
        if self._metadata is not None and self._metadata != metadata:
            raise exceptions.InconsistentMetaDataException(
                variantFile.filename)

    def _checkCallSetIds(self, variantFile):
        """
        Checks callSetIds for consistency
        """
        if len(self._callSetIdMap) > 0:
            callSetIds = set([
                self.getCallSetId(sample)
                for sample in variantFile.header.samples])
            if callSetIds != set(self._callSetIdMap.keys()):
                raise exceptions.InconsistentCallSetIdException(
                    variantFile.filename)

    def checkConsistency(self):
        """
        Perform consistency check on the variant set
        """
        filenames = self._getDataFilenames(self._dataDir, self._patterns)
        for filename in filenames:
            varFile = self.openFile(filename)
            for chrom in varFile.index:
                chrom, _, _ = self.sanitizeVariantFileFetch(chrom)
                if not isEmptyIter(varFile.fetch(chrom)):
                    self._checkMetadata(varFile)
                    self._checkCallSetIds(varFile)
            varFile.close()

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
        if len(self._callSetIdMap) == 0:
            for sample in variantFile.header.samples:
                self.addCallSet(sample)

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
                if isinstance(value, str):
                    value = value.split(',')
                variant.info[key] = _encodeValue(value)

        # NOTE: THE LABELED LINES SHOULD BE REMOVED ONCE PYSAM SUPPORTS
        # phaseset

        sampleData = record.__str__().split()[9:]  # REMOVAL
        variant.calls = []
        sampleIterator = 0  # REMOVAL
        if callSetIds is not None:
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
                description = value.description.strip('"')
                key = "{0}.{1}".format(prefix, value.name)
                if key != "FORMAT.GT":
                    ret.append(buildMetadata(
                        key=key, type_=value.type,
                        number="{}".format(value.number),
                        description=description))
        return ret

    def isAnnotated(self, path):
        # assumes that all files in the directory look like the first
        return self.hasAnnField(path) or self.hasConsequenceField(path)

    def hasAnnField(self, path):
        return 'ANN' in self._getHeaderItems(path)

    def hasConsequenceField(self, path):
        return 'CSQ' in self._getHeaderItems(path)

    def _getHeaderItems(self, path):
        pysamreader = self.openFile(
            glob.glob(os.path.join(path, "*.vcf.gz"))[0])
        return [x[0] for x in pysamreader.header.info.items()]


class HtslibVariantAnnotationSet(HtslibVariantSet):
    """
    Class representing a single variant annotation set backed by a directory of
    indexed VCF or BCF files.
    """
    def __init__(self, parentContainer, localId, dataDir,
                 backend, variantSet):
        super(HtslibVariantAnnotationSet, self).__init__(
            parentContainer, localId, dataDir, backend)
        self.compoundIdClass = datamodel.VariantAnnotationSetCompoundId
        self._variantSetId = variantSet.getCompoundId()
        self._variantSet = variantSet
        self._compoundId = datamodel.VariantAnnotationSetCompoundId(
            self.getCompoundId(), 'variantannotations')
        self._sequenceOntology = backend.getOntologyMap('sequence_ontology')
        self._creationTime = datetime.datetime.now().isoformat() + "Z"
        self._updatedTime = datetime.datetime.now().isoformat() + "Z"
        # Annotations are currently either from VEP or SNPEff. If they are
        # not from VEP we assume they are from SNPEff.
        # TODO Detect this more rigorously at import time and throw an
        # exception if we don't see the formats we're expecting.
        self._isVep = "VEP" in self.toProtocolElement().analysis.info
        self._isCsq = self.hasConsequenceField(dataDir)
        # Parse the annotation creation time out of the VCF header.
        # TODO Check this at import time, and raise an exception if the
        # time is not in the expected format.
        self._annotationCreatedDateTime = self._creationTime
        for r in self.getMetadata().records:
            # TODO handle more date formats
            if r.key == "created":
                self._annotationCreatedDateTime = datetime.datetime.strptime(
                    r.value, "%Y-%m-%d").isoformat() + "Z"

    def getVariantAnnotations(self, referenceName, startPosition, endPosition):
        """
        Generator for iterating through variant annotations in this
        variant annotation set.
        :param referenceName:
        :param startPosition:
        :param endPosition:
        :return: generator of protocol.VariantAnnotation
        """
        if referenceName in self._chromFileMap:
            varFileName = self._chromFileMap[referenceName]
            referenceName, startPosition, endPosition = \
                self.sanitizeVariantFileFetch(
                    referenceName, startPosition, endPosition)
            cursor = self.getFileHandle(varFileName).fetch(
                referenceName, startPosition, endPosition)
            transcriptConverter = self.convertTranscriptEffectSnpEff
            if self._isVep and not self._isCsq:
                transcriptConverter = self.convertTranscriptEffectVEP
            elif self._isCsq:
                transcriptConverter = self.convertTranscriptEffectCSQ
            for record in cursor:
                yield self.convertVariantAnnotation(
                    record, transcriptConverter)

    def convertLocation(self, pos):
        """
        Accepts a position string (start/length) and returns
        a GA4GH AlleleLocation with populated fields.
        :param pos:
        :return: protocol.AlleleLocation
        """
        if isUnspecified(pos):
            return None
        coordLen = pos.split('/')
        if len(coordLen) > 1:
            allLoc = self._createGaAlleleLocation()
            allLoc.start = int(coordLen[0]) - 1
            return allLoc
        return None

    def convertLocationHgvsC(self, hgvsc):
        """
        Accepts an annotation in HGVS notation and returns
        an AlleleLocation with populated fields.
        :param hgvsc:
        :return:
        """
        if isUnspecified(hgvsc):
            return None
        match = re.match(".*c.(\d+)(\D+)>(\D+)", hgvsc)
        if match:
            pos = int(match.group(1))
            if pos > 0:
                allLoc = self._createGaAlleleLocation()
                allLoc.start = pos - 1
                allLoc.referenceSequence = match.group(2)
                allLoc.alternateSequence = match.group(3)
                return allLoc
        return None

    def convertLocationHgvsP(self, hgvsp):
        """
        Accepts an annotation in HGVS notation and returns
        an AlleleLocation with populated fields.
        :param hgvsp:
        :return: protocol.AlleleLocation
        """
        if isUnspecified(hgvsp):
            return None
        match = re.match(".*p.(\D+)(\d+)(\D+)", hgvsp, flags=re.UNICODE)
        if match is not None:
            allLoc = self._createGaAlleleLocation()
            allLoc.referenceSequence = match.group(1)
            allLoc.start = int(match.group(2)) - 1
            allLoc.alternateSequence = match.group(3)
            return allLoc
        return None

    def addCDSLocation(self, effect, cdnaPos):
        hgvsC = effect.hgvsAnnotation.transcript
        if not isUnspecified(hgvsC):
            effect.CDSLocation = self.convertLocationHgvsC(hgvsC)
        if effect.CDSLocation is None:
            effect.CDSLocation = self.convertLocation(cdnaPos)
        else:
            # These are not stored in the VCF
            effect.CDSLocation.alternateSequence = None
            effect.CDSLocation.referenceSequence = None

    def addProteinLocation(self, effect, protPos):
        hgvsP = effect.hgvsAnnotation.protein
        if not isUnspecified(hgvsP):
            effect.proteinLocation = self.convertLocationHgvsP(hgvsP)
        if effect.proteinLocation is None:
            effect.proteinLocation = self.convertLocation(protPos)

    def addCDNALocation(self, effect, cdnaPos):
        hgvsC = effect.hgvsAnnotation.transcript
        effect.cDNALocation = self.convertLocation(cdnaPos)
        if self.convertLocationHgvsC(hgvsC):
            effect.cDNALocation.alternateSequence = \
                self.convertLocationHgvsC(hgvsC).alternateSequence
            effect.cDNALocation.referenceSequence = \
                self.convertLocationHgvsC(hgvsC).referenceSequence

    def addLocations(self, effect, protPos, cdnaPos):
        """
        Adds locations to a GA4GH transcript effect object
        by parsing HGVS annotation fields in concert with
        and supplied position values.
        :param effect: protocol.TranscriptEffect
        :param protPos: String representing protein position from VCF
        :param cdnaPos: String representing coding DNA location
        :return: effect protocol.TranscriptEffect
        """
        self.addCDSLocation(effect, cdnaPos)
        self.addCDNALocation(effect, cdnaPos)
        self.addProteinLocation(effect, protPos)
        return effect

    def convertTranscriptEffectCSQ(self, annStr, hgvsG):
        """
        Takes the consequence string of an annotated VCF using a
        CSQ field as opposed to ANN and returns an array of
        transcript effects.
        :param annStr: String
        :param hgvsG: String
        :return: [protocol.TranscriptEffect]
        """
        # Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|
        # CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|
        # DISTANCE|STRAND|SIFT|PolyPhen|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE

        (alt, gene, featureId, featureType, effects, cdnaPos,
         cdsPos, protPos, aminos, codons, existingVar, distance,
         strand, sift, polyPhen, motifName, motifPos,
         highInfPos, motifScoreChange) = annStr.split('|')
        terms = effects.split("&")
        transcriptEffects = []
        for term in terms:
            transcriptEffects.append(
                self._createCsqTranscriptEffect(
                    alt, term, protPos,
                    cdnaPos, featureId))
        return transcriptEffects

    def _createCsqTranscriptEffect(
            self, alt, term, protPos, cdnaPos, featureId):
        effect = self._createGaTranscriptEffect()
        effect.alternateBases = alt
        effect.effects = self.convertSeqOntology(term)
        effect.featureId = featureId
        effect.hgvsAnnotation = protocol.HGVSAnnotation()
        # These are not present in the data
        effect.hgvsAnnotation.genomic = None
        effect.hgvsAnnotation.transcript = None
        effect.hgvsAnnotation.protein = None
        self.addLocations(effect, protPos, cdnaPos)
        effect.id = self.getTranscriptEffectId(effect)
        effect.analysisResults = []
        return effect

    def convertTranscriptEffectVEP(self, annStr, hgvsG):
        """
        Takes the ANN string of a VEP generated VCF, splits it
        and returns a populated GA4GH transcript effect object.
        :param annStr: String
        :param hgvsG: String
        :return: effect protocol.TranscriptEffect
        """
        effect = self._createGaTranscriptEffect()
        (alt, effects, impact, symbol, geneName, featureType,
         featureId, trBiotype, exon, intron, hgvsC, hgvsP,
         cdnaPos, cdsPos, protPos, aminos, codons,
         existingVar, distance, strand, symbolSource,
         hgncId, hgvsOffset) = annStr.split('|')
        effect.alternateBases = alt
        effect.effects = self.convertSeqOntology(effects)
        effect.featureId = featureId
        effect.hgvsAnnotation = protocol.HGVSAnnotation()
        effect.hgvsAnnotation.genomic = hgvsG
        effect.hgvsAnnotation.transcript = hgvsC
        effect.hgvsAnnotation.protein = hgvsP
        self.addLocations(effect, protPos, cdnaPos)
        effect.id = self.getTranscriptEffectId(effect)
        effect.analysisResults = []
        return effect

    def convertTranscriptEffectSnpEff(self, annStr, hgvsG):
        """
        Takes the ANN string of a SnpEff generated VCF, splits it
        and returns a populated GA4GH transcript effect object.
        :param annStr: String
        :param hgvsG: String
        :return: effect protocol.TranscriptEffect()
        """
        effect = self._createGaTranscriptEffect()
        # SnpEff and VEP don't agree on this :)
        (alt, effects, impact, geneName, geneId, featureType,
            featureId, trBiotype, rank, hgvsC, hgvsP, cdnaPos,
            cdsPos, protPos, distance, errsWarns) = annStr.split('|')
        effect.alternateBases = alt
        effect.effects = self.convertSeqOntology(effects)
        effect.featureId = featureId
        effect.hgvsAnnotation = protocol.HGVSAnnotation()
        effect.hgvsAnnotation.genomic = hgvsG
        effect.hgvsAnnotation.transcript = hgvsC
        effect.hgvsAnnotation.protein = hgvsP
        self.addLocations(effect, protPos, cdnaPos)
        effect.id = self.getTranscriptEffectId(effect)
        effect.analysisResults = []
        return effect

    def convertSeqOntology(self, seqOntStr):
        """
        Splits a string of sequence ontology effects and creates
        an ontology term record for each, which are built into
        an array of return soTerms.
        :param seqOntStr:
        :return: [protocol.OntologyTerm]
        """
        seqOntTerms = seqOntStr.split('&')
        soTerms = []
        for soName in seqOntTerms:
            so = self._createGaOntologyTermSo()
            so.term = soName
            if self._sequenceOntology is not None:
                so.id = self._sequenceOntology.getId(soName, "")
            soTerms.append(so)
        return soTerms

    def convertVariantAnnotation(self, record, transcriptConverter):
        """
        Converts the specfied pysam variant record into a GA4GH variant
        annotation object using the specified function to convert the
        transcripts.
        """
        variant = self.convertVariant(record, None)
        annotation = self._createGaVariantAnnotation()
        annotation.start = variant.start
        annotation.end = variant.end
        annotation.createDateTime = self._annotationCreatedDateTime
        annotation.variantId = variant.id
        # Convert annotations from INFO field into TranscriptEffect
        transcriptEffects = []
        hgvsG = record.info.get('HGVS.g'.encode())
        if transcriptConverter != self.convertTranscriptEffectCSQ:
            annotations = record.info.get('ANN'.encode())
            transcriptEffects = self._convertAnnotations(
                annotations, variant, hgvsG, transcriptConverter)
        else:
            annotations = record.info.get('CSQ'.encode())
            transcriptEffects = self.convertTranscriptEffectCSQ(
                annotations[0], hgvsG)

        annotation.transcriptEffects = transcriptEffects
        annotation.id = self.getVariantAnnotationId(variant, annotation)
        return annotation

    def _convertAnnotations(
            self, annotations, variant, hgvsG, transcriptConverter):
        transcriptEffects = []
        if annotations is not None:
            for index, ann in enumerate(annotations):
                altshgvsG = ""
                if hgvsG is not None:
                    # The HGVS.g field contains an element for
                    # each alternate allele
                    altshgvsG = hgvsG[index % len(variant.alternateBases)]
                transcriptEffects.append(
                    transcriptConverter(ann, altshgvsG))
        return transcriptEffects

    def getVariantAnnotationId(self, gaVariant, gaAnnotation):
        """
        Produces a stringified compoundId representing a variant
        annotation.
        :param gaVariant:   protocol.Variant
        :param gaAnnotation: protocol.VariantAnnotation
        :return:  compoundId String
        """
        md5 = self.hashVariantAnnotation(gaVariant, gaAnnotation)
        compoundId = datamodel.VariantAnnotationCompoundId(
            self.getCompoundId(), gaVariant.referenceName,
            gaVariant.start, md5)
        return str(compoundId)

    def getVariantId(self, gaVariant):
        """
        Produces a variant ID for a variant annotated within this
        variant annotation set.
        :param gaVariant: protocol.Variant
        :return:  compoundId String
        """
        md5 = self.hashVariant(gaVariant)
        compoundId = datamodel.VariantCompoundId(
            self._variantSetId, gaVariant.referenceName,
            gaVariant.start, md5)
        return str(compoundId)

    def getAnalysis(self):
        """
        Assembles metadata within the VCF header into a GA4GH
        Analysis object.
        :return: protocol.Analysis
        """
        metadata = self.getMetadata()
        analysis = protocol.Analysis()
        formats = metadata.formats.items()
        infos = metadata.info.items()
        for prefix, content in [("FORMAT", formats), ("INFO", infos)]:
            for contentKey, value in content:
                key = "{0}.{1}".format(prefix, value.name)
                if key not in analysis.info:
                    analysis.info[key] = []
                if value.description is not None:
                    analysis.info[key].append(value.description)
        analysis.createDateTime = self._creationTime
        analysis.updateDateTime = self._updatedTime
        for r in metadata.records:
            # Don't add a key to info if there's nothing in the value
            if r.value is not None:
                if r.key not in analysis.info:
                    analysis.info[r.key] = []
                analysis.info[r.key].append(str(r.value))

            if r.key == "created":
                # TODO handle more date formats
                analysis.createDateTime = datetime.datetime.strptime(
                    r.value, "%Y-%m-%d").isoformat() + "Z"
            if r.key == "software":
                analysis.software.append(r.value)
            if r.key == "name":
                analysis.name = r.value
            if r.key == "description":
                analysis.description = r.value
        analysis.id = str(datamodel.VariantAnnotationSetAnalysisCompoundId(
            self._compoundId, "analysis"))
        return analysis

    def _getMetadataFromVcf(self, varFile):
        header = varFile.header
        return header

    def toProtocolElement(self):
        """
        Converts this VariantSet into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.VariantAnnotationSet()
        protocolElement.id = self.getId()
        protocolElement.variantSetId = str(self._variantSetId)
        protocolElement.name = self.getLocalId()
        protocolElement.analysis = self.getAnalysis()
        return protocolElement

    def getTranscriptEffectId(self, gaTranscriptEffect):
        effs = [eff.term for eff in gaTranscriptEffect.effects]
        return hashlib.md5(
            "{}\t{}\t{}\t{}".format(
                gaTranscriptEffect.alternateBases,
                gaTranscriptEffect.featureId,
                effs, gaTranscriptEffect.hgvsAnnotation)
            ).hexdigest()

    def hashVariantAnnotation(cls, gaVariant, gaVariantAnnotation):
        """
        Produces an MD5 hash of the gaVariant and gaVariantAnnotation objects
        """
        treffs = [treff.id for treff in gaVariantAnnotation.transcriptEffects]
        return hashlib.md5(
            "{}\t{}\t{}\t".format(
                gaVariant.referenceBases, tuple(gaVariant.alternateBases),
                treffs)
            ).hexdigest()
