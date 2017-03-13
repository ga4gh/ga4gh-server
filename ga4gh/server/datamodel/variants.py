"""
Module responsible for translating variant data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import glob
import hashlib
import json
import os
import random
import re

import pysam

import ga4gh.server.exceptions as exceptions
import ga4gh.server.datamodel as datamodel

import ga4gh.schemas.pb as pb
import ga4gh.schemas.protocol as protocol


ANNOTATIONS_VEP_V82 = "VEP_v82"
ANNOTATIONS_VEP_V77 = "VEP_v77"
ANNOTATIONS_SNPEFF = "SNPEff"


# Utility functions for module

def isUnspecified(str):
    """
    Checks whether a string is None or an
    empty string. Returns a boolean.
    """
    return str == "" or str is None


_nothing = object()


def isEmptyIter(it):
    """Return True iff the iterator is empty or exhausted"""
    return next(it, _nothing) is _nothing


class CallSet(datamodel.DatamodelObject):
    """
    Class representing a CallSet. A CallSet basically represents the
    metadata associated with a single VCF sample column.
    """
    compoundIdClass = datamodel.CallSetCompoundId

    def __init__(self, parentContainer, localId):
        super(CallSet, self).__init__(parentContainer, localId)
        self._info = {}
        self._biosampleId = None

    def populateFromRow(self, callSetRecord):
        """
        Populates this CallSet from the specified DB row.
        """
        self._biosampleId = callSetRecord.biosampleid
        self.setAttributesJson(callSetRecord.attributes)

    def toProtocolElement(self):
        """
        Returns the representation of this CallSet as the corresponding
        ProtocolElement.
        """
        variantSet = self.getParentContainer()
        gaCallSet = protocol.CallSet(
            biosample_id=self.getBiosampleId())
        if variantSet.getCreationTime():
            gaCallSet.created = variantSet.getCreationTime()
        if variantSet.getUpdatedTime():
            gaCallSet.updated = variantSet.getUpdatedTime()
        gaCallSet.id = self.getId()
        gaCallSet.name = self.getLocalId()
        gaCallSet.variant_set_ids.append(variantSet.getId())
        self.serializeAttributes(gaCallSet)
        return gaCallSet

    def getBiosampleId(self):
        """
        Returns the biosampleId for this CallSet.
        """
        return self._biosampleId

    def setBiosampleId(self, biosampleId):
        """
        Set the biosampleId for the current sample.
        """
        self._biosampleId = biosampleId

    def getSampleName(self):
        """
        Returns the sample name for this CallSet.
        """
        return self.getLocalId()

    def getInfo(self):
        """
        Returns info map for this CallSet
        """
        return self._info


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
        self._callSetIdToIndex = {}
        self._creationTime = None
        self._updatedTime = None
        self._referenceSet = None
        self._metadata = []
        self._variantAnnotationSetIds = []
        self._variantAnnotationSetIdMap = {}

    def addVariantAnnotationSet(self, variantAnnotationSet):
        """
        Adds the specified variantAnnotationSet to this dataset.
        """
        id_ = variantAnnotationSet.getId()
        self._variantAnnotationSetIdMap[id_] = variantAnnotationSet
        self._variantAnnotationSetIds.append(id_)

    def getVariantAnnotationSets(self):
        """
        Returns the list of VariantAnnotationSets in this dataset
        """
        return [
            self._variantAnnotationSetIdMap[id_] for id_ in
            self._variantAnnotationSetIds]

    def getVariantAnnotationSet(self, id_):
        """
        Returns the AnnotationSet in this dataset with the specified 'id'
        """
        if id_ not in self._variantAnnotationSetIdMap:
            raise exceptions.AnnotationSetNotFoundException(id_)
        return self._variantAnnotationSetIdMap[id_]

    def getNumVariantAnnotationSets(self):
        """
        Returns the number of variant annotation sets in this dataset.
        """
        return len(self._variantAnnotationSetIds)

    def getVariantAnnotationSetByIndex(self, index):
        """
        Returns the variant annotation set at the specified index in this
        dataset.
        """
        return self._variantAnnotationSetIdMap[
            self._variantAnnotationSetIds[index]]

    def setReferenceSet(self, referenceSet):
        """
        Sets the ReferenceSet for this VariantSet to the specified value.
        """
        self._referenceSet = referenceSet

    def getReferenceSet(self):
        """
        Returns the reference set associated with this VariantSet.
        """
        return self._referenceSet

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

    def addCallSet(self, callSet):
        """
        Adds the specfied CallSet to this VariantSet.
        """
        callSetId = callSet.getId()
        self._callSetIdMap[callSetId] = callSet
        self._callSetNameMap[callSet.getLocalId()] = callSet
        self._callSetIds.append(callSetId)
        self._callSetIdToIndex[callSet.getId()] = len(self._callSetIds) - 1

    def addCallSetFromName(self, sampleName):
        """
        Adds a CallSet for the specified sample name.
        """
        callSet = CallSet(self, sampleName)
        self.addCallSet(callSet)

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

    def getMetadata(self):
        """
        Returns Metdata associated with this VariantSet
        """
        return self._metadata

    def toProtocolElement(self):
        """
        Converts this VariantSet into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.VariantSet()
        protocolElement.id = self.getId()
        protocolElement.dataset_id = self.getParentContainer().getId()
        protocolElement.reference_set_id = self._referenceSet.getId()
        protocolElement.metadata.extend(self.getMetadata())
        protocolElement.dataset_id = self.getParentContainer().getId()
        protocolElement.reference_set_id = self._referenceSet.getId()
        protocolElement.name = self.getLocalId()
        self.serializeAttributes(protocolElement)
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
        if self._creationTime:
            ret.created = self._creationTime
        if self._updatedTime:
            ret.updated = self._updatedTime
        ret.variant_set_id = self.getId()
        return ret

    def getVariantId(self, gaVariant):
        """
        Returns an ID string suitable for the specified GA Variant
        object in this variant set.
        """
        md5 = self.hashVariant(gaVariant)
        compoundId = datamodel.VariantCompoundId(
            self.getCompoundId(), gaVariant.reference_name,
            str(gaVariant.start), md5)
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
        Produces an MD5 hash of the ga variant object to distinguish
        it from other variants at the same genomic coordinate.
        """
        hash_str = gaVariant.reference_bases + \
            str(tuple(gaVariant.alternate_bases))
        return hashlib.md5(hash_str).hexdigest()


class SimulatedVariantSet(AbstractVariantSet):
    """
    A variant set that doesn't derive from a data store.
    Used mostly for testing.
    """
    def __init__(
            self, parentContainer, referenceSet, localId, randomSeed=1,
            numCalls=1, variantDensity=1):
        super(SimulatedVariantSet, self).__init__(parentContainer, localId)
        self._referenceSet = referenceSet
        self._randomSeed = randomSeed
        self._numCalls = numCalls
        for i in range(numCalls):
            callSetName = "simCallSet_{}".format(i)
            self.addCallSetFromName(callSetName)
            callSet = self.getCallSetByName(callSetName)
            # build up infos of increasing size
            for j in range(i):
                callSet._info["key_{}".format(j)] = "value_{}".format(j)
        self._variantDensity = variantDensity
        self._metadata = self._createMetaData()
        now = protocol.convertDatetime(datetime.datetime.now())
        self._creationTime = now
        self._updatedTime = now

    def _createMetaData(self):
        metadata_1 = protocol.VariantSetMetadata()
        metadata_1.key = "version"
        metadata_1.value = "VCFv4.1"
        metadata_1.type = "String"
        metadata_1.number = "1"
        metadata_1.description = ""
        metadata_1.id = str(datamodel.VariantSetMetadataCompoundId(
            self.getCompoundId(), 'metadata:' + metadata_1.key))

        metadata_2 = protocol.VariantSetMetadata()
        metadata_2.key = "INFO.FIELD1"
        metadata_2.value = ""
        metadata_2.type = "String"
        metadata_2.number = ""
        metadata_2.description = "FIELD1 description"
        metadata_2.id = str(datamodel.VariantSetMetadataCompoundId(
            self.getCompoundId(), 'metadata:' + metadata_2.key))

        metadata_3 = protocol.VariantSetMetadata()
        metadata_3.key = "INFO.FIELD2"
        metadata_3.value = ""
        metadata_3.type = "Integer"
        metadata_3.number = "1"
        metadata_3.description = "FIELD2 description"
        metadata_3.id = str(datamodel.VariantSetMetadataCompoundId(
            self.getCompoundId(), 'metadata:' + metadata_3.key))

        return [metadata_1, metadata_2, metadata_3]

    def getNumVariants(self):
        return 0

    def getVariant(self, compoundId):
        randomNumberGenerator = random.Random()
        start = int(compoundId.start)
        randomNumberGenerator.seed(self._randomSeed + start)
        variant = self.generateVariant(
            compoundId.reference_name, start, randomNumberGenerator)
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
        variant.reference_name = referenceName
        variant.start = position
        variant.end = position + 1  # SNPs only for now
        bases = ["A", "C", "G", "T"]
        ref = randomNumberGenerator.choice(bases)
        variant.reference_bases = ref
        alt = randomNumberGenerator.choice(
            [base for base in bases if base != ref])
        variant.alternate_bases.append(alt)
        randChoice = randomNumberGenerator.randint(0, 2)
        if randChoice == 0:
            variant.filters_applied = False
        elif randChoice == 1:
            variant.filters_applied = True
            variant.filters_passed = True
        else:
            variant.filters_applied = True
            variant.filters_passed = False
            variant.filters_failed.append('q10')
        for callSet in self.getCallSets():
            call = variant.calls.add()
            call.call_set_id = callSet.getId()
            # for now, the genotype is either [0,1], [1,1] or [1,0] with equal
            # probability; probably will want to do something more
            # sophisticated later.
            randomChoice = randomNumberGenerator.choice(
                [[0, 1], [1, 0], [1, 1]])
            call.genotype.extend(randomChoice)
            # TODO What is a reasonable model for generating these likelihoods?
            # Are these log-scaled? Spec does not say.
            call.genotype_likelihood.extend([-100, -100, -100])
        variant.id = self.getVariantId(variant)
        return variant


class HtslibVariantSet(datamodel.PysamDatamodelMixin, AbstractVariantSet):
    """
    Class representing a single variant set backed by a directory of indexed
    VCF or BCF files.
    """
    def __init__(self, parentContainer, localId):
        super(HtslibVariantSet, self).__init__(parentContainer, localId)
        self._chromFileMap = {}
        self._metadata = None

    def isAnnotated(self):
        """
        Returns True if there is a VariantAnnotationSet associated with this
        VariantSet.
        """
        return len(self._variantAnnotationSetIdMap) > 0

    def getReferenceToDataUrlIndexMap(self):
        """
        Returns the map of Reference names to the (dataUrl, indexFile) pairs.
        """
        return self._chromFileMap

    def getDataUrlIndexPairs(self):
        """
        Returns the set of (dataUrl, indexFile) pairs.
        """
        return set(self._chromFileMap.values())

    def populateFromRow(self, variantSetRecord):
        """
        Populates this VariantSet from the specified DB row.
        """
        self._created = variantSetRecord.created
        self._updated = variantSetRecord.updated
        self.setAttributesJson(variantSetRecord.attributes)
        self._chromFileMap = {}
        # We can't load directly as we want tuples to be stored
        # rather than lists.
        for key, value in json.loads(variantSetRecord.dataurlindexmap).items():
            self._chromFileMap[key] = tuple(value)
        self._metadata = []
        for jsonDict in json.loads(variantSetRecord.metadata):
            metadata = protocol.fromJson(json.dumps(jsonDict),
                                         protocol.VariantSetMetadata)
            self._metadata.append(metadata)

    def populateFromFile(self, dataUrls, indexFiles):
        """
        Populates this variant set using the specified lists of data
        files and indexes. These must be in the same order, such that
        the jth index file corresponds to the jth data file.
        """
        assert len(dataUrls) == len(indexFiles)
        for dataUrl, indexFile in zip(dataUrls, indexFiles):
            varFile = pysam.VariantFile(dataUrl, index_filename=indexFile)
            try:
                self._populateFromVariantFile(varFile, dataUrl, indexFile)
            finally:
                varFile.close()

    def populateFromDirectory(self, vcfDirectory):
        """
        Populates this VariantSet by examing all the VCF files in the
        specified directory. This is mainly used for as a convenience
        for testing purposes.
        """
        pattern = os.path.join(vcfDirectory, "*.vcf.gz")
        dataFiles = []
        indexFiles = []
        for vcfFile in glob.glob(pattern):
            dataFiles.append(vcfFile)
            indexFiles.append(vcfFile + ".tbi")
        self.populateFromFile(dataFiles, indexFiles)

    def getVcfHeaderReferenceSetName(self):
        """
        Returns the name of the reference set from the VCF header.
        """
        # TODO implemenent
        return None

    def checkConsistency(self):
        """
        Perform consistency check on the variant set
        """
        for referenceName, (dataUrl, indexFile) in self._chromFileMap.items():
            varFile = pysam.VariantFile(dataUrl, index_filename=indexFile)
            try:
                for chrom in varFile.index:
                    chrom, _, _ = self.sanitizeVariantFileFetch(chrom)
                    if not isEmptyIter(varFile.fetch(chrom)):
                        self._checkMetadata(varFile)
                        self._checkCallSetIds(varFile)
            finally:
                varFile.close()

    def _populateFromVariantFile(self, varFile, dataUrl, indexFile):
        """
        Populates the instance variables of this VariantSet from the specified
        pysam VariantFile object.
        """
        if varFile.index is None:
            raise exceptions.NotIndexedException(dataUrl)
        for chrom in varFile.index:
            # Unlike Tabix indices, CSI indices include all contigs defined
            # in the BCF header.  Thus we must test each one to see if
            # records exist or else they are likely to trigger spurious
            # overlapping errors.
            chrom, _, _ = self.sanitizeVariantFileFetch(chrom)
            if not isEmptyIter(varFile.fetch(chrom)):
                if chrom in self._chromFileMap:
                    raise exceptions.OverlappingVcfException(dataUrl, chrom)
            self._chromFileMap[chrom] = dataUrl, indexFile
        self._updateMetadata(varFile)
        self._updateCallSetIds(varFile)
        self._updateVariantAnnotationSets(varFile, dataUrl)

    def _updateVariantAnnotationSets(self, variantFile, dataUrl):
        """
        Updates the variant annotation set associated with this variant using
        information in the specified pysam variantFile.
        """
        # TODO check the consistency of this between VCF files.
        if not self.isAnnotated():
            annotationType = None
            for record in variantFile.header.records:
                if record.type == "GENERIC":
                    if record.key == "SnpEffVersion":
                        annotationType = ANNOTATIONS_SNPEFF
                    elif record.key == "VEP":
                        version = record.value.split()[0]
                        # TODO we need _much_ more sophisticated processing
                        # of VEP versions here. When do they become
                        # incompatible?
                        if version == "v82":
                            annotationType = ANNOTATIONS_VEP_V82
                        elif version == "v77":
                            annotationType = ANNOTATIONS_VEP_V77
                        else:
                            # TODO raise a proper typed exception there with
                            # the file name as an argument.
                            raise ValueError(
                                "Unsupported VEP version {} in '{}'".format(
                                    version, dataUrl))
            if annotationType is None:
                infoKeys = variantFile.header.info.keys()
                if 'CSQ' in infoKeys or 'ANN' in infoKeys:
                    # TODO likewise, we want a properly typed exception that
                    # we can throw back to the repo manager UI and display
                    # as an import error.
                    raise ValueError(
                        "Unsupported annotations in '{}'".format(dataUrl))
            if annotationType is not None:
                vas = HtslibVariantAnnotationSet(self, self.getLocalId())
                vas.populateFromFile(variantFile, annotationType)
                self.addVariantAnnotationSet(vas)

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
                self.addCallSetFromName(sample)

    def openFile(self, dataUrlIndexFilePair):
        dataUrl, indexFile = dataUrlIndexFilePair
        return pysam.VariantFile(dataUrl, index_filename=indexFile)

    def _convertGaCall(self, callSet, pysamCall):
        phaseset = None
        if pysamCall.phased:
            phaseset = str(pysamCall.phased)
        genotypeLikelihood = []
        info = {}
        for key, value in pysamCall.iteritems():
            if key == 'GL' and value is not None:
                genotypeLikelihood = list(value)
            elif key != 'GT':
                info[key] = protocol.encodeValue(value)
        call = protocol.Call()
        call.call_set_name = callSet.getSampleName()
        call.call_set_id = callSet.getId()
        call.genotype.extend(list(pysamCall.allele_indices))
        call.phaseset = pb.string(phaseset)
        call.genotype_likelihood.extend(genotypeLikelihood)
        for key in info:
            call.attributes.attr[key].values.extend(info[key])
        return call

    def convertVariant(self, record, callSetIds):
        """
        Converts the specified pysam variant record into a GA4GH Variant
        object. Only calls for the specified list of callSetIds will
        be included.
        """
        variant = self._createGaVariant()
        variant.reference_name = record.contig
        if record.id is not None:
            variant.names.extend(record.id.split(';'))
        variant.start = record.start          # 0-based inclusive
        variant.end = record.stop             # 0-based exclusive
        variant.reference_bases = record.ref
        if record.alts is not None:
            variant.alternate_bases.extend(list(record.alts))
        filterKeys = record.filter.keys()
        if len(filterKeys) == 0:
            variant.filters_applied = False
        else:
            variant.filters_applied = True
            if len(filterKeys) == 1 and filterKeys[0] == 'PASS':
                variant.filters_passed = True
            else:
                variant.filters_passed = False
                variant.filters_failed.extend(filterKeys)
        # record.qual is also available, when supported by GAVariant.
        for key, value in record.info.iteritems():
            if value is None:
                continue
            if key == 'SVTYPE':
                variant.variant_type = value
            elif key == 'SVLEN':
                variant.svlen = int(value[0])
            elif key == 'CIPOS':
                variant.cipos.extend(value)
            elif key == 'CIEND':
                variant.ciend.extend(value)
            elif isinstance(value, str):
                value = value.split(',')
            protocol.setAttribute(
                variant.attributes.attr[key].values, value)
        for callSetId in callSetIds:
            callSet = self.getCallSet(callSetId)
            pysamCall = record.samples[str(callSet.getSampleName())]
            variant.calls.add().CopyFrom(
                self._convertGaCall(callSet, pysamCall))
        variant.id = self.getVariantId(variant)
        return variant

    def getVariant(self, compoundId):
        if compoundId.reference_name in self._chromFileMap:
            varFileName = self._chromFileMap[compoundId.reference_name]
        else:
            raise exceptions.ObjectNotFoundException(compoundId)
        start = int(compoundId.start)
        referenceName, startPosition, endPosition = \
            self.sanitizeVariantFileFetch(
                compoundId.reference_name, start, start + 1)
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

    def getPysamVariants(self, referenceName, startPosition, endPosition):
        """
        Returns an iterator over the pysam VCF records corresponding to the
        specified query.
        """
        if referenceName in self._chromFileMap:
            varFileName = self._chromFileMap[referenceName]
            referenceName, startPosition, endPosition = \
                self.sanitizeVariantFileFetch(
                    referenceName, startPosition, endPosition)
            cursor = self.getFileHandle(varFileName).fetch(
                referenceName, startPosition, endPosition)
            for record in cursor:
                yield record

    def getVariants(self, referenceName, startPosition, endPosition,
                    callSetIds=[]):
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
        for record in self.getPysamVariants(
                referenceName, startPosition, endPosition):
            yield self.convertVariant(record, callSetIds)

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
        filters = header.filters.items()
        # TODO: currently ALT field is not implemented through pysam
        # NOTE: contigs field is different between vcf files,
        # so it's not included in metadata
        for prefix, content in [("FORMAT", formats), ("INFO", infos),
                                ("FILTER", filters)]:
            for contentKey, value in content:
                description = value.description.strip('"')
                key = "{0}.{1}".format(prefix, value.name)
                if prefix == "FILTER":
                    ret.append(buildMetadata(
                        key=key,
                        description=description))
                elif key != "FORMAT.GT":
                    ret.append(buildMetadata(
                        key=key, type_=value.type,
                        number="{}".format(value.number),
                        description=description))
        return ret

#############################################

# Variant Annotations.

#############################################


class AbstractVariantAnnotationSet(datamodel.DatamodelObject):
    """
    Class representing a variant annotation set derived from an
    annotated variant set.
    """
    compoundIdClass = datamodel.VariantAnnotationSetCompoundId

    def __init__(self, variantSet, localId):
        super(AbstractVariantAnnotationSet, self).__init__(variantSet, localId)
        self._variantSet = variantSet
        self._ontology = None
        self._analysis = None
        self._creationTime = ''
        self._updatedTime = ''

    def setOntology(self, ontology):
        """
        Sets the Ontology used in this VariantAnnotationSet to
        translate sequence ontology term names into IDs to the
        specified value.
        """
        self._ontology = ontology

    def getCreationTime(self):
        """
        Returns the creation time for this VariantAnnotationSet
        """
        return self._creationTime

    def getUpdatedTime(self):
        """
        Returns the update time for this VariantAnnotationSet
        """
        return self._updatedTime

    def getOntology(self):
        """
        Returns the ontology term map used in this VariantAnnotationSet.
        """
        return self._ontology

    def getAnalysis(self):
        """
        Returns the Analysis object associated with this VariantAnnotationSet.
        """
        return self._analysis

    def getVariantSet(self):
        """
        Returns the VariantSet that this VariantAnnotationSet refers to.
        """
        return self._variantSet

    def _createGaVariantAnnotation(self):
        """
        Convenience method to set the common fields in a GA VariantAnnotation
        object from this variant set.
        """
        ret = protocol.VariantAnnotation()
        ret.created = self._creationTime
        ret.variant_annotation_set_id = self.getId()
        return ret

    def _createGaTranscriptEffect(self):
        """
        Convenience method to set the common fields in a GA TranscriptEffect
        object.
        """
        ret = protocol.TranscriptEffect()
        return ret

    def _createGaAlleleLocation(self):
        """
        Convenience method to set the common fields in a AlleleLocation
        object.
        """
        ret = protocol.AlleleLocation()
        return ret

    def toProtocolElement(self):
        """
        Converts this VariantAnnotationSet into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.VariantAnnotationSet()
        protocolElement.id = self.getId()
        protocolElement.variant_set_id = self._variantSet.getId()
        protocolElement.name = self.getLocalId()
        protocolElement.analysis.CopyFrom(self.getAnalysis())
        self.serializeAttributes(protocolElement)
        return protocolElement

    def getTranscriptEffectId(self, gaTranscriptEffect):
        effs = [eff.term for eff in gaTranscriptEffect.effects]
        return hashlib.md5(
            "{}\t{}\t{}\t{}".format(
                gaTranscriptEffect.alternate_bases,
                gaTranscriptEffect.feature_id,
                effs, gaTranscriptEffect.hgvs_annotation)
            ).hexdigest()

    def hashVariantAnnotation(cls, gaVariant, gaVariantAnnotation):
        """
        Produces an MD5 hash of the gaVariant and gaVariantAnnotation objects
        """
        treffs = [treff.id for treff in gaVariantAnnotation.transcript_effects]
        return hashlib.md5(
            "{}\t{}\t{}\t".format(
                gaVariant.reference_bases, tuple(gaVariant.alternate_bases),
                treffs)
            ).hexdigest()

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
            self.getCompoundId(), gaVariant.reference_name,
            str(gaVariant.start), md5)
        return str(compoundId)


class SimulatedVariantAnnotationSet(AbstractVariantAnnotationSet):
    """
    A variant annotation set that doesn't derive from a data store.
    Used mostly for testing.
    """
    def __init__(self, variantSet, localId, randomSeed):
        super(SimulatedVariantAnnotationSet, self).__init__(
            variantSet, localId)
        self._randomSeed = randomSeed
        self._analysis = self._createAnalysis()

    def _createAnalysis(self):
        analysis = protocol.Analysis()
        analysis.created = datetime.datetime.now().isoformat() + "Z"
        analysis.updated = datetime.datetime.now().isoformat() + "Z"
        analysis.software.append("software")
        analysis.name = "name"
        analysis.description = "description"
        analysis.id = str(datamodel.VariantAnnotationSetAnalysisCompoundId(
            self._compoundId, "analysis"))
        return analysis

    def getVariantAnnotation(self, variant, randomNumberGenerator):
        ann = self.generateVariantAnnotation(variant, randomNumberGenerator)
        return ann

    def getVariantAnnotations(self, referenceName, start, end):
        for variant in self._variantSet.getVariants(referenceName, start, end):
            yield variant, self.generateVariantAnnotation(variant)

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
        ann.variant_annotation_set_id = str(self.getCompoundId())
        ann.variant_id = variant.id
        ann.created = datetime.datetime.now().isoformat() + "Z"
        # make a transcript effect for each alternate base element
        # multiplied by a random integer (1,5)
        for base in variant.alternate_bases:
            ann.transcript_effects.add().CopyFrom(
                self.generateTranscriptEffect(
                    variant, ann, base, randomNumberGenerator))
        ann.id = self.getVariantAnnotationId(variant, ann)
        return ann

    def _addTranscriptEffectLocations(self, effect, ann, variant):
        # TODO Make these valid HGVS values
        effect.hgvs_annotation.genomic = str(variant.start)
        effect.hgvs_annotation.transcript = str(variant.start)
        effect.hgvs_annotation.protein = str(variant.start)
        effect.protein_location.start = variant.start
        effect.cds_location.start = variant.start
        effect.cdna_location.start = variant.start
        return effect

    def _getRandomOntologyTerm(self, randomNumberGenerator):
        # TODO more mock options from simulated seqOnt?
        ontologyTuples = [
            ("intron_variant", "SO:0001627", "LOW"),
            ("exon_variant", "SO:0001791", "LOW"),
            ("misense_variant", "SO:0001583", "MODERATE"),
            ("stop_gained", "SO:0001587", "HIGH")]
        term = protocol.OntologyTerm()
        picked = randomNumberGenerator.choice(ontologyTuples)
        term.term, term.term_id, impact = picked
        return term, impact

    def _addTranscriptEffectOntologyTerm(self, effect, randomNumberGenerator):
        term, impact = self._getRandomOntologyTerm(randomNumberGenerator)
        effect.effects.add().CopyFrom(term)
        protocol.setAttribute(effect.attributes.attr["impact"].values, impact)
        return effect

    def _generateAnalysisResult(self, effect, ann, randomNumberGenerator):
        # TODO make these sensible
        analysisResult = protocol.AnalysisResult()
        analysisResult.analysis_id = "analysisId"
        analysisResult.result = "result string"
        analysisResult.score = randomNumberGenerator.randint(0, 100)
        return analysisResult

    def _addAnalysisResult(self, effect, ann, randomNumberGenerator):
        effect.analysis_results.add().CopyFrom(
            self._generateAnalysisResult(
                effect, ann, randomNumberGenerator))
        return effect

    def generateTranscriptEffect(
            self, variant, ann, alts, randomNumberGenerator):
        effect = self._createGaTranscriptEffect()
        effect.alternate_bases = alts
        # TODO how to make these featureIds sensical?
        effect.feature_id = "E4TB33F"
        effect = self._addTranscriptEffectLocations(effect, ann, variant)
        effect = self._addTranscriptEffectOntologyTerm(
            effect, randomNumberGenerator)
        effect = self._addTranscriptEffectOntologyTerm(
            effect, randomNumberGenerator)
        effect.id = self.getTranscriptEffectId(effect)
        effect = self._addAnalysisResult(effect, ann, randomNumberGenerator)
        return effect


class HtslibVariantAnnotationSet(AbstractVariantAnnotationSet):
    """
    Class representing a single variant annotation derived from an
    annotated variant set.
    """
    CSQ_FIELDS = ("alternate_bases", "gene", "feature_id",
                  "featureType", "effects", "cdnaPos", "cdsPos",
                  "protPos", "aminos", "codons", "existingVar",
                  "distance", "strand", "sift", "polyPhen",
                  "motifName", "motifPos", "highInfPos",
                  "motifScoreChange")
    VEP_FIELDS = ("alternate_bases", "effects", "impact", "symbol",
                  "geneName", "featureType", "feature_id",
                  "trBiotype", "exon", "intron",
                  "hgvs_annotation.transcript",
                  "hgvs_annotation.protein", "cdnaPos", "cdsPos",
                  "protPos", "aminos", "codons", "existingVar",
                  "distance", "strand", "symbolSource", "hgncId",
                  "hgvsOffset")
    SNPEFF_FIELDS = ("alternate_bases", "effects", "impact",
                     "geneName", "geneId", "featureType",
                     "feature_id", "trBiotype", "rank",
                     "hgvs_annotation.transcript",
                     "hgvs_annotation.protein", "cdnaPos", "cdsPos",
                     "protPos", "distance", "errsWarns")
    EXCLUDED_FIELDS = ("effects", "geneName", "gene", "geneId", "featureType",
                       "trBiotype", "rank", "cdnaPos", "cdsPos",
                       "protPos", "distance", "exon", "intron",
                       "aminos", "codons", "existingVar", "distance",
                       "strand", "symbolSource", "hgncId",
                       "hgvsOffset")

    def __init__(self, variantSet, localId):
        super(HtslibVariantAnnotationSet, self).__init__(variantSet, localId)

    def populateFromFile(self, varFile, annotationType):
        self._annotationType = annotationType
        self._analysis = self._getAnnotationAnalysis(varFile)
        self._creationTime = self._analysis.created
        self._updatedTime = datetime.datetime.now().isoformat() + "Z"

    def populateFromRow(self, annotationSetRecord):
        """
        Populates this VariantAnnotationSet from the specified DB row.
        """
        self._annotationType = annotationSetRecord.annotationtype
        self._analysis = protocol.fromJson(
            annotationSetRecord.analysis, protocol.Analysis)
        self._creationTime = annotationSetRecord.created
        self._updatedTime = annotationSetRecord.updated
        self.setAttributesJson(annotationSetRecord.attributes)

    def getAnnotationType(self):
        """
        Returns the type of variant annotations, allowing us to determine
        how to interpret the annotations within the VCF file.
        """
        return self._annotationType

    def _getAnnotationAnalysis(self, varFile):
        """
        Assembles metadata within the VCF header into a GA4GH Analysis object.

        :return: protocol.Analysis
        """
        header = varFile.header
        analysis = protocol.Analysis()
        formats = header.formats.items()
        infos = header.info.items()
        filters = header.filters.items()
        for prefix, content in [("FORMAT", formats), ("INFO", infos),
                                ("FILTER", filters)]:
            for contentKey, value in content:
                key = "{0}.{1}".format(prefix, value.name)
                if key not in analysis.attributes.attr:
                    analysis.attributes.attr[key].Clear()
                if value.description is not None:
                    analysis.attributes.attr[
                        key].values.add().string_value = value.description
        analysis.created = self._creationTime
        analysis.updated = self._updatedTime
        for r in header.records:
            # Don't add a key to info if there's nothing in the value
            if r.value is not None:
                if r.key not in analysis.attributes.attr:
                    analysis.attributes.attr[r.key].Clear()
                analysis.attributes.attr[r.key] \
                    .values.add().string_value = str(r.value)
            if r.key == "created" or r.key == "fileDate":
                # TODO handle more date formats
                try:
                    if '-' in r.value:
                        fmtStr = "%Y-%m-%d"
                    else:
                        fmtStr = "%Y%m%d"
                    analysis.created = datetime.datetime.strptime(
                        r.value, fmtStr).isoformat() + "Z"
                except ValueError:
                    # is there a logger we should tell?
                    # print("INFO: Could not parse variant annotation time")
                    pass  # analysis.create_date_time remains datetime.now()
            if r.key == "software":
                analysis.software.append(r.value)
            if r.key == "name":
                analysis.name = r.value
            if r.key == "description":
                analysis.description = r.value
        analysis.id = str(datamodel.VariantAnnotationSetAnalysisCompoundId(
            self._compoundId, "analysis"))
        return analysis

    def getVariantAnnotations(self, referenceName, startPosition, endPosition):
        """
        Generator for iterating through variant annotations in this
        variant annotation set.
        :param referenceName:
        :param startPosition:
        :param endPosition:
        :return: generator of protocol.VariantAnnotation
        """
        variantIter = self._variantSet.getPysamVariants(
            referenceName, startPosition, endPosition)
        for record in variantIter:
            yield self.convertVariantAnnotation(record)

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
                allLoc.reference_sequence = match.group(2)
                allLoc.alternate_sequence = match.group(3)
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
            allLoc.reference_sequence = match.group(1)
            allLoc.start = int(match.group(2)) - 1
            allLoc.alternate_sequence = match.group(3)
            return allLoc
        return None

    def addCDSLocation(self, effect, cdnaPos):
        hgvsC = effect.hgvs_annotation.transcript
        allele_location = None
        if not isUnspecified(hgvsC):
            allele_location = self.convertLocationHgvsC(hgvsC)
            if allele_location:
                effect.cds_location.CopyFrom(self.convertLocationHgvsC(hgvsC))
        if allele_location is None and self.convertLocation(cdnaPos):
                effect.cds_location.CopyFrom(self.convertLocation(cdnaPos))
        else:
            # These are not stored in the VCF
            effect.cds_location.alternate_sequence = ""
            effect.cds_location.reference_sequence = ""

    def addProteinLocation(self, effect, protPos):
        hgvsP = effect.hgvs_annotation.protein
        protein_location = None
        if not isUnspecified(hgvsP):
            protein_location = self.convertLocationHgvsP(hgvsP)
            if protein_location:
                effect.protein_location.CopyFrom(
                    self.convertLocationHgvsP(hgvsP))
        if protein_location is None and self.convertLocation(protPos):
            effect.protein_location.CopyFrom(self.convertLocation(protPos))

    def addCDNALocation(self, effect, cdnaPos):
        hgvsC = effect.hgvs_annotation.transcript
        if self.convertLocation(cdnaPos):
            effect.cdna_location.CopyFrom(self.convertLocation(cdnaPos))
        if self.convertLocationHgvsC(hgvsC):
            effect.cdna_location.alternate_sequence = \
                self.convertLocationHgvsC(hgvsC).alternate_sequence
            effect.cdna_location.reference_sequence = \
                self.convertLocationHgvsC(hgvsC).reference_sequence

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

    def convertTranscriptEffect(self, annStr, hgvsG):
        """
        Takes the ANN string of a SnpEff generated VCF, splits it
        and returns a populated GA4GH transcript effect object.
        :param annStr: String
        :param hgvsG: String
        :return: effect protocol.TranscriptEffect()
        """
        effect = self._createGaTranscriptEffect()
        effect.hgvs_annotation.CopyFrom(protocol.HGVSAnnotation())
        annDict = dict()
        if self._annotationType == ANNOTATIONS_SNPEFF:
            annDict = dict(zip(self. SNPEFF_FIELDS, annStr.split("|")))
        elif self._annotationType == ANNOTATIONS_VEP_V82:
            annDict = dict(zip(self.VEP_FIELDS, annStr.split("|")))
        else:
            annDict = dict(zip(self.CSQ_FIELDS, annStr.split("|")))
        annDict["hgvs_annotation.genomic"] = hgvsG if hgvsG else u''
        for key, val in annDict.items():
            try:
                protocol.deepSetAttr(effect, key, val)
            except AttributeError:
                if val and key not in self.EXCLUDED_FIELDS:
                    protocol.setAttribute(
                        effect.attributes.attr[key].values, val)
        effect.effects.extend(self.convertSeqOntology(annDict.get('effects')))
        self.addLocations(
            effect, annDict.get('protPos'), annDict.get('cdnaPos'))
        effect.id = self.getTranscriptEffectId(effect)
        return effect

    def convertSeqOntology(self, seqOntStr):
        """
        Splits a string of sequence ontology effects and creates
        an ontology term record for each, which are built into
        an array of return soTerms.
        :param seqOntStr:
        :return: [protocol.OntologyTerm]
        """
        return [
            self._ontology.getGaTermByName(soName)
            for soName in seqOntStr.split('&')]

    def convertVariantAnnotation(self, record):
        """
        Converts the specfied pysam variant record into a GA4GH variant
        annotation object using the specified function to convert the
        transcripts.
        """
        variant = self._variantSet.convertVariant(record, [])
        annotation = self._createGaVariantAnnotation()
        annotation.variant_id = variant.id
        gDots = record.info.get(b'HGVS.g')
        # Convert annotations from INFO field into TranscriptEffect
        transcriptEffects = []
        annotations = record.info.get(b'ANN') or record.info.get(b'CSQ')
        for i, ann in enumerate(annotations):
            hgvsG = gDots[i % len(variant.alternate_bases)] if gDots else None
            transcriptEffects.append(self.convertTranscriptEffect(ann, hgvsG))
        annotation.transcript_effects.extend(transcriptEffects)
        annotation.id = self.getVariantAnnotationId(variant, annotation)
        return variant, annotation
