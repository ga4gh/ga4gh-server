"""
Module responsible for translating reference sequence data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import hashlib
import json
import os
import random

import pysam

import ga4gh.datamodel as datamodel
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


DEFAULT_REFERENCESET_NAME = "Default"
"""
This is the name used for any reference set referred to in a BAM
file that does not provide the 'AS' tag in the @SQ header.
"""


class AbstractReferenceSet(datamodel.DatamodelObject):
    """
    Class representing ReferenceSets. A ReferenceSet is a set of
    References which typically comprise a reference assembly, such as
    GRCh38.
    """
    compoundIdClass = datamodel.ReferenceSetCompoundId

    def __init__(self, localId):
        super(AbstractReferenceSet, self).__init__(None, localId)
        self._referenceIdMap = {}
        self._referenceNameMap = {}
        self._referenceIds = []
        self._assemblyId = None
        self._description = None
        self._isDerived = False
        self._ncbiTaxonId = None
        self._sourceAccessions = []
        self._sourceUri = None

    def addReference(self, reference):
        """
        Adds the specified reference to this ReferenceSet.
        """
        id_ = reference.getId()
        self._referenceIdMap[id_] = reference
        self._referenceNameMap[reference.getLocalId()] = reference
        self._referenceIds.append(id_)

    def getReferences(self):
        """
        Returns the References in this ReferenceSet.
        """
        return [self._referenceIdMap[id_] for id_ in self._referenceIds]

    def getNumReferences(self):
        """
        Returns the number of references in this ReferenceSet.
        """
        return len(self._referenceIds)

    def getReferenceByIndex(self, index):
        """
        Returns the reference at the specified index in this ReferenceSet.
        """
        return self._referenceIdMap[self._referenceIds[index]]

    def getReferenceByName(self, name):
        """
        Returns the reference with the specified name.
        """
        if name not in self._referenceNameMap:
            raise exceptions.ReferenceNameNotFoundException(name)
        return self._referenceNameMap[name]

    def getReference(self, id_):
        """
        Returns the Reference with the specified ID or raises a
        ReferenceNotFoundException if it does not exist.
        """
        if id_ not in self._referenceIdMap:
            raise exceptions.ReferenceNotFoundException(id_)
        return self._referenceIdMap[id_]

    def getMd5Checksum(self):
        """
        Returns the MD5 checksum for this reference set. This checksum is
        calculated by making a list of `Reference.md5checksum` for all
        `Reference`s in this set. We then sort this list, and take the
        MD5 hash of all the strings concatenated together.
        """
        references = sorted(
            self.getReferences(),
            key=lambda ref: ref.getMd5Checksum())
        checksums = ''.join([ref.getMd5Checksum() for ref in references])
        md5checksum = hashlib.md5(checksums).hexdigest()
        return md5checksum

    def getAssemblyId(self):
        """
        Returns the assembly ID for this reference set.
        This is the public id of this reference set, such as `GRCh37`
        """
        return self._assemblyId

    def getDescription(self):
        """
        Returns the free text description of this reference set.
        """
        return self._description

    def getIsDerived(self):
        """
        Returns True if this ReferenceSet is derived. A ReferenceSet
        may be derived from a source if it contains additional sequences,
        or some of the sequences within it are derived.
        """
        return self._isDerived

    def getSourceAccessions(self):
        """
        Returns the list of source accession strings. These are all known
        corresponding accession IDs in INSDC (GenBank/ENA/DDBJ) ideally
        with a version number, e.g. `NC_000001.11`.
        """
        return self._sourceAccessions

    def getSourceUri(self):
        """
        Returns the sourceURI for this ReferenceSet.
        """
        return self._sourceUri

    def getNcbiTaxonId(self):
        """
        Returns the NCBI Taxon ID for this reference set. This is the
        ID from http://www.ncbi.nlm.nih.gov/taxonomy (e.g. 9606->human)
        indicating the species which this assembly is intended to model.
        Note that contained `Reference`s may specify a different
        `ncbiTaxonId`, as assemblies may contain reference sequences
        which do not belong to the modeled species, e.g.  EBV in a
        human reference genome.
        """
        return self._ncbiTaxonId

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReferenceSet.
        """
        ret = protocol.ReferenceSet()
        ret.assemblyId = self.getAssemblyId()
        ret.description = self.getDescription()
        ret.id = self.getId()
        ret.isDerived = self.getIsDerived()
        ret.md5checksum = self.getMd5Checksum()
        ret.ncbiTaxonId = self.getNcbiTaxonId()
        ret.referenceIds = self._referenceIds
        ret.sourceAccessions = self.getSourceAccessions()
        ret.sourceURI = self.getSourceUri()
        ret.name = self.getLocalId()
        return ret


class AbstractReference(datamodel.DatamodelObject):
    """
    Class representing References. A Reference is a canonical
    assembled contig, intended to act as a reference coordinate space
    for other genomic annotations. A single Reference might represent
    the human chromosome 1, for instance.
    """
    compoundIdClass = datamodel.ReferenceCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractReference, self).__init__(parentContainer, localId)
        self._length = -1
        self._md5checksum = ""
        self._sourceUri = None
        self._sourceAccessions = []
        self._isDerived = False
        self._sourceDivergence = None
        self._ncbiTaxonId = None

    def getLength(self):
        """
        Returns the length of this reference's sequence string.
        """
        return self._length

    def getName(self):
        """
        Returns the name of this reference, e.g., '22'.
        """
        return self.getLocalId()

    def getIsDerived(self):
        """
        Returns True if this Reference is derived. A sequence X is said to be
        derived from source sequence Y, if X and Y are of the same length and
        the per-base sequence divergence at A/C/G/T bases is sufficiently
        small. Two sequences derived from the same official sequence share the
        same coordinates and annotations, and can be replaced with the official
        sequence for certain use cases.
        """
        return self._isDerived

    def getSourceDivergence(self):
        """
        Returns the source divergence for this reference.  The sourceDivergence
        is the fraction of non-indel bases that do not match the
        reference this record was derived from.
        """
        return self._sourceDivergence

    def getSourceAccessions(self):
        """
        Returns the list of source accession strings. These are all known
        corresponding accession IDs in INSDC (GenBank/ENA/DDBJ) ideally
        with a version number, e.g. `NC_000001.11`.
        """
        return self._sourceAccessions

    def getSourceUri(self):
        """
        The URI from which the sequence was obtained. Specifies a FASTA format
        file/string with one name, sequence pair.
        """
        return self._sourceUri

    def getNcbiTaxonId(self):
        """
        Returns the NCBI Taxon ID for this reference. This is the
        ID from http://www.ncbi.nlm.nih.gov/taxonomy (e.g. 9606->human)
        indicating the species which this assembly is intended to model.
        Note that contained `Reference`s may specify a different
        `ncbiTaxonId`, as assemblies may contain reference sequences
        which do not belong to the modeled species, e.g.  EBV in a
        human reference genome.
        """
        return self._ncbiTaxonId

    def getMd5Checksum(self):
        """
        Returns the MD5 checksum uniquely representing this `Reference` as a
        lower-case hexadecimal string, calculated as the MD5 of the upper-case
        sequence excluding all whitespace characters.
        """
        return self._md5checksum

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this Reference.
        """
        reference = protocol.Reference()
        reference.id = self.getId()
        reference.isDerived = self.getIsDerived()
        reference.length = self.getLength()
        reference.md5checksum = self.getMd5Checksum()
        reference.name = self.getName()
        reference.ncbiTaxonId = self.getNcbiTaxonId()
        reference.sourceAccessions = self.getSourceAccessions()
        reference.sourceDivergence = self.getSourceDivergence()
        reference.sourceURI = self.getSourceUri()
        return reference

    def checkQueryRange(self, start, end):
        """
        Checks to ensure that the query range is valid within this reference.
        If not, raise ReferenceRangeErrorException.
        """
        condition = (
            (start < 0 or end > self.getLength()) or
            start > end)
        if condition:
            raise exceptions.ReferenceRangeErrorException(
                self.getId(), start, end)

    def getBases(self, start, end):
        """
        Returns the string representing the bases of this reference from
        start (inclusive) to end (exclusive).
        """
        raise NotImplemented()

##################################################################
#
# Simulated references
#
##################################################################


class SimulatedReferenceSet(AbstractReferenceSet):
    """
    A simulated referenceSet
    """
    def __init__(self, localId, randomSeed=0, numReferences=1):
        super(SimulatedReferenceSet, self).__init__(localId)
        self._randomSeed = randomSeed
        self._randomGenerator = random.Random()
        self._randomGenerator.seed(self._randomSeed)
        self._description = "Simulated reference set"
        self._assemblyId = str(random.randint(0, 2**32))
        self._isDerived = bool(random.randint(0, 1))
        self._ncbiTaxonId = random.randint(0, 2**16)
        self._sourceAccessions = []
        for i in range(random.randint(1, 3)):
                self._sourceAccessions.append("sim_accession_{}".format(
                    random.randint(1, 2**32)))
        self._sourceUri = "http://example.com/reference.fa"
        for i in range(numReferences):
            referenceSeed = self._randomGenerator.getrandbits(32)
            referenceLocalId = "srs{}".format(i)
            reference = SimulatedReference(
                self, referenceLocalId, referenceSeed)
            self.addReference(reference)


class SimulatedReference(AbstractReference):
    """
    A simulated reference. Stores a random sequence of a given length, and
    generates remaining attributes randomly.
    """

    def __init__(self, parentContainer, localId, randomSeed=0, length=200):
        super(SimulatedReference, self).__init__(parentContainer, localId)
        rng = random.Random()
        rng.seed(randomSeed)
        self._length = length
        bases = [rng.choice('ACGT') for _ in range(self._length)]
        self._bases = ''.join(bases)
        self._md5checksum = hashlib.md5(self._bases).hexdigest()
        self._isDerived = bool(rng.randint(0, 1))
        self._sourceDivergence = 0
        if self._isDerived:
            self._sourceDivergence = rng.uniform(0, 0.1)
        self._ncbiTaxonId = random.randint(0, 2**16)
        self._sourceAccessions = []
        for i in range(random.randint(1, 3)):
                self._sourceAccessions.append("sim_accession_{}".format(
                    random.randint(1, 2**32)))
        self._sourceUri = "http://example.com/reference.fa"

    def getBases(self, start, end):
        self.checkQueryRange(start, end)
        return self._bases[start:end]

##################################################################
#
# References based on htslib's FASTA file handling.
#
##################################################################


class HtslibReferenceSet(datamodel.PysamDatamodelMixin, AbstractReferenceSet):
    """
    A referenceSet based on data on a file system
    """
    def __init__(self, localId, dataDir, backend):
        super(HtslibReferenceSet, self).__init__(localId)
        self._dataDir = dataDir
        self._setMetadata()
        self._scanDataFiles(dataDir, ["*.fa.gz"])

    def _setMetadata(self):
        metadataFileName = '{}.json'.format(self._dataDir)
        with open(metadataFileName) as metadataFile:
            metadata = json.load(metadataFile)
            try:
                self._assemblyId = metadata['assemblyId']
                self._description = metadata['description']
                self._isDerived = metadata['isDerived']
                self._ncbiTaxonId = metadata['ncbiTaxonId']
                self._sourceAccessions = metadata['sourceAccessions']
                self._sourceUri = metadata['sourceUri']
            except KeyError as err:
                raise exceptions.MissingReferenceSetMetadata(
                    metadataFileName, str(err))

    def _addDataFile(self, path):
        dirname, filename = os.path.split(path)
        localId = filename.split(".")[0]
        metadataFileName = os.path.join(dirname, "{}.json".format(localId))
        with open(metadataFileName) as metadataFile:
            metadata = json.load(metadataFile)
        reference = HtslibReference(self, localId, path, metadata)
        self.addReference(reference)


class HtslibReference(datamodel.PysamDatamodelMixin, AbstractReference):
    """
    A reference based on data stored in a file on the file system
    """
    def __init__(self, parentContainer, localId, dataFile, metadata):
        super(HtslibReference, self).__init__(parentContainer, localId)
        self._fastaFilePath = dataFile
        fastaFile = self.getFileHandle(dataFile)
        numReferences = len(fastaFile.references)
        if numReferences != 1:
            raise exceptions.NotExactlyOneReferenceException(
                self._fastaFilePath, numReferences)
        if fastaFile.references[0] != localId:
            raise exceptions.InconsistentReferenceNameException(
                self._fastaFilePath)
        self._length = fastaFile.lengths[0]
        try:
            self._md5checksum = metadata["md5checksum"]
            self._sourceUri = metadata["sourceUri"]
            self._ncbiTaxonId = metadata["ncbiTaxonId"]
            self._isDerived = metadata["isDerived"]
            self._sourceDivergence = metadata["sourceDivergence"]
            self._sourceAccessions = metadata["sourceAccessions"]
        except KeyError as err:
            raise exceptions.MissingReferenceMetadata(dataFile, str(err))

    def getFastaFilePath(self):
        """
        Returns the fasta file that this reference is derived from.
        """
        return self._fastaFilePath

    def openFile(self, dataFile):
        return pysam.FastaFile(dataFile)

    def getBases(self, start, end):
        self.checkQueryRange(start, end)
        fastaFile = self.getFileHandle(self._fastaFilePath)
        # TODO we should have some error checking here...
        bases = fastaFile.fetch(self.getLocalId(), start, end)
        return bases
