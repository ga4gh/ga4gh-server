"""
Module responsible for translating reference sequence data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import hashlib
import json
import random

import pysam

import ga4gh.server.datamodel as datamodel
import ga4gh.server.exceptions as exceptions

import ga4gh.schemas.pb as pb
import ga4gh.schemas.protocol as protocol


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
        self._species = None
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

    def setDescription(self, description):
        """
        Sets the description to the specified value.
        """
        self._description = description

    def setSpeciesFromJson(self, speciesJson):
        """
        Sets the species, an OntologyTerm, to the specified value, given as
        a JSON string.

        See the documentation for details of this field.
        """
        try:
            parsed = protocol.fromJson(speciesJson, protocol.OntologyTerm)
        except:
            raise exceptions.InvalidJsonException(speciesJson)
        self._species = protocol.toJsonDict(parsed)

    def setIsDerived(self, isDerived):
        """
        Sets the isDerived attribute
        """
        self._isDerived = isDerived

    def setAssemblyId(self, assemblyId):
        """
        Sets the assemblyId attribute
        """
        self._assemblyId = assemblyId

    def setSourceAccessions(self, sourceAccessions):
        """
        Sets the sourceAccessions attribute
        """
        self._sourceAccessions = sourceAccessions

    def setSourceUri(self, sourceUri):
        """
        Sets the sourceUri attribute
        """
        self._sourceUri = sourceUri

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

    def getSpecies(self):
        """
        Returns the species for this reference set. This is the
        ontology term with data from
        www.obofoundry.org/ontology/ncbitaxon.html
        (e.g. 9606 for human)
        Note that contained `Reference`s may specify a different
        species, as assemblies may contain reference sequences
        which do not belong to the modeled species, e.g. EBV in a
        human reference genome.
        """
        if self._species is not {}:
            return self._species
        else:
            return None

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReferenceSet.
        """
        ret = protocol.ReferenceSet()
        ret.assembly_id = pb.string(self.getAssemblyId())
        ret.description = pb.string(self.getDescription())
        ret.id = self.getId()
        ret.is_derived = self.getIsDerived()
        ret.md5checksum = self.getMd5Checksum()
        if self.getSpecies():
            term = protocol.fromJson(
                json.dumps(self.getSpecies()), protocol.OntologyTerm)
            ret.species.term_id = term.term_id
            ret.species.term = term.term
        ret.source_accessions.extend(self.getSourceAccessions())
        ret.source_uri = pb.string(self.getSourceUri())
        ret.name = self.getLocalId()
        self.serializeAttributes(ret)
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
        self._sourceUri = ""
        self._sourceAccessions = []
        self._isDerived = False
        self._sourceDivergence = pb.DEFAULT_INT
        self._species = None

    def setMd5checksum(self, md5checksum):
        """
        Sets the md5checksum to the specified value.
        """
        self._md5checksum = md5checksum

    def setSpeciesFromJson(self, speciesJson):
        """
        Sets the species, an OntologyTerm, to the specified value, given as
        a JSON string.

        See the documentation for details of this field.
        """
        try:
            parsed = protocol.fromJson(speciesJson, protocol.OntologyTerm)
        except:
            raise exceptions.InvalidJsonException(speciesJson)
        self._species = protocol.toJsonDict(parsed)

    def setSourceAccessions(self, sourceAccessions):
        """
        Sets the sourceAccessions to the specified list. See the documentation
        for getSourceAccessions for details on this field.
        """
        self._sourceAccessions = sourceAccessions

    def setLength(self, length):
        """
        Sets the length of this referefnce to the specified value.
        """
        self._length = length

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

    def getSpecies(self):
        """
        Returns the species for this reference set. This is the
        ontology term with data from
        www.obofoundry.org/ontology/ncbitaxon.html
        (e.g. 9606 for human)
        Note that contained `Reference`s may specify a different
        species, as assemblies may contain reference sequences
        which do not belong to the modeled species, e.g. EBV in a
        human reference genome.
        """
        if self._species is not {}:
            return self._species
        else:
            return None

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
        reference.is_derived = self.getIsDerived()
        reference.length = self.getLength()
        reference.md5checksum = self.getMd5Checksum()
        reference.name = self.getName()
        if self.getSpecies():
            term = protocol.fromJson(
                json.dumps(self.getSpecies()), protocol.OntologyTerm)
            reference.species.term_id = term.term_id
            reference.species.term = term.term
        reference.source_accessions.extend(self.getSourceAccessions())
        reference.source_divergence = pb.int(self.getSourceDivergence())
        reference.source_uri = self.getSourceUri()
        self.serializeAttributes(reference)
        return reference

    def checkQueryRange(self, start, end):
        """
        Checks to ensure that the query range is valid within this reference.
        If not, raise ReferenceRangeErrorException.
        """
        condition = (
            (start < 0 or end > self.getLength()) or
            start > end or start == end)
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
        self._species = json.loads(
                    '{"term": "Homo sapiens", "termId": "9606"}')
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
        self._species = json.loads(
                            '{"term": "Homo sapiens", "termId": "9606"}')
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
    def __init__(self, localId):
        super(HtslibReferenceSet, self).__init__(localId)
        self._dataUrl = None

    def populateFromFile(self, dataUrl):
        """
        Populates the instance variables of this ReferencSet from the
        data URL.
        """
        self._dataUrl = dataUrl
        fastaFile = self.getFastaFile()
        for referenceName in fastaFile.references:
            reference = HtslibReference(self, referenceName)
            # TODO break this up into chunks and calculate the MD5
            # in bits (say, 64K chunks?)
            bases = fastaFile.fetch(referenceName)
            md5checksum = hashlib.md5(bases).hexdigest()
            reference.setMd5checksum(md5checksum)
            reference.setLength(len(bases))
            self.addReference(reference)

    def populateFromRow(self, referenceSetRecord):
        """
        Populates this reference set from the values in the specified DB
        row.
        """
        self._dataUrl = referenceSetRecord.dataurl
        self._description = referenceSetRecord.description
        self._assemblyId = referenceSetRecord.assemblyid
        self._isDerived = bool(referenceSetRecord.isderived)
        self._md5checksum = referenceSetRecord.md5checksum
        species = referenceSetRecord.species
        if species is not None and species != 'null':
            self.setSpeciesFromJson(species)
        self._sourceAccessions = json.loads(
            referenceSetRecord.sourceaccessions)
        self._sourceUri = referenceSetRecord.sourceuri

    def getDataUrl(self):
        """
        Returns the path of the the data URL for this ReferenceSet.
        """
        return self._dataUrl

    def openFile(self, dataFile):
        return pysam.FastaFile(dataFile)

    def getFastaFile(self):
        """
        Returns a reference to the Fasta file instance used to read the
        data in this reference set.
        """
        return self.getFileHandle(self._dataUrl)


class HtslibReference(datamodel.PysamDatamodelMixin, AbstractReference):
    """
    A reference based on data stored in a file on the file system
    """
    def __init__(self, parentContainer, localId):
        super(HtslibReference, self).__init__(parentContainer, localId)

    def populateFromRow(self, referenceRecord):
        """
        Populates this reference from the values in the specified DB row.
        """
        self._length = referenceRecord.length
        self._isDerived = bool(referenceRecord.isderived)
        self._md5checksum = referenceRecord.md5checksum
        species = referenceRecord.species
        if species is not None and species != 'null':
            self.setSpeciesFromJson(species)
        self._sourceAccessions = json.loads(referenceRecord.sourceaccessions)
        self._sourceDivergence = referenceRecord.sourcedivergence
        self._sourceUri = referenceRecord.sourceuri

    def getBases(self, start, end):
        self.checkQueryRange(start, end)
        fastaFile = self._parentContainer.getFastaFile()
        localId = self.getLocalId().encode()
        # TODO we should have some error checking here...
        bases = fastaFile.fetch(localId, start, end)
        return bases
