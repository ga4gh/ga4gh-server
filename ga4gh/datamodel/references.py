"""
Module responsible for translating reference sequence data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import random
import hashlib

import pysam

import ga4gh.datamodel as datamodel
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


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
        self._referenceIds = []

    def getReferences(self):
        """
        Returns the References in this ReferenceSet.
        """
        return [self._referenceIdMap[id_] for id_ in self._referenceIds]

    def getReferenceIdMap(self):
        return self._referenceIdMap

    def getReferenceIds(self):
        return self._referenceIds

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReferenceSet.
        """
        ret = protocol.ReferenceSet()
        ret.assemblyId = None
        ret.description = None
        ret.id = self.getId()
        ret.isDerived = False
        ret.md5checksum = self._generateMd5Checksum()
        ret.ncbiTaxonId = None
        ret.referenceIds = self._referenceIds
        ret.sourceAccessions = []
        ret.sourceURI = None
        return ret

    def _generateMd5Checksum(self):
        return "TODO"
        """
        references = sorted(
            self.getReferences(),
            key=lambda ref: ref.getMd5Checksum())
        checksums = [ref.getMd5Checksum() for ref in references]
        checksumsString = ''.join(checksums)
        md5checksum = hashlib.md5(checksumsString).hexdigest()
        return md5checksum
        """


class SimulatedReferenceSet(AbstractReferenceSet):
    """
    A simulated referenceSet
    """
    def __init__(self, localId, randomSeed=0, numReferences=1):
        super(SimulatedReferenceSet, self).__init__(localId)
        self._randomSeed = randomSeed
        self._randomGenerator = random.Random()
        self._randomGenerator.seed(self._randomSeed)
        for i in range(numReferences):
            referenceSeed = self._randomGenerator.getrandbits(32)
            referenceLocalId = "srs{}".format(i)
            reference = SimulatedReference(
                self, referenceLocalId, referenceSeed)
            self._referenceIdMap[reference.getId()] = reference
        self._referenceIds = sorted(self._referenceIdMap.keys())


class HtslibReferenceSet(datamodel.PysamDatamodelMixin, AbstractReferenceSet):
    """
    A referenceSet based on data on a file system
    """
    def __init__(self, localId, dataDir):
        super(HtslibReferenceSet, self).__init__(localId)
        self._dataDir = dataDir
        # TODO get metadata from a file within dataDir? How else will we
        # fill in the fields like ncbiTaxonId etc?
        self._scanDataFiles(dataDir, ["*.fa.gz"])
        self._referenceIds = sorted(self._referenceIdMap.keys())

    def _addDataFile(self, path):
        filename = os.path.split(path)[1]
        localId = filename.split(".")[0]
        reference = HtslibReference(self, localId, path)
        self._referenceIdMap[reference.getId()] = reference


class AbstractReference(datamodel.DatamodelObject):
    """
    Class representing References. A Reference is a canonical
    assembled contig, intended to act as a reference coordinate space
    for other genomic annotations. A single Reference might represent
    the human chromosome 1, for instance.
    """
    compoundIdClass = datamodel.ReferenceCompoundId

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this Reference.
        """
        reference = protocol.Reference()
        reference.id = self.getId()
        reference.isDerived = False
        reference.length = self.getLength()
        reference.md5checksum = self.getMd5Checksum()
        reference.name = self.getName()
        reference.ncbiTaxonId = None
        reference.sourceAccessions = []
        reference.sourceDivergence = None
        reference.sourceURI = None
        return reference

    def getMd5Checksum(self):
        """
        Returns the md5 checksum
        """
        return self._md5checksum


class SimulatedReference(AbstractReference):
    """
    A simulated reference
    """
    choices = 'AGCT'

    def __init__(self, parentContainer, name, randomSeed=0, length=200):
        super(SimulatedReference, self).__init__(parentContainer, name)
        self._randomSeed = randomSeed
        self._randomGenerator = random.Random()
        self._randomGenerator.seed(self._randomSeed)
        self._length = length
        bases = []
        for _ in range(length):
            choice = self._randomGenerator.choice(self.choices)
            bases.append(choice)
        self.bases = ''.join(bases)
        self._md5checksum = hashlib.md5(self.bases).hexdigest()

    def getBases(self, start=None, end=None):
        return self.bases[start:end]

    def getLength(self):
        return len(self.bases)

    def getName(self):
        return self.getLocalId()


class HtslibReference(datamodel.PysamDatamodelMixin, AbstractReference):
    """
    A reference based on data stored in a file on the file system
    """
    def __init__(self, parentContainer, name, dataFile):
        super(HtslibReference, self).__init__(parentContainer, name)
        self._fastaFilePath = dataFile
        fastaFile = self.openFile(dataFile)
        numReferences = len(fastaFile.references)
        if numReferences != 1:
            raise exceptions.NotExactlyOneReferenceException(
                self.getId(), numReferences)
        self._refName = fastaFile.references[0]
        self._length = fastaFile.lengths[0]
        # refData = fastaFile.fetch(self._refName)
        self._md5checksum = "TODO"  # hashlib.md5(refData).hexdigest()
        fastaFile.close()

    def openFile(self, dataFile):
        return pysam.FastaFile(dataFile)

    def getFastaFilePath(self):
        """
        Returns the file path of the fasta file
        """
        return self._fastaFilePath

    def getBases(self, start=None, end=None):
        start, end = self.sanitizeFastaFileFetch(start, end)
        bases = self.getFileHandle(self.getFastaFilePath()).fetch(
                      self._refName, start, end)
        return bases

    def getLength(self):
        return self._length

    def getName(self):
        return self._refName
