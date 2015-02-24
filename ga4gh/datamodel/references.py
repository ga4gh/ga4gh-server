"""
Module responsible for translating reference sequence data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob

import pysam

import ga4gh.protocol as protocol


class ReferenceSet(object):
    """
    Class representing ReferenceSets. A ReferenceSet is a set of
    References which typically comprise a reference assembly, such as
    GRCh38.
    """
    def __init__(self, id_, dataDir):
        self._id = id_
        self._dataDir = dataDir
        self._referenceIdMap = {}
        # TODO get metadata from a file within dataDir? How else will we
        # fill in the fields like ncbiTaxonId etc?
        for relativePath in glob.glob(os.path.join(self._dataDir, "*.fa.gz")):
            filename = os.path.split(relativePath)[1]
            localId = filename.split(".")[0]
            referenceId = "{}:{}".format(self._id, localId)
            reference = Reference(referenceId, relativePath)
            self._referenceIdMap[referenceId] = reference
        self._referenceIds = sorted(self._referenceIdMap.keys())

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReferenceSet.
        """
        referenceSet = protocol.GAReferenceSet()
        referenceSet.id = self._id
        # TODO fill out details
        return referenceSet


class Reference(object):
    """
    Class representing References. A Reference is a canonical
    assembled contig, intended to act as a reference coordinate space
    for other genomic annotations. A single Reference might represent
    the human chromosome 1, for instance.
    """
    def __init__(self, id_, dataFile):
        self._id = id_
        self._fastaFile = pysam.FastaFile(dataFile)

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this Reference.
        """
        reference = protocol.GAReference()
        reference.id = self._id
        # TODO fill out details
        return reference
