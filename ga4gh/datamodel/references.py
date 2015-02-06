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
    Class representing a reference set in the GA4GH datamodel. A reference
    set consists of a directory containing a set of indexed FASTA files.
    """
    def __init__(self, id_, dataDir):
        self._id = id_
        self._dataDir = dataDir
        self._referenceIdMap = {}
        self._referenceIds = []
        for relPath in glob.glob(os.path.join(dataDir, "*.fa.gz")):
            base = os.path.basename(relPath).split(".")[0]
            referenceId = "{}:{}".format(self._id, base)
            reference = Reference(referenceId, relPath)
            self._referenceIdMap[referenceId] = reference
            self._referenceIds.append(referenceId)

    def toProtocolElement(self):
        """
        Returns a representation of this ReferenceSet as a GA4GH protocol
        element.
        """
        protocolElement = protocol.GAReferenceSet()
        protocolElement.id = self._id
        protocolElement.referenceIds = list(self._referenceIdMap.keys())
        protocolElement.md5checksum = "TODO"
        # TODO fill in the remaining fields.
        return protocolElement


    def getReferences(self, startIndex, md5checksums, accessions):
        """
        Returns an iterator over the references in this reference set.

        TODO document
        """
        # TODO paging is broken here; FIXME!
        for offset, id_ in enumerate(self._referenceIds[startIndex:-1]):
            reference = self._referenceIdMap[id_]
            yield reference.toProtocolElement(), startIndex + offset + 1
        reference = self._referenceIdMap[-1]
        yield reference.toProtocolElement(), None


class Reference(object):
    """
    Class reprenting a single Reference in the GA4GH data model. This
    corresponds to a single assembled chromosome.

    TODO document this more accurately based on the official docs.
    """
    def __init__(self, id_, dataFile):
        self._id = id_
        self._fastaFile = pysam.FastaFile(dataFile)

    def toProtocolElement(self):
        """
        Returns a representation of this Reference as a GA4GH protocol
        element.
        """
        protocolElement = protocol.GAReference()
        protocolElement.id = self._id
        # TODO fill in the remaining fields.
        return protocolElement
