"""
Module responsible for translating read data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob

import pysam

import ga4gh.protocol as protocol


class ReadGroupSet(object):
    """
    Class representing a logical collection ReadGroups.
    """
    def __init__(self, id_, dataDir):
        self._id = id_
        self._dataDir = dataDir
        self._readGroups = []
        for fileType in ["sam", "bam"]:
            pattern = "*.{}".format(fileType)
            for path in glob.glob(os.path.join(self._dataDir, pattern)):
                filename = os.path.split(path)[1]
                localId = filename.split(".")[0]
                readGroupId = "{}:{}".format(self._id, localId)
                readGroup = ReadGroup(readGroupId, path)
                self._readGroups.append(readGroup)

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReadGroupSet.
        """
        readGroupSet = protocol.GAReadGroupSet()
        readGroupSet.id = self._id
        readGroupSet.readGroups = [
            readGroup.toProtocolElement() for readGroup in self._readGroups]
        # TODO fill out details
        return readGroupSet


class ReadGroup(object):
    """
    Class representing a ReadGroup. A ReadGroup is all the data that's
    processed the same way by the sequencer.  There are typically 1-10
    ReadGroups in a ReadGroupSet.
    """
    def __init__(self, id_, dataFile):
        self._id = id_
        self._samFile = pysam.AlignmentFile(dataFile)

    def toProtocolElement(self):
        """
        Returns the GA4GH protocol representation of this ReadGroup.
        """
        readGroup = protocol.GAReadGroup()
        readGroup.id = self._id
        # TODO fill out details
        return readGroup
