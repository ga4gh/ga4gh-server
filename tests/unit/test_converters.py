"""
Tests the converters
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import tempfile
import unittest

import pysam

import ga4gh.protocol as protocol
import ga4gh.converters as converters


class TestSamConverter(unittest.TestCase):
    """
    Base class for SAM converter tests.
    Provides common  methods.
    """
    readsFilePath = 'tests/unit/reads.dat'

    def getReads(self):
        readsFile = file(self.readsFilePath)
        lines = readsFile.readlines()
        reads = []
        for line in lines:
            read = protocol.ReadAlignment.fromJsonString(line)
            reads.append(read)
        return reads


class TestSamConverterRoundTrip(TestSamConverter):
    """
    Write a sam file and see if pysam can read it
    """
    def _testRoundTrip(self, binaryOutput):
        with tempfile.NamedTemporaryFile() as fileHandle:
            # write SAM file
            filePath = fileHandle.name
            samConverter = converters.SamConverter(
                None, self.getReads(), filePath, binaryOutput)
            samConverter.convert()

            # read SAM file
            samfile = pysam.AlignmentFile(filePath, "r")
            reads = list(samfile.fetch())
            self.assertEqual(reads[0].query_name, "SRR622461.77861202")
            # TODO more in-depth testing
            samfile.close()

    def testPlainText(self):
        self._testRoundTrip(False)

    def testBinary(self):
        self._testRoundTrip(True)
