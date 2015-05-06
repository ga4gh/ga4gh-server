"""
Tests the converters
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import mock
import tempfile
import unittest

import pysam

import ga4gh.protocol as protocol
import ga4gh.converters as converters
import tests.utils as utils


class TestSamConverter(unittest.TestCase):
    """
    Base class for SAM converter tests.
    Provides common  methods.
    """
    readsFilePath = 'tests/unit/reads.dat'

    def _getReads(self):
        readsFile = file(self.readsFilePath)
        lines = readsFile.readlines()
        reads = []
        for line in lines:
            read = protocol.ReadAlignment.fromJsonString(line)
            reads.append(read)
        return reads

    def getHttpClient(self):
        httpClient = utils.makeHttpClient()
        httpClient.searchReads = self.getSearchReadsResponse
        return httpClient

    def getSearchReadsRequest(self):
        request = protocol.SearchReadsRequest()
        return request

    def getSearchReadsResponse(self, request):
        return self._getReads()


class TestSamConverterLogic(TestSamConverter):
    """
    Test the SamConverter logic
    """
    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testSamConverter(self):
        mockPysam = mock.Mock()
        with mock.patch('pysam.AlignmentFile', mockPysam):
            httpClient = self.getHttpClient()
            searchReadsRequest = self.getSearchReadsRequest()
            outputFile = None
            binaryOutput = False
            samConverter = converters.SamConverter(
                httpClient, searchReadsRequest, outputFile, binaryOutput)
            samConverter.convert()


class TestSamConverterRoundTrip(TestSamConverter):
    """
    Write a sam file and see if pysam can read it
    """
    def _testRoundTrip(self, binaryOutput):
        with tempfile.NamedTemporaryFile() as fileHandle:
            # write SAM file
            httpClient = self.getHttpClient()
            searchReadsRequest = self.getSearchReadsRequest()
            filePath = fileHandle.name
            samConverter = converters.SamConverter(
                httpClient, searchReadsRequest, filePath, binaryOutput)
            samConverter.convert()

            # read SAM file
            samfile = pysam.AlignmentFile(filePath, "r")
            reads = list(samfile.fetch())
            self.assertEqual(reads[0].query_name, "SRR622461.77861202")
            # TODO more in-depth testing
            samfile.close()

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testPlainText(self):
        self._testRoundTrip(False)

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testBinary(self):
        self._testRoundTrip(True)
