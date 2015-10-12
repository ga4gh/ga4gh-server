"""
Tests the converters
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import tempfile
import unittest

import pysam

import ga4gh.backend as backend
import ga4gh.client as client
import ga4gh.converters as converters


class TestSamConverter(unittest.TestCase):
    """
    Tests for the GA4GH reads API -> SAM conversion.
    """
    def setUp(self):
        self._backend = backend.FileSystemBackend("tests/data")
        self._client = client.LocalClient(self._backend)

    def verifySamRecordsEqual(self, sourceReads, convertedReads):
        """
        Verify that a read from pysam matches a read from the reference server
        """
        self.assertEqual(len(sourceReads), len(convertedReads))
        for source, converted in zip(sourceReads, convertedReads):
            self.assertEqual(source.query_name, converted.query_name)
            self.assertEqual(source.query_sequence, converted.query_sequence)
            self.assertEqual(source.qual, converted.qual)
            self.assertEqual(source.flag, converted.flag)
            # TODO add more comparisons.

    def verifyFullConversion(self, readGroupSet, readGroup, reference):
        """
        Verify that the conversion of the specified readGroup in the
        specified readGroupSet for the specified reference is correct.
        This involves pulling out the reads from the original BAM file
        and comparing these with the converted SAM records.
        """
        with tempfile.NamedTemporaryFile() as fileHandle:
            converter = converters.SamConverter(
                self._client, readGroup.getId(), reference.getId(),
                outputFileName=fileHandle.name)
            converter.convert()
            samFile = pysam.AlignmentFile(fileHandle.name, "r")
            try:
                convertedReads = list(samFile.fetch())
            finally:
                samFile.close()
            samFile = pysam.AlignmentFile(readGroupSet.getSamFilePath(), "rb")
            try:
                sourceReads = []
                referenceName = reference.getName().encode()
                readGroupName = readGroup.getLocalId().encode()
                for readAlignment in samFile.fetch(referenceName):
                    tags = dict(readAlignment.tags)
                    if 'RG' in tags and tags['RG'] == readGroupName:
                        sourceReads.append(readAlignment)
            finally:
                samFile.close()
            self.verifySamRecordsEqual(sourceReads, convertedReads)

    def testSamConversion(self):
        # TODO generalise to all readGroupSets in the test data backend.
        dataset = self._backend.getDatasetByIndex(0)
        readGroupSet = dataset.getReadGroupSetByIndex(0)
        referenceSet = readGroupSet.getReferenceSet()
        for reference in referenceSet.getReferences():
            for readGroup in readGroupSet.getReadGroups():
                self.verifyFullConversion(readGroupSet, readGroup, reference)
