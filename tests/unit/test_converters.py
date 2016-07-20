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
import ga4gh.registry as registry
import tests.paths as paths


class TestSamConverter(unittest.TestCase):
    """
    Tests for the GA4GH reads API -> SAM conversion.
    """
    def setUp(self):
        self._registry_db = registry.RegistryDb(paths.testDataRepoUrl)
        self._registry_db.open()
        self._backend = backend.Backend(self._registry_db)
        self._client = client.LocalClient(self._backend)

    def verifySamRecordsEqual(self, sourceReads, convertedReads):
        """
        Verify that a read from pysam matches a read from the reference server
        """
        self.assertEqual(len(sourceReads), len(convertedReads))
        for source, converted in zip(sourceReads, convertedReads):
            self.assertEqual(source.query_name, converted.query_name)
            self.assertEqual(source.query_sequence, converted.query_sequence)
            self.assertEqual(source.flag, converted.flag)
            self.assertEqual(source.reference_id, converted.reference_id)
            self.assertEqual(
                source.mapping_quality,
                converted.mapping_quality)
            self.assertEqual(
                source.template_length,
                converted.template_length)
            self.assertEqual(
                source.query_qualities, converted.query_qualities)
            # TODO the below fields can not be tested since we don't
            # encode them in the case that either the read is not mapped
            # or the read pair is not mapped
            # self.assertEqual(
            #     source.reference_start,
            #     converted.reference_start)
            # self.assertEqual(source.cigar, converted.cigar)
            # self.assertEqual(
            #     source.next_reference_id,
            #     converted.next_reference_id)
            # self.assertEqual(
            #     source.next_reference_start,
            #     converted.next_reference_start)
            # TODO can't uncomment until round trip tags are fixed;
            # see schemas issue 758
            # self.assertEqual(
            #     source.tags,
            #     converted.tags)

    def verifyFullConversion(self, readGroupSet, readGroup, reference):
        """
        Verify that the conversion of the specified readGroup in the
        specified readGroupSet for the specified reference is correct.
        This involves pulling out the reads from the original BAM file
        and comparing these with the converted SAM records.
        """
        with tempfile.NamedTemporaryFile() as fileHandle:
            converter = converters.SamConverter(
                self._client, str(readGroup.id), str(reference.id),
                outputFileName=fileHandle.name)
            converter.convert()
            samFile = pysam.AlignmentFile(fileHandle.name, "r")
            try:
                convertedReads = list(samFile.fetch())
            finally:
                samFile.close()
            samFile = pysam.AlignmentFile(readGroupSet.data_url, "rb")
            try:
                sourceReads = []
                referenceName = reference.name.encode()
                readGroupName = readGroup.name.encode()
                for readAlignment in samFile.fetch(referenceName):
                    tags = dict(readAlignment.tags)
                    if 'RG' in tags and tags['RG'] == readGroupName:
                        sourceReads.append(readAlignment)
            finally:
                samFile.close()
            self.verifySamRecordsEqual(sourceReads, convertedReads)

    def testSamConversion(self):
        for dataset in self._registry_db.get_datasets():
            for read_group_set in dataset.read_group_sets:
                reference_set = read_group_set.reference_set
                for reference in reference_set.references:
                    for read_group in read_group_set.read_groups:
                        self.verifyFullConversion(
                            read_group_set, read_group, reference)
