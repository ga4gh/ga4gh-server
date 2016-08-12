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
import ga4gh.datarepo as datarepo
import tests.paths as paths
import tests.utils as utils


class TestSamConverter(unittest.TestCase):
    """
    Tests for the GA4GH reads API -> SAM conversion.
    """
    def setUp(self):
        dataRepository = datarepo.SqlDataRepository(paths.testDataRepo)
        dataRepository.open(datarepo.MODE_READ)
        self._backend = backend.Backend(dataRepository)
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
                self._client, readGroup.getId(), reference.getId(),
                outputFileName=fileHandle.name)
            converter.convert()
            samFile = pysam.AlignmentFile(fileHandle.name, "r")
            try:
                # TODO suppressed because of pysam output:
                # [W::sam_parse1] mapped mate cannot have zero coordinate;
                # treated as unmapped
                # and
                # [W::sam_parse1] mapped mate cannot have zero coordinate;
                # treated as unmapped
                # see discussion in https://github.com/ga4gh/server/pull/789
                with utils.suppressOutput():
                    convertedReads = list(samFile.fetch())
            finally:
                samFile.close()
            samFile = pysam.AlignmentFile(readGroupSet.getDataUrl(), "rb")
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
        datasets = self._backend.getDataRepository().getDatasets()
        for dataset in datasets:
            readGroupSets = dataset.getReadGroupSets()
            for readGroupSet in readGroupSets:
                referenceSet = readGroupSet.getReferenceSet()
                for reference in referenceSet.getReferences():
                    for readGroup in readGroupSet.getReadGroups():
                        self.verifyFullConversion(
                            readGroupSet, readGroup, reference)
