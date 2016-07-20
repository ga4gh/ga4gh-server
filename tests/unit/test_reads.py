"""
Unit tests for reads objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datasource.htslib as htslib
import ga4gh.protocol as protocol


class TestParseMalformedBamHeader(unittest.TestCase):
    """
    Tests for parsing of malformed bam headers.

    htslib.parseMalformedBamHeader should not modify correct parsed
    headers and should parse out additional fields separated by spaces
    (instead of tabs as defined in the SAM spec).
    """

    def testGoodHeaderUnmodified(self):
        header = {'SO': 'coordinate', 'VN': '1.0'}
        self.assertEqual(header, htslib.parseMalformedBamHeader(header))

    def testOriginalTypesUnmodified(self):
        # note real field tags, just checking that types are preserved
        header = {'int': 2845856850,
                  'float': 206.6,
                  'bool': True,
                  'string': '123'}
        self.assertEqual(header, htslib.parseMalformedBamHeader(header))

    def testCommandsWithSpacesNotParsed(self):
        header = {'CL': 'bwa aln -q 15 -f $sai_file ' +
                        '$reference_fasta $fastq_file\tPP:bwa_index',
                  'ID': 'bwa_aln_fastq',
                  'PN': 'bwa',
                  'VN': '0.5.9-r16'}
        self.assertEqual(header, htslib.parseMalformedBamHeader(header))

    def testSpaceSeparatedUnparsedFieldsParsed(self):
        header = {'LN': 249250621,
                  'M5': '1b22b98cdeb4a9304cb5d48026a85128',
                  'SN': '1',
                  'UR': 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/' +
                        'technical/reference/phase2_reference_assembly' +
                        '_sequence/hs37d5.fa.gz        AS:NCBI37' +
                        '       SP:Human'}
        expected = {'LN': 249250621,
                    'M5': '1b22b98cdeb4a9304cb5d48026a85128',
                    'SN': '1',
                    'UR': 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/' +
                          'technical/reference/phase2_reference_assembly' +
                          '_sequence/hs37d5.fa.gz',
                    'AS': 'NCBI37',
                    'SP': 'Human'}
        self.assertEqual(expected, htslib.parseMalformedBamHeader(header))


class TestSamCigar(unittest.TestCase):
    """
    Test Sam Cigar class handles Cigar mappings correctly

    The integer codes are defined in the SAM spec. Thus, the ordering of
    SamCigar.cigarStrings implicitly implements this spec.
    """

    def testAlignmentMatch(self):
        self.assertEqual(0, htslib.SamCigar.ga2int(
            protocol.CigarUnit.ALIGNMENT_MATCH))

        self.assertEqual(protocol.CigarUnit.ALIGNMENT_MATCH,
                         htslib.SamCigar.int2ga(0))

    def testInsertion(self):
        self.assertEqual(1, htslib.SamCigar.ga2int(
            protocol.CigarUnit.INSERT))

        self.assertEqual(protocol.CigarUnit.INSERT,
                         htslib.SamCigar.int2ga(1))

    def testDeletion(self):
        self.assertEqual(2, htslib.SamCigar.ga2int(
            protocol.CigarUnit.DELETE))

        self.assertEqual(protocol.CigarUnit.DELETE,
                         htslib.SamCigar.int2ga(2))

    def testSkipped(self):
        self.assertEqual(3, htslib.SamCigar.ga2int(
            protocol.CigarUnit.SKIP))

        self.assertEqual(protocol.CigarUnit.SKIP,
                         htslib.SamCigar.int2ga(3))

    def testSoftClipping(self):
        self.assertEqual(4, htslib.SamCigar.ga2int(
            protocol.CigarUnit.CLIP_SOFT))

        self.assertEqual(protocol.CigarUnit.CLIP_SOFT,
                         htslib.SamCigar.int2ga(4))

    def testHardClipping(self):
        self.assertEqual(5, htslib.SamCigar.ga2int(
            protocol.CigarUnit.CLIP_HARD))

        self.assertEqual(protocol.CigarUnit.CLIP_HARD,
                         htslib.SamCigar.int2ga(5))

    def testPadding(self):
        self.assertEqual(6, htslib.SamCigar.ga2int(
            protocol.CigarUnit.PAD))

        self.assertEqual(protocol.CigarUnit.PAD,
                         htslib.SamCigar.int2ga(6))

    def testSequenceMatch(self):
        self.assertEqual(7, htslib.SamCigar.ga2int(
            protocol.CigarUnit.SEQUENCE_MATCH))

        self.assertEqual(protocol.CigarUnit.SEQUENCE_MATCH,
                         htslib.SamCigar.int2ga(7))

    def testSequenceMismatch(self):
        self.assertEqual(8, htslib.SamCigar.ga2int(
            protocol.CigarUnit.SEQUENCE_MISMATCH))

        self.assertEqual(protocol.CigarUnit.SEQUENCE_MISMATCH,
                         htslib.SamCigar.int2ga(8))


class TestSamFlags(unittest.TestCase):
    """
    Tests SamFlags utilities for checking the status of and
    setting flags.

    Flags are defined by the SAM spec.
    """

    def setUp(self):
        self.flag = 0x0

    def testPairedReadFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.READ_PAIRED)
        self.assertEqual(0x1, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.READ_PAIRED))

    def testProperPairReadFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.READ_PROPER_PAIR)
        self.assertEqual(0x2, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.READ_PROPER_PAIR))

    def testUnmappedReadFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.READ_UNMAPPED)
        self.assertEqual(0x4, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.READ_UNMAPPED))

    def testUnmappedMateFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.MATE_UNMAPPED)
        self.assertEqual(0x8, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.MATE_UNMAPPED))

    def testReverseStrandReadFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.READ_REVERSE_STRAND)
        self.assertEqual(0x10, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.READ_REVERSE_STRAND))

    def testReverseStrandMateFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.MATE_REVERSE_STRAND)
        self.assertEqual(0x20, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.MATE_REVERSE_STRAND))

    def testFirstPairFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.FIRST_IN_PAIR)
        self.assertEqual(0x40, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.FIRST_IN_PAIR))

    def testSecondPairFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.SECOND_IN_PAIR)
        self.assertEqual(0x80, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.SECOND_IN_PAIR))

    def testSecondaryAlignmentFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.SECONDARY_ALIGNMENT)
        self.assertEqual(0x100, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.SECONDARY_ALIGNMENT))

    def testFailedQualityCheckFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.FAILED_QUALITY_CHECK)
        self.assertEqual(0x200, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.FAILED_QUALITY_CHECK))

    def testDuplicateReadFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.DUPLICATE_READ)
        self.assertEqual(0x400, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.DUPLICATE_READ))

    def testSupplementaryAlignmentFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.SUPPLEMENTARY_ALIGNMENT)
        self.assertEqual(0x800, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.SUPPLEMENTARY_ALIGNMENT))

    def testFlagNotSet(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.READ_PAIRED)
        self.assertFalse(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.READ_REVERSE_STRAND))

    def testComboFlag(self):
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.READ_PAIRED)
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.FIRST_IN_PAIR)
        self.flag = htslib.SamFlags.setFlag(
            self.flag, htslib.SamFlags.FAILED_QUALITY_CHECK)
        self.assertEqual(0x241, self.flag)
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.READ_PAIRED))
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.FIRST_IN_PAIR))
        self.assertTrue(htslib.SamFlags.isFlagSet(
            self.flag, htslib.SamFlags.FAILED_QUALITY_CHECK))
