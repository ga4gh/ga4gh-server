"""
Unit tests for variant objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel.variants as variants


class TestGenotypes(unittest.TestCase):
    """
    Unit tests for genotype conversion.
    """
    def verifyGenotypeConversion(
            self, vcfGenotype, vcfPhaseset, callGenotype, callPhaseset):
        """
        Verifies that the convertGenotype function properly converts the vcf
        genotype and phaseset values into the desired call genotype and
        phaseset values.
        """
        self.assertEqual(
            (callGenotype, callPhaseset),
            variants.convertVCFGenotype(vcfGenotype, vcfPhaseset))

    def testGenotypeUnphasedNoCall(self):
        self.verifyGenotypeConversion("./.", "0", [-1], None)

    def testGenotypeUnphasedSecondHalfCall(self):
        self.verifyGenotypeConversion("./0", "25", [-1], None)

    def testGenotypeUnphasedFirstHalfCall(self):
        self.verifyGenotypeConversion("0/.", "", [-1], None)

    def testGenotypeUnphasedRefRef(self):
        self.verifyGenotypeConversion("0/0", "3124234", [0, 0], None)

    def testGenotypeUnphasedAltRef(self):
        self.verifyGenotypeConversion("1/0", "-56809", [1, 0], None)

    def testGenotypeUnphasedRefAlt(self):
        self.verifyGenotypeConversion("0/1", "134965", [0, 1], None)

    def testGenotypePhasedNoCall(self):
        self.verifyGenotypeConversion(".|.", "36", [-1], "36")

    def testGenotypePhasedSecondHalfCall(self):
        self.verifyGenotypeConversion(".|0", "45032", [-1], "45032")

    def testGenotypePhasedFirstHalfCall(self):
        self.verifyGenotypeConversion("0|.", "645", [-1], "645")

    def testGenotypePhasedRefRef(self):
        self.verifyGenotypeConversion("0|0", ".", [0, 0], "*")

    def testGenotypePhasedRefAlt(self):
        self.verifyGenotypeConversion("0|1", "45", [0, 1], "45")

    def testGenotypePhasedAltAlt(self):
        self.verifyGenotypeConversion("1|1", ".", [1, 1], "*")

    def testGenotypePhasedDiffAlt(self):
        self.verifyGenotypeConversion("2|1", "245624", [2, 1], "245624")

    def testPhasesetZero(self):
        self.verifyGenotypeConversion("3|0", "0", [3, 0], "0")

    def testGenotypeHaploid(self):
        self.verifyGenotypeConversion("1", "376", [1], None)
