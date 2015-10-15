"""
Unit tests for variant objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.exceptions as exceptions
import ga4gh.backend as backend
import ga4gh.datamodel.variants as variants
import ga4gh.datamodel.datasets as datasets


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


class TestAbstractVariantSet(unittest.TestCase):
    """
    Unit tests for the abstract variant set.
    """
    def setUp(self):
        self._variantSetName = "testVariantSet"
        self._backend = backend.AbstractBackend()
        self._dataset = datasets.AbstractDataset(self._backend)
        self._variantSet = variants.AbstractVariantSet(
            self._dataset, self._variantSetName)

    def testAddOneCallSet(self):
        self.assertEqual(self._variantSet.getNumCallSets(), 0)
        callSetName = "callSetName"
        self._variantSet.addCallSet(callSetName)
        self.assertEqual(self._variantSet.getNumCallSets(), 1)
        callSet = self._variantSet.getCallSetByIndex(0)
        self.assertEqual(
            self._variantSet.getCallSetByName(callSetName), callSet)
        self.assertEqual(
            self._variantSet.getCallSet(callSet.getId()), callSet)
        self.assertEqual(self._variantSet.getCallSets(), [callSet])

    def testAddMultipleCallSets(self):
        self.assertEqual(self._variantSet.getNumCallSets(), 0)
        callSetName = "callSetName{}"
        callSetCount = 10
        allCallSets = []
        for i in range(callSetCount):
            self._variantSet.addCallSet(callSetName.format(i))
            callSet = self._variantSet.getCallSetByIndex(i)
            self.assertEqual(
                self._variantSet.getCallSetByName(
                    callSetName.format(i)), callSet)
            self.assertEqual(
                self._variantSet.getCallSet(callSet.getId()), callSet)
            allCallSets.append(callSet)
        self.assertEqual(self._variantSet.getNumCallSets(), callSetCount)
        self.assertEqual(self._variantSet.getCallSets(), allCallSets)

    def testCallSetMethods(self):
        # test invalid lookups
        self.assertRaises(IndexError, self._variantSet.getCallSetByIndex, 10)
        self.assertRaises(exceptions.CallSetNameNotFoundException,
                          self._variantSet.getCallSetByName, 'noname')
        self.assertRaises(exceptions.CallSetNameNotFoundException,
                          self._variantSet.getCallSetByName, None)
        self.assertRaises(exceptions.CallSetNotFoundException,
                          self._variantSet.getCallSet, 617)
        self.assertRaises(exceptions.CallSetNotFoundException,
                          self._variantSet.getCallSet, None)
        self.assertRaises(TypeError, self._variantSet.addCallSet,
                          ['a list of', 2])
        self.assertRaises(NotImplementedError,
                          self._variantSet.getNumVariants)

    def testGetVariantId(self):
        self.assertRaises(AttributeError,
                          self._variantSet.getVariantId, None)
        self.assertRaises(AttributeError,
                          self._variantSet.getVariantId, "hola")

    def testHashVariant(self):
        self.assertRaises(AttributeError,
                          self._variantSet.hashVariant, None)
        self.assertRaises(AttributeError,
                          self._variantSet.hashVariant, "hi")

    def testVariantSetProtocolElement(self):
        # 'AbstractVariantSet' object has no attribute 'getMetadata'
        self.assertRaises(AttributeError,
                          self._variantSet.toProtocolElement)
