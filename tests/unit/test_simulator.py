"""
Tests for simulated data
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.variants as variants


class TestSimulatedVariantSet(unittest.TestCase):
    """
    Test properties of the SimulatedVariantSet
    """
    def setUp(self):
        self.randomSeed = 0
        self.numCalls = 2
        # ensure variantDensity is >= 1 so we get deterministic behavoir
        self.variantDensity = 1
        self.variantSetId = 3
        self.simulatedVariantSet = self._getSimulatedVariantSet()
        self.referenceName = 'ref'
        self.startPosition = 100
        self.endPosition = 103
        self.variantName = 'unused'
        self.callSetIds = ['unused']
        self.bases = ["A", "C", "G", "T"]

    def _getSimulatedVariantSet(self):
        simulatedVariantSet = variants.SimulatedVariantSet(
            self.randomSeed, self.numCalls,
            self.variantDensity, self.variantSetId)
        return simulatedVariantSet

    def _getSimulatedVariantsList(self, simulatedVariantSet=None):
        if simulatedVariantSet is None:
            simulatedVariantSet = self.simulatedVariantSet
        simulatedVariants = simulatedVariantSet.getVariants(
            self.referenceName, self.startPosition, self.endPosition,
            self.variantName, self.callSetIds)
        variantList = list(simulatedVariants)
        return variantList

    def testConstruction(self):
        # initializing SimulatedVariantSet should store values correctly
        self.assertEqual(
            self.randomSeed, self.simulatedVariantSet._randomSeed)
        self.assertEqual(
            self.numCalls, self.simulatedVariantSet._numCalls)
        self.assertEqual(
            self.variantDensity, self.simulatedVariantSet._variantDensity)
        self.assertEqual(
            self.variantSetId, self.simulatedVariantSet._id)
        self.assertEqual(
            self.simulatedVariantSet.getCreationTime(),
            self.simulatedVariantSet.getUpdatedTime())

    def testGetVariants(self):
        # calling getVariants should produce the expected results
        variantList = self._getSimulatedVariantsList()
        self.assertEqual(
            len(variantList), self.endPosition - self.startPosition)
        simulatedVariant = variantList[0]
        self.assertIn(str(self.variantSetId), simulatedVariant.id)
        self.assertIn(self.referenceName, simulatedVariant.id)
        self.assertIn(str(self.startPosition), simulatedVariant.id)
        self.assertEqual(simulatedVariant.variantSetId, self.variantSetId)
        self.assertEqual(simulatedVariant.referenceName, self.referenceName)
        self.assertEqual(simulatedVariant.created, simulatedVariant.updated)
        self.assertEqual(simulatedVariant.start, self.startPosition)
        self.assertEqual(simulatedVariant.end, self.startPosition + 1)
        self.assertIn(simulatedVariant.referenceBases, self.bases)
        self.assertIn(
            simulatedVariant.alternateBases[0], self.bases)
        self.assertEqual(len(simulatedVariant.calls), self.numCalls)

    def testConsistency(self):
        # two SimulatedBackend objects given the same parameters
        # should produce identical variant lists
        variantListOne = self._getSimulatedVariantsList()
        simulatedVariantSetTwo = self._getSimulatedVariantSet()
        variantListTwo = self._getSimulatedVariantsList(
            simulatedVariantSetTwo)
        self._assertEqualVariantLists(variantListOne, variantListTwo)

    def testSelfConsistent(self):
        # the same SimulatedBackend should produce identical
        # variant lists given the same parameters
        variantListOne = self._getSimulatedVariantsList()
        variantListTwo = self._getSimulatedVariantsList()
        self.assertEqual(variantListOne, variantListTwo)

    def _assertEqualVariantLists(self, variantListOne, variantListTwo):
        # need to make time-dependent fields equal before the comparison,
        # otherwise we're introducing a race condition
        timeDependentFields = ['created', 'updated']
        for variantList in (variantListOne, variantListTwo):
            for variant in variantList:
                for field in timeDependentFields:
                    setattr(variant, field, 0)
        self.assertEqual(variantListOne, variantListTwo)


class TestSimulatedReadGroupSet(unittest.TestCase):
    """
    Test properties of the simulated ReadGroupSet
    """
    def testCreation(self):
        readGroupSetId = "readGroupSetId"
        simulatedReadGroupSet = reads.SimulatedReadGroupSet(readGroupSetId)
        for readGroup in simulatedReadGroupSet.getReadGroups():
            alignments = list(readGroup.getReadAlignments())
            self.assertGreater(len(alignments), 0)
