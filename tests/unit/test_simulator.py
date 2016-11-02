"""
Tests for simulated data
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import datetime

import ga4gh.server.datamodel as datamodel
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.datamodel.reads as reads
import ga4gh.server.datamodel.references as references
import ga4gh.server.datamodel.variants as variants


class TestSimulatedVariantSet(unittest.TestCase):
    """
    Test properties of the SimulatedVariantSet
    """
    def setUp(self):
        self.randomSeed = 0
        self.numCalls = 2
        # ensure variantDensity is >= 1 so we get deterministic behavoir
        self.variantDensity = 1
        self.simulatedVariantSet = self._getSimulatedVariantSet()
        self.referenceName = 'ref'
        self.startPosition = 100
        self.endPosition = 103
        self.callSetIds = ['unused']
        self.bases = ["A", "C", "G", "T"]

    def _getSimulatedVariantSet(self):
        dataset = datasets.Dataset('dataset1')
        referenceSet = references.SimulatedReferenceSet("srs1")
        simulatedVariantSet = variants.SimulatedVariantSet(
            dataset, referenceSet, 'variantSet1', randomSeed=self.randomSeed,
            numCalls=self.numCalls, variantDensity=self.variantDensity)
        return simulatedVariantSet

    def _getSimulatedVariantsList(self, simulatedVariantSet=None):
        if simulatedVariantSet is None:
            simulatedVariantSet = self.simulatedVariantSet
        simulatedVariants = simulatedVariantSet.getVariants(
            self.referenceName, self.startPosition, self.endPosition,
            self.callSetIds)
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
            self.simulatedVariantSet.getCreationTime(),
            self.simulatedVariantSet.getUpdatedTime())

    def testGetVariants(self):
        # calling getVariants should produce the expected results
        variantList = self._getSimulatedVariantsList()
        self.assertEqual(
            len(variantList), self.endPosition - self.startPosition)
        for offset, simulatedVariant in enumerate(variantList):
            start = self.startPosition + offset
            variantSetCompoundId = self.simulatedVariantSet.getCompoundId()
            variantCompoundId = datamodel.VariantCompoundId.parse(
                simulatedVariant.id)
            self.assertEqual(
                variantSetCompoundId.variant_set_id,
                self.simulatedVariantSet.getId())
            self.assertEqual(
                variantSetCompoundId.variant_set_id,
                variantCompoundId.variant_set_id)
            self.assertEqual(
                variantCompoundId.reference_name, self.referenceName)
            self.assertEqual(
                variantCompoundId.start, str(simulatedVariant.start))
            self.assertEqual(
                simulatedVariant.variant_set_id,
                self.simulatedVariantSet.getId())
            self.assertEqual(
                simulatedVariant.reference_name, self.referenceName)
            self.assertEqual(
                simulatedVariant.created, simulatedVariant.updated)
            self.assertEqual(simulatedVariant.start, start)
            self.assertEqual(simulatedVariant.end, start + 1)
            self.assertIn(simulatedVariant.reference_bases, self.bases)
            self.assertIn(
                simulatedVariant.alternate_bases[0], self.bases)
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


class TestSimulatedVariantAnnotationSet(unittest.TestCase):
    def setUp(self):
        self.randomSeed = 1
        self.numCalls = 2
        # ensure variantDensity is >= 1 so we get deterministic behavoir
        self.variantDensity = 1
        self.referenceName = 'ref'
        self.startPosition = 100
        self.endPosition = 120
        self.callSetIds = ['unused']
        self.bases = ["A", "C", "G", "T"]

    def testCreation(self):
        dataset = datasets.Dataset('dataset1')
        referenceSet = references.SimulatedReferenceSet("srs1")
        localId = "variantAnnotationSetId"
        simulatedVariantSet = variants.SimulatedVariantSet(
            dataset, referenceSet, 'variantSet1', randomSeed=self.randomSeed,
            numCalls=self.numCalls, variantDensity=self.variantDensity)
        simulatedVariantAnnotationSet = variants.SimulatedVariantAnnotationSet(
            simulatedVariantSet, localId, self.randomSeed)
        annotations = simulatedVariantAnnotationSet.getVariantAnnotations(
                    self.referenceName, self.startPosition, self.endPosition)
        self.assertEquals(
            simulatedVariantSet.toProtocolElement().id,
            simulatedVariantAnnotationSet.toProtocolElement().variant_set_id,
            "Variant Set ID should match the annotation's variant set ID")
        for variant, ann in annotations:
            self.assertEquals(datetime.datetime.strptime(
                ann.created, "%Y-%m-%dT%H:%M:%S.%fZ").strftime(
                "%Y-%m-%dT%H:%M:%S.%fZ"), ann.created,
                "Expect time format to be in ISO8601")
            self.assertEqual(variant.id, ann.variant_id)


class TestSimulatedReadGroupSet(unittest.TestCase):
    """
    Test properties of the simulated ReadGroupSet
    """
    def testCreation(self):
        dataset = datasets.Dataset('dataset1')
        localId = "readGroupSetId"
        referenceSet = references.SimulatedReferenceSet("srs1")
        simulatedReadGroupSet = reads.SimulatedReadGroupSet(
            dataset, localId, referenceSet)
        for readGroup in simulatedReadGroupSet.getReadGroups():
            alignments = list(readGroup.getReadAlignments())
            self.assertGreater(len(alignments), 0)
