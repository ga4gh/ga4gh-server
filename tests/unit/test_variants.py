"""
Unit tests for variant objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.server.exceptions as exceptions
import ga4gh.server.datamodel.variants as variants
import ga4gh.server.datamodel.datasets as datasets


class TestAbstractVariantSet(unittest.TestCase):
    """
    Unit tests for the abstract variant set.
    """
    def setUp(self):
        self._variantSetName = "testVariantSet"
        self._dataset = datasets.Dataset("datasetId")
        self._variantSet = variants.AbstractVariantSet(
            self._dataset, self._variantSetName)

    def testAddOneCallSet(self):
        self.assertEqual(self._variantSet.getNumCallSets(), 0)
        callSetName = "callSetName"
        self._variantSet.addCallSetFromName(callSetName)
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
            self._variantSet.addCallSetFromName(callSetName.format(i))
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
        self.assertRaises(AttributeError,
                          self._variantSet.toProtocolElement)
