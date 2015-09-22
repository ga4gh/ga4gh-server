"""
Tests the backend response generators
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.backend as backend
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.variants as variants
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol


def generateVariant():
    variant = protocol.Variant()
    return variant


class MockVariantSet(variants.AbstractVariantSet):

    def __init__(self, parentContainer, localId, numVariants):
        super(MockVariantSet, self).__init__(parentContainer, localId)
        self.numVariants = numVariants

    def getVariants(self, referenceName, startPosition, endPosition,
                    variantName=None, callSetIds=None):
        for i in range(self.numVariants):
            yield generateVariant()


class TestVariantsGenerator(unittest.TestCase):
    """
    Tests the logic of variantsGenerator
    """
    def setUp(self):
        self.request = protocol.SearchVariantsRequest()
        self.backend = backend.SimulatedBackend()
        self.dataset = self.backend.getDatasets()[0]

    def testNonexistantVariantSet(self):
        # a request for a variant set that doesn't exist should throw an error
        variantSet = variants.AbstractVariantSet(
            self.dataset, 'notFound')
        self.request.variantSetId = variantSet.getId()
        with self.assertRaises(exceptions.VariantSetNotFoundException):
            self.backend.variantsGenerator(self.request)

    def testVariantSetEmpty(self):
        # a variant set with no variants should return none
        self._initVariantSet(0)
        iterator = self.backend.variantsGenerator(self.request)
        self.assertIsNone(next(iterator, None))

    def testVariantSetOneVariant(self):
        # a variant set with one variant should return it and a null pageToken
        self._initVariantSet(1)
        iterator = self.backend.variantsGenerator(self.request)
        variant, nextPageToken = next(iterator)
        self.assertIsNotNone(variant)
        self.assertIsNone(nextPageToken)
        self.assertIsNone(next(iterator, None))

    def testVariantSetTwoVariants(self):
        # a variant set with two variants should return the first with
        # a non-null pageToken and the second with a null pageToken
        self._initVariantSet(2)
        iterator = self.backend.variantsGenerator(self.request)
        variant, nextPageToken = next(iterator)
        self.assertIsNotNone(variant)
        self.assertIsNotNone(nextPageToken)
        variant, nextPageToken = next(iterator)
        self.assertIsNotNone(variant)
        self.assertIsNone(nextPageToken)
        self.assertIsNone(next(iterator, None))

    def _initVariantSet(self, numVariants):
        variantSet = MockVariantSet(
            self.dataset, "mockvs", numVariants)
        self.dataset.addVariantSet(variantSet)
        self.request.variantSetId = variantSet.getId()


def generateReadAlignment(position=0, sequence='abc'):
    alignment = protocol.ReadAlignment()
    alignment.alignment = protocol.LinearAlignment()
    alignment.alignment.position = protocol.Position()
    alignment.alignment.position.position = position
    alignment.alignedSequence = sequence
    return alignment


class MockReadGroup(reads.AbstractReadGroup):

    def __init__(self, parentContainer, localId, numAlignments):
        super(MockReadGroup, self).__init__(parentContainer, localId)
        self.numAlignments = numAlignments

    def getReadAlignments(self, referenceName=None, referenceId=None,
                          start=None, end=None):
        for i in range(self.numAlignments):
            yield generateReadAlignment(i)


class TestReadsGenerator(unittest.TestCase):
    """
    Tests the logic of readsGenerator
    """
    def setUp(self):
        self.request = protocol.SearchReadsRequest()
        self.backend = backend.SimulatedBackend(numAlignments=0)
        referenceSet = self.backend.getReferenceSetByIndex(0)
        reference = referenceSet.getReferenceByIndex(0)
        self.request.referenceId = reference.getId()
        self.dataset = self.backend.getDatasets()[0]
        self.readGroupSet = self.dataset.getReadGroupSets()[0]

    def testNoReadGroupsNotSupported(self):
        # a request for no read groups should throw an exception
        self.request.readGroupIds = []
        with self.assertRaises(exceptions.NotImplementedException):
            self.backend.readsGenerator(self.request)

    def testMultipleReadGroupsNotSupported(self):
        # a request for multiple read groups should throw an exception
        self.request.readGroupIds = ["1", "2"]
        with self.assertRaises(exceptions.NotImplementedException):
            self.backend.readsGenerator(self.request)

    def testNonexistantReadGroup(self):
        # a request for a readGroup that doesn't exist should throw an error
        readGroup = reads.AbstractReadGroup(self.readGroupSet, 'notFound')
        self.request.readGroupIds = [readGroup.getId()]
        with self.assertRaises(exceptions.ReadGroupNotFoundException):
            self.backend.readsGenerator(self.request)

    def testReadGroupEmpty(self):
        # a readGroup with no reads should return none
        self._initReadGroup(0)
        iterator = self.backend.readsGenerator(self.request)
        self.assertIsNone(next(iterator, None))

    def testReadGroupOneRead(self):
        # a readGroup with one read should return it and a null nextPageToken
        self._initReadGroup(1)
        iterator = self.backend.readsGenerator(self.request)
        alignment, nextPageToken = next(iterator)
        self.assertIsNotNone(alignment)
        self.assertIsNone(nextPageToken)
        self.assertIsNone(next(iterator, None))

    def testReadGroupTwoReads(self):
        # a readGroup with two reads should return the first with
        # a non-null pageToken and the second with a null pageToken
        self._initReadGroup(2)
        iterator = self.backend.readsGenerator(self.request)
        alignment, nextPageToken = next(iterator)
        self.assertIsNotNone(alignment)
        self.assertIsNotNone(nextPageToken)
        alignment, nextPageToken = next(iterator)
        self.assertIsNotNone(alignment)
        self.assertIsNone(nextPageToken)
        self.assertIsNone(next(iterator, None))

    def _initReadGroup(self, numAlignments):
        readGroup = MockReadGroup(
            self.readGroupSet, "mockrg", numAlignments)
        self.readGroupSet.addReadGroup(readGroup)
        self.request.readGroupIds = [readGroup.getId()]


class TestVariantsIntervalIteratorClassMethods(unittest.TestCase):
    """
    Test the variants interval iterator class methods
    """
    def setUp(self):
        self.variant = protocol.Variant()
        self.variant.start = 4
        self.variant.end = 6
        self.intervalIterator = backend.VariantsIntervalIterator

    def testGetVariantStart(self):
        result = self.intervalIterator._getStart(self.variant)
        self.assertEqual(self.variant.start, result)

    def testGetVariantEnd(self):
        result = self.intervalIterator._getEnd(self.variant)
        self.assertEqual(self.variant.end, result)


class TestReadsIntervalIteratorClassMethods(unittest.TestCase):
    """
    Test the variants interval iterator class methods
    """
    def setUp(self):
        self.read = generateReadAlignment(5)
        self.intervalIterator = backend.ReadsIntervalIterator

    def testGetReadStart(self):
        result = self.intervalIterator._getStart(self.read)
        self.assertEqual(self.read.alignment.position.position, result)

    def testGetReadEnd(self):
        result = self.intervalIterator._getEnd(self.read)
        self.assertEqual(
            self.intervalIterator._getStart(self.read) +
            len(self.read.alignedSequence), result)
