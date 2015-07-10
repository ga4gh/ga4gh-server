"""
Tests the backend response generators
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.backend as backend
import ga4gh.exceptions as exceptions
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.variants as variants
import ga4gh.protocol as protocol


def generateVariant():
    variant = protocol.Variant()
    return variant


class MockVariantSet(variants.AbstractVariantSet):

    def __init__(self, id_, numVariants):
        super(MockVariantSet, self).__init__(id_)
        self.numVariants = numVariants

    def getVariants(self, reference_name, startPosition, endPosition,
                    variant_name=None, call_set_ids=None):
        for i in range(self.numVariants):
            yield generateVariant()


class TestVariantsGenerator(unittest.TestCase):
    """
    Tests the logic of variantsGenerator
    """
    def setUp(self):
        self.request = protocol.SearchVariantsRequest()
        self.backend = backend.SimulatedBackend()
        self.variant_set_id = "variant_set_id"
        self.dataset_id = self.backend.getDatasetIds()[0]

    def testNoVariantSetsNotSupported(self):
        # a request for no variant sets should throw an exception
        self.request.variant_set_ids = []
        with self.assertRaises(exceptions.NotImplementedException):
            self.backend.variantsGenerator(self.request)

    def testMultipleVariantSetsNotSupported(self):
        # a request for multiple variant sets should throw an exception
        self.request.variant_set_ids = ["1", "2"]
        with self.assertRaises(exceptions.NotImplementedException):
            self.backend.variantsGenerator(self.request)

    def testNonexistantVariantSet(self):
        # a request for a variant set that doesn't exist should throw an error
        self.request.variant_set_ids = ["notFound"]
        with self.assertRaises(exceptions.VariantSetNotFoundException):
            self.backend.variantsGenerator(self.request)

    def testVariantSetEmpty(self):
        # a variant set with no variants should return none
        self._initVariantSet(0)
        iterator = self.backend.variantsGenerator(self.request)
        self.assertIsNone(next(iterator, None))

    def testVariantSetOneVariant(self):
        # a variant set with one variant should return it and a null page_token
        self._initVariantSet(1)
        iterator = self.backend.variantsGenerator(self.request)
        variant, next_page_token = next(iterator)
        self.assertIsNotNone(variant)
        self.assertIsNone(next_page_token)
        self.assertIsNone(next(iterator, None))

    def testVariantSetTwoVariants(self):
        # a variant set with two variants should return the first with
        # a non-null page_token and the second with a null page_token
        self._initVariantSet(2)
        iterator = self.backend.variantsGenerator(self.request)
        variant, next_page_token = next(iterator)
        self.assertIsNotNone(variant)
        self.assertIsNotNone(next_page_token)
        variant, next_page_token = next(iterator)
        self.assertIsNotNone(variant)
        self.assertIsNone(next_page_token)
        self.assertIsNone(next(iterator, None))

    def _initVariantSet(self, numVariants):
        variantSet = MockVariantSet(self.variant_set_id, numVariants)
        self.backend.getDataset(self.dataset_id)._variantSetIdMap = {
            self.variant_set_id: variantSet}
        self.request.variant_set_ids = [self.variant_set_id]


def generateReadAlignment(position=0, sequence='abc'):
    alignment = protocol.ReadAlignment()
    alignment.alignment = protocol.LinearAlignment()
    alignment.alignment.position = protocol.Position()
    alignment.alignment.position.position = position
    alignment.aligned_sequence = sequence
    return alignment


class MockReadGroup(reads.AbstractReadGroup):

    def __init__(self, id_, numAlignments):
        super(MockReadGroup, self).__init__(id_)
        self.numAlignments = numAlignments

    def getReadAlignments(self, reference_name=None, reference_id=None,
                          start=None, end=None):
        for i in range(self.numAlignments):
            yield generateReadAlignment(i)


class TestReadsGenerator(unittest.TestCase):
    """
    Tests the logic of readsGenerator
    """
    def setUp(self):
        self.request = protocol.SearchReadsRequest()
        self.backend = backend.SimulatedBackend()
        self.read_group_id = "read_group_id"
        self.dataset_id = self.backend.getDatasetIds()[0]

    def testNoReadGroupsNotSupported(self):
        # a request for no read groups should throw an exception
        self.request.read_group_ids = []
        with self.assertRaises(exceptions.NotImplementedException):
            self.backend.readsGenerator(self.request)

    def testMultipleReadGroupsNotSupported(self):
        # a request for multiple read groups should throw an exception
        self.request.read_group_ids = ["1", "2"]
        with self.assertRaises(exceptions.NotImplementedException):
            self.backend.readsGenerator(self.request)

    def testNonexistantReadGroup(self):
        # a request for a readGroup that doesn't exist should throw an error
        self.request.read_group_ids = ["notFound"]
        with self.assertRaises(exceptions.ReadGroupNotFoundException):
            self.backend.readsGenerator(self.request)

    def testReadGroupEmpty(self):
        # a readGroup with no reads should return none
        self._initReadGroup(0)
        iterator = self.backend.readsGenerator(self.request)
        self.assertIsNone(next(iterator, None))

    def testReadGroupOneRead(self):
        # a readGroup with one read should return it and a null next_page_token
        self._initReadGroup(1)
        iterator = self.backend.readsGenerator(self.request)
        alignment, next_page_token = next(iterator)
        self.assertIsNotNone(alignment)
        self.assertIsNone(next_page_token)
        self.assertIsNone(next(iterator, None))

    def testReadGroupTwoReads(self):
        # a readGroup with two reads should return the first with
        # a non-null page_token and the second with a null page_token
        self._initReadGroup(2)
        iterator = self.backend.readsGenerator(self.request)
        alignment, next_page_token = next(iterator)
        self.assertIsNotNone(alignment)
        self.assertIsNotNone(next_page_token)
        alignment, next_page_token = next(iterator)
        self.assertIsNotNone(alignment)
        self.assertIsNone(next_page_token)
        self.assertIsNone(next(iterator, None))

    def _initReadGroup(self, numAlignments):
        readGroup = MockReadGroup(self.read_group_id, numAlignments)
        self.backend.getDataset(self.dataset_id)._readGroupIdMap = {
            self.read_group_id: readGroup}
        self.request.read_group_ids = [self.read_group_id]


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
            len(self.read.aligned_sequence), result)
