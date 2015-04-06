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


def generateVariant(position=0):
    variant = protocol.GAVariant()
    return variant


class MockVariantSet(variants.AbstractVariantSet):

    def __init__(self, id_, numVariants):
        super(MockVariantSet, self).__init__(id_)
        self.numVariants = numVariants

    def getVariants(self, referenceName, startPosition, endPosition,
                    variantName=None, callSetIds=None):
        for i in range(self.numVariants):
            yield generateVariant(i)


class TestVariantsGenerator(unittest.TestCase):
    """
    Tests the logic of variantsGenerator
    """
    def setUp(self):
        self.request = protocol.GASearchVariantsRequest()
        self.backend = backend.AbstractBackend()
        self.variantSetId = "variantSetId"

    def testNoVariantsNotSupported(self):
        # a request for no variant sets should throw an exception
        self.request.variantSetIds = []
        with self.assertRaises(exceptions.NotImplementedException):
            iterator = self.backend.variantsGenerator(self.request)
            iterator.next()

    def testMultipleVariantsNotSupported(self):
        # a request for multiple variant sets should throw an exception
        self.request.variantSetIds = ["1", "2"]
        with self.assertRaises(exceptions.NotImplementedException):
            iterator = self.backend.variantsGenerator(self.request)
            iterator.next()

    def testNonexistantVariantSet(self):
        # a request for a variant set that doesn't exist should throw an error
        self.request.variantSetIds = ["notFound"]
        with self.assertRaises(exceptions.VariantSetNotFoundException):
            iterator = self.backend.variantsGenerator(self.request)
            iterator.next()

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
        variantSet = MockVariantSet(self.variantSetId, numVariants)
        self.backend._variantSetIdMap = {self.variantSetId: variantSet}
        self.request.variantSetIds = [self.variantSetId]


def generateReadAlignment(position=0, sequence='abc'):
    alignment = protocol.GAReadAlignment()
    alignment.alignment = protocol.GALinearAlignment()
    alignment.alignment.position = protocol.GAPosition()
    alignment.alignment.position.position = position
    alignment.alignedSequence = sequence
    return alignment


class MockReadGroup(reads.AbstractReadGroup):

    def __init__(self, id_, numAlignments):
        super(MockReadGroup, self).__init__(id_)
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
        self.request = protocol.GASearchReadsRequest()
        self.backend = backend.AbstractBackend()
        self.readGroupId = "readGroupId"

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
        self.request.readGroupIds = ["notFound"]
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
        readGroup = MockReadGroup(self.readGroupId, numAlignments)
        self.backend._readGroupIdMap = {self.readGroupId: readGroup}
        self.request.readGroupIds = [self.readGroupId]


class TestReadsGeneratorMethods(unittest.TestCase):
    """
    Test the helper methods for readsGenerator
    """
    def setUp(self):
        self.request = protocol.GASearchReadsRequest()
        self.request.start = 3
        self.request.pageToken = "notNone"
        self.backend = backend.AbstractBackend()
        self.alignmentPosition = 5
        self.startPosition = 1
        self.equalPositionsToSkip = 1

    def testGetPosition(self):
        alignment = generateReadAlignment(self.alignmentPosition)
        position = self.backend._getPosition(alignment)
        self.assertEqual(position, self.alignmentPosition)

    def testGetEndPosition(self):
        alignment = generateReadAlignment(self.alignmentPosition)
        endPosition = self.backend._getEndPosition(alignment)
        self.assertEqual(
            self.alignmentPosition + 3, endPosition)

    def testGetCounters(self):
        # with pageToken == None
        startPos = 9
        self.request.start = startPos
        self.request.pageToken = None
        startPosition, equalPositionsToSkip = self.backend._getCounters(
            self.request)
        self.assertEqual(startPosition, startPos)
        self.assertEqual(equalPositionsToSkip, 0)

        # with pageToken != None
        self.request.pageToken = "1:2"
        startPosition, equalPositionsToSkip = self.backend._getCounters(
            self.request)
        self.assertEqual(startPosition, 1)
        self.assertEqual(equalPositionsToSkip, 2)

    def testPositionIterator(self):
        def getIterator():
            # skipped because position < startPosition
            yield generateReadAlignment(0)
            # skipped because endPosition < request.start
            yield generateReadAlignment(1, 'a')
            # skipped because equalPositionsSkipped
            yield generateReadAlignment(1, 'abcdefg')
            # this one should be returned
            yield generateReadAlignment(1, 'abcdefghijklmnop')
        iterator = getIterator()
        readAlignment = self.backend._positionIterator(
            self.request, iterator, self.startPosition,
            self.equalPositionsToSkip)
        self.assertEqual(readAlignment.alignedSequence, 'abcdefghijklmnop')

    def testPositionIteratorException1(self):
        # test that the first exception in _positionIterator is hit
        def getIterator():
            # skipped because position < startPosition
            yield generateReadAlignment(0)
            # error because no more alignments
        self._assertBadPageTokenExceptionRaised(getIterator())

    def testPositionIteratorException2(self):
        # test that the second exception in _positionIterator is hit
        def getIterator():
            # error because position != startPosition
            yield generateReadAlignment(2)
        self._assertBadPageTokenExceptionRaised(getIterator())

    def testPositionIteratorException3(self):
        # test that the third exception in _positionIterator is hit
        def getIterator():
            # skipped because equalPositionsToSkip == 1
            yield generateReadAlignment(1)
            # error because no more alignments
        self._assertBadPageTokenExceptionRaised(getIterator())

    def _assertBadPageTokenExceptionRaised(self, iterator):
        with self.assertRaises(exceptions.BadPageTokenException):
            self.backend._positionIterator(
                self.request, iterator, self.startPosition,
                self.equalPositionsToSkip)

    def testYieldReads(self):
        firstAlignedSequence = 'abc'
        secondAlignedSequence = 'def'

        def getIterator():
            yield generateReadAlignment(2, firstAlignedSequence)

        iterator = getIterator()
        readAlignment = generateReadAlignment(2, secondAlignedSequence)
        equalPositionsToSkip = 0
        readIterator = self.backend._yieldReads(
            iterator, readAlignment, equalPositionsToSkip)
        readAlignment, nextPageToken = next(readIterator)
        self.assertEqual(secondAlignedSequence, readAlignment.alignedSequence)
        self.assertEqual(nextPageToken, "2:1")
        readAlignment, nextPageToken = next(readIterator)
        self.assertIsNone(next(readIterator, None))
