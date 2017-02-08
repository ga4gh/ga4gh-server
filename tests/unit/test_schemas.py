from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.server.response_builder as response_builder
import ga4gh.schemas.protocol as protocol


def getValueListName(cls):
    return [field for field in cls.DESCRIPTOR.fields_by_name
            if field != 'next_page_token'][0]


class SearchResponseBuilderTest(unittest.TestCase):
    """
    Tests the SearchResponseBuilder class to ensure that it behaves
    correctly.
    """
    def testIntegrity(self):
        # Verifies that the values we put in are exactly what we get
        # back across all subclasses of SearchResponse
        for class_ in [responseClass for _, _, responseClass in
                       protocol.postMethods]:
            instance = class_()
            valueList = getattr(instance, getValueListName(class_))
            valueList.add()
            builder = response_builder.SearchResponseBuilder(
                class_, len(valueList), 2 ** 32)
            for value in valueList:
                builder.addValue(value)
            builder.setNextPageToken(instance.next_page_token)
            otherInstance = protocol.fromJson(
                builder.getSerializedResponse(), class_)
            self.assertEqual(instance, otherInstance)

    def testPageSizeOverflow(self):
        # Verifies that the page size behaviour is correct when we keep
        # filling after full is True.
        responseClass = protocol.SearchVariantsResponse
        valueClass = protocol.Variant
        for pageSize in range(1, 10):
            builder = response_builder.SearchResponseBuilder(
                responseClass, pageSize, 2 ** 32)
            self.assertEqual(builder.getPageSize(), pageSize)
            self.assertFalse(builder.isFull())
            for listLength in range(1, 2 * pageSize):
                builder.addValue(valueClass())
                instance = protocol.fromJson(
                    builder.getSerializedResponse(), responseClass)

                valueList = getattr(instance, getValueListName(responseClass))
                self.assertEqual(len(valueList), listLength)
                if listLength < pageSize:
                    self.assertFalse(builder.isFull())
                else:
                    self.assertTrue(builder.isFull())

    def testPageSizeExactFill(self):
        responseClass = protocol.SearchVariantsResponse
        valueClass = protocol.Variant
        for pageSize in range(1, 10):
            builder = response_builder.SearchResponseBuilder(
                responseClass, pageSize, 2 ** 32)
            self.assertEqual(builder.getPageSize(), pageSize)
            while not builder.isFull():
                builder.addValue(valueClass())
            instance = protocol.fromJson(builder.getSerializedResponse(),
                                         responseClass)
            valueList = getattr(instance, getValueListName(responseClass))
            self.assertEqual(len(valueList), pageSize)

    def testMaxBufferSizeOverridesPageSize(self):
        responseClass = protocol.SearchVariantsResponse
        typicalValue = protocol.Variant()
        # We have to put some values in here or it will have zero length.
        typicalValue.start = 1
        typicalValue.end = 2
        typicalValue.reference_bases = "AAAAAAAA"
        typicalValueLength = typicalValue.ByteSize()
        for numValues in range(1, 10):
            maxBufferSize = numValues * typicalValueLength
            builder = response_builder.SearchResponseBuilder(
                responseClass, 1000, maxBufferSize)
            self.assertEqual(
                maxBufferSize, builder.getMaxBufferSize())
            while not builder.isFull():
                builder.addValue(typicalValue)
            instance = protocol.fromJson(builder.getSerializedResponse(),
                                         responseClass)
            valueList = getattr(instance, getValueListName(responseClass))
            self.assertEqual(len(valueList), numValues)

    def testNextPageToken(self):
        responseClass = protocol.SearchVariantsResponse
        builder = response_builder.SearchResponseBuilder(
            responseClass, 100, 2 ** 32)
        # If not set, pageToken should be empty string
        self.assertIsNone(builder.getNextPageToken())
        instance = protocol.fromJson(builder.getSerializedResponse(),
                                     responseClass)
        self.assertEqual(instance.next_page_token, "")
        # page tokens can be any string.
        for nextPageToken in ["", "string"]:
            builder.setNextPageToken(nextPageToken)
            self.assertEqual(nextPageToken, builder.getNextPageToken())
            instance = protocol.fromJson(builder.getSerializedResponse(),
                                         responseClass)
            self.assertEqual(nextPageToken, instance.next_page_token)
