"""
Tests the pysam sanitizer
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel as datamodel
import ga4gh.exceptions as exceptions


class TestPysamSanitizer(datamodel.PysamDatamodelMixin, unittest.TestCase):
    """
    Test the pysam sanitizer
    """
    def testSanitizeString(self):
        # only accept strings as params
        with self.assertRaises(exceptions.DatamodelValidationException):
            self.sanitizeString(1, 'longString')

        # return same string if string is short enough
        shortString = 'shortString'
        self.assertIsInstance(shortString, unicode)
        result = self.sanitizeString(shortString, shortString)
        self.assertEqual(shortString, result)
        self.assertNotIsInstance(result, unicode)

        # shorten string length if string too long
        longString = 'x' * (self.maxStringLength + 1)
        result = self.sanitizeString(longString, 'longString')
        self.assertEqual(longString[:self.maxStringLength], result)

    def testSanitizeInt(self):
        # only accept ints as params
        with self.assertRaises(exceptions.DatamodelValidationException):
            self.sanitizeInt('whoops', 0, 0, 0)

        # set to minVal if less than minVal
        minVal = 0
        result = self.sanitizeInt(-10, minVal, 100, 'int')
        self.assertEqual(minVal, result)

        # set to maxVal if more than maxVal
        maxVal = 100
        result = self.sanitizeInt(200, 0, maxVal, 'int')
        self.assertEqual(maxVal, result)

        # don't change if within range
        val = 50
        result = self.sanitizeInt(val, 0, 100, 'int')
        self.assertEqual(result, val)

    def testAssertValidRange(self):
        # wrong range should throw error
        with self.assertRaises(exceptions.DatamodelValidationException):
            self.assertValidRange(100, 0, 'start', 'end')

        # correct range should not throw error
        self.assertValidRange(0, 100, 'start', 'end')

    def testSanitizeVariantFileFetch(self):
        contigArg = 'x' * (self.maxStringLength + 1)
        startArg = self.vcfMin - 1
        stopArg = self.vcfMax + 1
        contig, start, stop = self.sanitizeVariantFileFetch(
            contigArg, startArg, stopArg)
        self.assertEqual(contigArg[:self.maxStringLength], contig)
        self.assertEqual(start, self.vcfMin)
        self.assertEqual(stop, self.vcfMax)

    def testSanitizeAlignmentFileFetch(self):
        referenceNameArg = 'x' * (self.maxStringLength + 1)
        startArg = self.samMin - 1
        endArg = self.samMaxEnd + 1
        referenceName, start, end = self.sanitizeAlignmentFileFetch(
            referenceNameArg, startArg, endArg)
        self.assertEqual(
            referenceNameArg[:self.maxStringLength], referenceName)
        self.assertEqual(start, self.samMin)
        self.assertEqual(end, self.samMaxEnd)
