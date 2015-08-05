"""
Tests the compound ids
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel as datamodel
import ga4gh.exceptions as exceptions
import ga4gh.datamodel.variants as variants


class ExampleCompoundId(datamodel.CompoundId):
    separator = ';'
    fields = ['foo', 'bar', 'baz']
    comboFields = {
        'bazfoo': [2, 0],
    }


class TestCompoundIds(unittest.TestCase):
    """
    Test the compound ids
    """
    def testBadInit(self):
        with self.assertRaises(exceptions.BadIdentifierException):
            ExampleCompoundId(5)
        with self.assertRaises(exceptions.BadIdentifierException):
            ExampleCompoundId(None)
        with self.assertRaises(exceptions.BadIdentifierException):
            ExampleCompoundId('a;b')

    def testAttrs(self):
        compoundId = ExampleCompoundId('a;b;c')
        self.assertEqual(compoundId.foo, 'a')
        self.assertEqual(compoundId.bar, 'b')
        self.assertEqual(compoundId.baz, 'c')
        self.assertEqual(compoundId.bazfoo, 'c;a')

    def testCompose(self):
        compoundId = ExampleCompoundId.compose(bar=5, bazfoo='c;a')
        compoundIdStr = str(compoundId)
        self.assertEqual(compoundIdStr, 'a;5;c')
        self.assertEqual(compoundId.__class__, ExampleCompoundId)

    def testComposeFail(self):
        with self.assertRaises(ValueError):
            ExampleCompoundId.compose(notValid='invalidKey')
        with self.assertRaises(ValueError):
            ExampleCompoundId.compose(bazfoo='too;many;entries')

    def testCompoundVariantId(self):
        compoundId = variants.CompoundVariantId('a:b:c:d:e')
        self.assertEqual(compoundId.datasetId, 'a')
        self.assertEqual(compoundId.vsId, 'b')
        self.assertEqual(compoundId.referenceName, 'c')
        self.assertEqual(compoundId.start, 'd')
        self.assertEqual(compoundId.md5, 'e')
        self.assertEqual(compoundId.variantSetId, 'a:b')
        self.assertEqual(compoundId.variantId, 'c:d:e')
