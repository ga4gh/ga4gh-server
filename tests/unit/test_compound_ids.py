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
import ga4gh.datamodel.references as references
import ga4gh.datamodel.reads as reads


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
        compoundId = variants.CompoundVariantId('a:b:c:d')
        self.assertEqual(compoundId.datasetId, 'a')
        self.assertEqual(compoundId.vsId, 'b')
        self.assertEqual(compoundId.referenceName, 'c')
        self.assertEqual(compoundId.start, 'd')
        self.assertEqual(compoundId.variantSetId, 'a:b')
        self.assertEqual(compoundId.variantId, 'c:d')

    def testCompoundVariantSetId(self):
        compoundId = variants.CompoundVariantSetId('a:b')
        self.assertEqual(compoundId.datasetId, 'a')
        self.assertEqual(compoundId.vsId, 'b')

    def testCompoundCallsetId(self):
        compoundId = variants.CompoundCallsetId('a:b:c')
        self.assertEqual(compoundId.datasetId, 'a')
        self.assertEqual(compoundId.vsId, 'b')
        self.assertEqual(compoundId.csId, 'c')

    def testCompoundReadGroupSetId(self):
        compoundId = reads.CompoundReadGroupSetId('a:b')
        self.assertEqual(compoundId.datasetId, 'a')
        self.assertEqual(compoundId.rgsId, 'b')

    def testCompoundReadGroupId(self):
        compoundId = reads.CompoundReadGroupId('a:b:c')
        self.assertEqual(compoundId.datasetId, 'a')
        self.assertEqual(compoundId.rgsId, 'b')
        self.assertEqual(compoundId.rgId, 'c')
        self.assertEqual(compoundId.readGroupSetId, 'a:b')

    def testCompoundReadAlignmentId(self):
        compoundId = reads.CompoundReadAlignmentId('a:b:c:d')
        self.assertEqual(compoundId.datasetId, 'a')
        self.assertEqual(compoundId.rgsId, 'b')
        self.assertEqual(compoundId.rgId, 'c')
        self.assertEqual(compoundId.raId, 'd')
        self.assertEqual(compoundId.readGroupSetId, 'a:b')
        self.assertEqual(compoundId.readGroupId, 'a:b:c')

    def testCompoundReferenceId(self):
        compoundId = references.CompoundReferenceId('a:b')
        self.assertEqual(compoundId.referenceSetId, 'a')
        self.assertEqual(compoundId.rId, 'b')
