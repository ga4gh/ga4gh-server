"""
Tests the compound ids
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel as datamodel
import ga4gh.exceptions as exceptions
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.variants as variants
import ga4gh.datamodel.references as references
import ga4gh.datamodel.reads as reads


class ExampleCompoundId(datamodel.CompoundId):
    separator = ';'
    fields = ['foo', 'bar', 'baz']
    containerIds = [('foobar', 1), ('foobarbaz', 2)]


class TestCompoundIds(unittest.TestCase):
    """
    Test the compound ids
    """
    def testBadParse(self):
        for badId in [5, None, 'a;b', 'a;b;c;d', 'a;b;sd;', ';;;;']:
            with self.assertRaises(exceptions.BadIdentifierException):
                ExampleCompoundId.parse(badId)

    def verifyParseFailure(self, idStr, compoundIdClass):
        """
        Verifies that substrings and superstrings of the specified parsing
        ID string correctly raise parse failures.
        """
        # first, check if we really can parse the string
        cid = compoundIdClass.parse(idStr)
        self.assertIsNotNone(cid)
        # Now, check for substrings
        for j in range(len(idStr) - 1):
            self.assertRaises(
                exceptions.BadIdentifierException, compoundIdClass.parse,
                idStr[:j])
        # Adding on an extra field should also provoke a parse error.
        self.assertRaises(
            exceptions.BadIdentifierException, compoundIdClass.parse,
            idStr + ":b")

    def testAttrs(self):
        compoundId = ExampleCompoundId.parse('a;b;c')
        self.assertEqual(compoundId.foo, 'a')
        self.assertEqual(compoundId.bar, 'b')
        self.assertEqual(compoundId.baz, 'c')
        self.assertEqual(compoundId.foobar, 'a;b')
        self.assertEqual(compoundId.foobarbaz, 'a;b;c')

    def testInstantiate(self):
        compoundId = ExampleCompoundId(None, "a", "5", "c")
        compoundIdStr = str(compoundId)
        self.assertEqual(compoundIdStr, 'a;5;c')
        self.assertEqual(compoundId.__class__, ExampleCompoundId)

    def getDataset(self):
        return datasets.AbstractDataset("dataset")

    def getReferenceSet(self):
        return references.AbstractReferenceSet("referenceSet")

    def getVariantSet(self):
        return variants.AbstractVariantSet(self.getDataset(), "variantSet")

    def getReadGroupSet(self):
        return reads.AbstractReadGroupSet(self.getDataset(), "readGroupSet")

    def getReadGroup(self):
        return reads.AbstractReadGroup(self.getReadGroupSet(), "readGroup")

    def testDataset(self):
        localId = "dataset"
        cid = datamodel.DatasetCompoundId(None, localId)
        self.assertRaises(
            ValueError, datamodel.DatasetCompoundId, None)
        self.assertEqual(cid.dataset, localId)

    def testDatasetParse(self):
        idStr = "a"
        cid = datamodel.DatasetCompoundId.parse(idStr)
        self.assertEqual(cid.dataset, "a")
        self.verifyParseFailure(idStr, datamodel.DatasetCompoundId)

    def testVariantSet(self):
        dataset = self.getDataset()
        localId = "variantSet"
        cid = datamodel.VariantSetCompoundId(dataset.getCompoundId(), localId)
        self.assertRaises(
            ValueError, datamodel.VariantCompoundId, dataset.getCompoundId())
        self.assertEqual(cid.dataset, dataset.getLocalId())
        self.assertEqual(cid.variantSet, localId)
        self.assertEqual(cid.datasetId, dataset.getId())

    def testVariantSetParse(self):
        idStr = "a:b"
        cid = datamodel.VariantSetCompoundId.parse(idStr)
        self.assertEqual(cid.dataset, "a")
        self.assertEqual(cid.variantSet, "b")
        self.verifyParseFailure(idStr, datamodel.VariantSetCompoundId)

    def testCallSet(self):
        name = "sampleName"
        variantSet = self.getVariantSet()
        dataset = variantSet.getParentContainer()
        cid = datamodel.CallSetCompoundId(variantSet.getCompoundId(), name)
        self.assertRaises(
            ValueError, datamodel.CallSetCompoundId,
            variantSet.getCompoundId())
        self.assertEqual(cid.dataset, dataset.getLocalId())
        self.assertEqual(cid.variantSet, variantSet.getLocalId())
        self.assertEqual(cid.name, name)
        self.assertEqual(cid.datasetId, dataset.getId())
        self.assertEqual(cid.variantSetId, variantSet.getId())

    def testCallSetParse(self):
        idStr = "a:b:c"
        cid = datamodel.CallSetCompoundId.parse(idStr)
        self.assertEqual(cid.dataset, "a")
        self.assertEqual(cid.variantSet, "b")
        self.assertEqual(cid.name, "c")
        self.verifyParseFailure(idStr, datamodel.CallSetCompoundId)

    def testVariant(self):
        referenceName = "referenceName"
        start = "start"
        variantSet = self.getVariantSet()
        dataset = variantSet.getParentContainer()
        cid = datamodel.VariantCompoundId(
            variantSet.getCompoundId(), referenceName, start)
        self.assertRaises(
            ValueError, datamodel.VariantCompoundId,
            variantSet.getCompoundId())
        self.assertRaises(
            ValueError, datamodel.VariantCompoundId,
            variantSet.getCompoundId(), referenceName)
        self.assertEqual(cid.dataset, dataset.getLocalId())
        self.assertEqual(cid.variantSet, variantSet.getLocalId())
        self.assertEqual(cid.referenceName, referenceName)
        self.assertEqual(cid.start, start)
        self.assertEqual(cid.datasetId, dataset.getId())
        self.assertEqual(cid.variantSetId, variantSet.getId())

    def testVariantParse(self):
        idStr = "a:b:c:d"
        cid = datamodel.VariantCompoundId.parse(idStr)
        self.assertEqual(cid.dataset, "a")
        self.assertEqual(cid.variantSet, "b")
        self.assertEqual(cid.referenceName, "c")
        self.assertEqual(cid.start, "d")
        self.verifyParseFailure(idStr, datamodel.VariantCompoundId)

    def testReferenceSet(self):
        localId = "referenceSet"
        cid = datamodel.ReferenceSetCompoundId(None, localId)
        self.assertRaises(
            ValueError, datamodel.ReferenceSetCompoundId, None)
        self.assertEqual(cid.referenceSet, localId)

    def testReferenceSetParse(self):
        idStr = "a"
        cid = datamodel.ReferenceSetCompoundId.parse(idStr)
        self.assertEqual(cid.referenceSet, "a")
        self.verifyParseFailure(idStr, datamodel.ReferenceSetCompoundId)

    def testReference(self):
        referenceSet = self.getReferenceSet()
        localId = "reference"
        cid = datamodel.ReferenceCompoundId(
            referenceSet.getCompoundId(), localId)
        self.assertRaises(
            ValueError, datamodel.ReferenceCompoundId,
            referenceSet.getCompoundId())
        self.assertEqual(cid.referenceSet, referenceSet.getLocalId())
        self.assertEqual(cid.reference, localId)
        self.assertEqual(cid.referenceSetId, referenceSet.getId())

    def testReferenceParse(self):
        idStr = "a:b"
        cid = datamodel.ReferenceCompoundId.parse(idStr)
        self.assertEqual(cid.referenceSet, "a")
        self.assertEqual(cid.reference, "b")
        self.verifyParseFailure(idStr, datamodel.ReferenceCompoundId)

    def testReadGroupSet(self):
        dataset = self.getDataset()
        localId = "readGroupSet"
        cid = datamodel.ReadGroupSetCompoundId(
            dataset.getCompoundId(), localId)
        self.assertRaises(
            ValueError, datamodel.ReadGroupSetCompoundId,
            dataset.getCompoundId())
        self.assertEqual(cid.dataset, dataset.getLocalId())
        self.assertEqual(cid.readGroupSet, localId)

    def testReadGroupSetParse(self):
        idStr = "a:b"
        cid = datamodel.ReadGroupSetCompoundId.parse(idStr)
        self.assertEqual(cid.dataset, "a")
        self.assertEqual(cid.readGroupSet, "b")
        self.verifyParseFailure(idStr, datamodel.ReadGroupSetCompoundId)

    def testReadGroup(self):
        readGroupSet = self.getReadGroupSet()
        dataset = readGroupSet.getParentContainer()
        localId = "readGroup"
        cid = datamodel.ReadGroupCompoundId(
            readGroupSet.getCompoundId(), localId)
        self.assertRaises(
            ValueError, datamodel.ReadGroupCompoundId,
            readGroupSet.getCompoundId())
        self.assertEqual(cid.dataset, dataset.getLocalId())
        self.assertEqual(cid.readGroupSet, readGroupSet.getLocalId())
        self.assertEqual(cid.readGroup, localId)
        self.assertEqual(cid.datasetId, dataset.getId())
        self.assertEqual(cid.readGroupSetId, readGroupSet.getId())

    def testReadGroupParse(self):
        idStr = "a:b:c"
        cid = datamodel.ReadGroupCompoundId.parse(idStr)
        self.assertEqual(cid.dataset, "a")
        self.assertEqual(cid.readGroupSet, "b")
        self.assertEqual(cid.readGroup, "c")
        self.verifyParseFailure(idStr, datamodel.ReadGroupCompoundId)

    def testReadAlignment(self):
        readGroup = self.getReadGroup()
        readGroupSet = readGroup.getParentContainer()
        dataset = readGroupSet.getParentContainer()
        localId = "alignment"
        cid = datamodel.ReadAlignmentCompoundId(
            readGroup.getCompoundId(), localId)
        self.assertRaises(
            ValueError, datamodel.ReadAlignmentCompoundId,
            dataset.getCompoundId())
        self.assertEqual(cid.dataset, dataset.getLocalId())
        self.assertEqual(cid.readGroupSet, readGroupSet.getLocalId())
        self.assertEqual(cid.readGroup, readGroup.getLocalId())
        self.assertEqual(cid.readAlignment, localId)
        self.assertEqual(cid.datasetId, dataset.getId())
        self.assertEqual(cid.readGroupSetId, readGroupSet.getId())
        self.assertEqual(cid.readGroupId, readGroup.getId())

    def testReadAlignmentParse(self):
        idStr = "a:b:c:d"
        cid = datamodel.ReadAlignmentCompoundId.parse(idStr)
        self.assertEqual(cid.dataset, "a")
        self.assertEqual(cid.readGroupSet, "b")
        self.assertEqual(cid.readGroup, "c")
        self.assertEqual(cid.readAlignment, "d")
        self.verifyParseFailure(idStr, datamodel.ReadAlignmentCompoundId)
