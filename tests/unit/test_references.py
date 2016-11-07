"""
Unit tests for reference objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.server.backend as backend
import ga4gh.server.datamodel.references as references
import ga4gh.server.exceptions as exceptions
import ga4gh.server.datarepo as datarepo


class TestAbstractReferenceSet(unittest.TestCase):
    """
    Unit tests for the abstract reference set.
    """
    def setUp(self):
        self._backend = backend.Backend(datarepo.AbstractDataRepository())
        self._referenceSet = references.AbstractReferenceSet(
            'refSetId')

    def testAddOneReference(self):
        self.assertEqual(self._referenceSet.getNumReferences(), 0)
        referenceName = "ref"
        reference = references.AbstractReference(
            self._referenceSet, referenceName)
        self._referenceSet.addReference(reference)
        self.assertEqual(self._referenceSet.getNumReferences(), 1)
        self.assertEqual(
            self._referenceSet.getReferenceByIndex(0), reference)
        self.assertEqual(
            self._referenceSet.getReferenceByName(referenceName), reference)
        self.assertEqual(
            self._referenceSet.getReference(reference.getId()), reference)
        self.assertEqual(self._referenceSet.getReferences(), [reference])

    def testAddMultipleReference(self):
        referenceList = []
        referenceBaseName = "ref"
        referenceCount = 10

        self.assertEqual(self._referenceSet.getNumReferences(), 0)

        for i in range(referenceCount):
            referenceName = referenceBaseName + str(i)
            reference = references.AbstractReference(
                self._referenceSet, referenceName)
            referenceList.append(reference)
            self._referenceSet.addReference(reference)
            self.assertEqual(self._referenceSet.getNumReferences(), i + 1)
            self.assertEqual(
                self._referenceSet.getReferenceByIndex(i), reference)
            self.assertEqual(
                self._referenceSet.getReferenceByName(referenceName),
                reference)
            self.assertEqual(
                self._referenceSet.getReference(reference.getId()), reference)
            self.assertEqual(self._referenceSet.getReferences(), referenceList)

    def testReferenceNameNotFound(self):
        for badName in ["", None, "NO SUCH NAME"]:
            self.assertRaises(
                exceptions.ReferenceNameNotFoundException,
                self._referenceSet.getReferenceByName, badName)

    def testReferenceNotFound(self):
        for badId in ["", None, "NO SUCH ID"]:
            self.assertRaises(
                exceptions.ReferenceNotFoundException,
                self._referenceSet.getReference, badId)


class TestAbstractReference(unittest.TestCase):
    """
    Unit tests for the abstract reference object.
    """
    def setUp(self):
        self._backend = backend.Backend(
            datarepo.AbstractDataRepository())
        self._referenceSet = references.AbstractReferenceSet(
            'refSetId')
        self._reference = references.AbstractReference(
            self._referenceSet, "ref")

    def testReferenceRangeError(self):
        for badRange in [(-1, 3), (100, 50), (100, 1000)]:
            self.assertRaises(
                exceptions.ReferenceRangeErrorException,
                self._reference.checkQueryRange, badRange[0], badRange[1])
