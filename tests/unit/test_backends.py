"""
Tests for the backend objects. We instantiate local copies of
the backends and invoke the entry points for the protocol methods.
We do not set up any server processes or communicate over sockets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import unittest

import ga4gh.exceptions as exceptions
import ga4gh.backend as backend
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.references as references


class TestAbstractBackend(unittest.TestCase):
    """
    Tests for shared functionality between backends.
    """
    def setUp(self):
        self._backend = backend.AbstractBackend()

    def testAddOneDataset(self):
        datasetName = "ds"
        dataset = datasets.AbstractDataset(datasetName)
        self.assertEqual(self._backend.getNumDatasets(), 0)
        self.assertEqual(self._backend.getDatasets(), [])
        self._backend.addDataset(dataset)
        self.assertEqual(self._backend.getNumDatasets(), 1)
        self.assertEqual(self._backend.getDatasets(), [dataset])
        self.assertEqual(self._backend.getDatasetByIndex(0), dataset)
        self.assertEqual(self._backend.getDatasetByName(datasetName), dataset)
        self.assertEqual(self._backend.getDataset(dataset.getId()), dataset)

    def testAddMultipleDatasets(self):
        firstDatasetName = "ds1"
        firstDataset = datasets.AbstractDataset(firstDatasetName)
        secondDatasetName = "ds2"
        secondDataset = datasets.AbstractDataset(secondDatasetName)
        self.assertEqual(self._backend.getNumDatasets(), 0)
        self.assertEqual(self._backend.getDatasets(), [])
        self._backend.addDataset(firstDataset)
        self._backend.addDataset(secondDataset)
        self.assertEqual(self._backend.getNumDatasets(), 2)
        self.assertEqual(self._backend.getDatasets(),
                         [firstDataset, secondDataset])
        self.assertEqual(self._backend.getDatasetByIndex(0), firstDataset)
        self.assertEqual(self._backend.getDatasetByIndex(1), secondDataset)
        self.assertEqual(self._backend.getDatasetByName(firstDatasetName),
                         firstDataset)
        self.assertEqual(self._backend.getDatasetByName(secondDatasetName),
                         secondDataset)
        self.assertEqual(self._backend.getDataset(firstDataset.getId()),
                         firstDataset)
        self.assertEqual(self._backend.getDataset(secondDataset.getId()),
                         secondDataset)

    def testAddOneReferenceSet(self):
        referenceSetLocalId = "id"
        referenceSet = references.AbstractReferenceSet(referenceSetLocalId)
        self.assertEqual(self._backend.getNumReferenceSets(), 0)
        self.assertEqual(self._backend.getReferenceSets(), [])
        self._backend.addReferenceSet(referenceSet)
        self.assertEqual(self._backend.getNumReferenceSets(), 1)
        self.assertEqual(self._backend.getReferenceSets(), [referenceSet])
        self.assertEqual(self._backend.getReferenceSetByIndex(0), referenceSet)
        self.assertEqual(
            self._backend.getReferenceSetByName(referenceSet.getLocalId()),
            referenceSet)
        self.assertEqual(self._backend.getReferenceSet(referenceSet.getId()),
                         referenceSet)

    def testAddMultipleReferenceSet(self):
        firstRSLocalId = "id1"
        firstRS = references.AbstractReferenceSet(firstRSLocalId)
        secondRSLocalId = "id2"
        secondRS = references.AbstractReferenceSet(secondRSLocalId)
        self.assertEqual(self._backend.getNumReferenceSets(), 0)
        self.assertEqual(self._backend.getReferenceSets(), [])
        self._backend.addReferenceSet(firstRS)
        self._backend.addReferenceSet(secondRS)
        self.assertEqual(self._backend.getNumReferenceSets(), 2)
        self.assertEqual(self._backend.getReferenceSets(),
                         [firstRS, secondRS])
        self.assertEqual(self._backend.getReferenceSetByIndex(0),
                         firstRS)
        self.assertEqual(self._backend.getReferenceSetByIndex(1),
                         secondRS)
        self.assertEqual(
            self._backend.getReferenceSetByName(firstRS.getLocalId()),
            firstRS)
        self.assertEqual(
            self._backend.getReferenceSetByName(secondRS.getLocalId()),
            secondRS)
        self.assertEqual(self._backend.getReferenceSet(firstRS.getId()),
                         firstRS)
        self.assertEqual(self._backend.getReferenceSet(secondRS.getId()),
                         secondRS)

    def testGetDatasetBadId(self):
        for badId in ["", None, "NO SUCH ID"]:
            self.assertRaises(
                exceptions.DatasetNotFoundException,
                self._backend.getDataset, badId)

    def testGetReferenceSetBadId(self):
        for badId in ["", None, "NO SUCH ID"]:
            self.assertRaises(
                exceptions.ReferenceSetNotFoundException,
                self._backend.getReferenceSet, badId)

    def testGetDatasetBadName(self):
        for badName in ["", None, "NO SUCH NAME"]:
            self.assertRaises(
                exceptions.DatasetNameNotFoundException,
                self._backend.getDatasetByName, badName)

    def testGetReferenceSetBadName(self):
        for badName in ["", None, "NO SUCH NAME"]:
            self.assertRaises(
                exceptions.ReferenceSetNameNotFoundException,
                self._backend.getReferenceSetByName, badName)

    def testGetDatasetByIndexBadIndex(self):
        self.assertRaises(IndexError, self._backend.getDatasetByIndex, 0)
        self.assertRaises(TypeError, self._backend.getDatasetByIndex, None)
        self.assertRaises(TypeError, self._backend.getDatasetByIndex, "")
        datasetName = "ds"
        dataset = datasets.AbstractDataset(datasetName)
        self._backend.addDataset(dataset)
        self.assertRaises(IndexError, self._backend.getDatasetByIndex, 1)

    def testGetReferenceSetByIndexBadIndex(self):
        self.assertRaises(IndexError, self._backend.getReferenceSetByIndex, 0)
        self.assertRaises(TypeError,
                          self._backend.getReferenceSetByIndex, None)
        self.assertRaises(TypeError, self._backend.getReferenceSetByIndex, "")
        referenceSetName = "id"
        referenceSet = references.AbstractReferenceSet(referenceSetName)
        self._backend.addReferenceSet(referenceSet)
        self.assertRaises(IndexError, self._backend.getReferenceSetByIndex, 1)


class TestFileSystemBackend(unittest.TestCase):
    """
    Tests proper initialization of the filesystem backend using indexed
    files in the tests/data directory.
    """
    def setUp(self):
        self._dataDir = os.path.join("tests", "data")
        self._backend = backend.FileSystemBackend(self._dataDir)

    def testDatasets(self):
        self.assertEqual(self._backend.getNumDatasets(), 1)
        dataset = self._backend.getDatasetByIndex(0)
        self.assertEqual(dataset.getLocalId(), "dataset1")
        self.assertEqual(self._backend.getDatasetByName("dataset1"), dataset)

    def testReferenceSets(self):
        self.assertEqual(self._backend.getNumReferenceSets(), 4)
        referenceSets = enumerate(self._backend.getReferenceSets())
        referenceSetsByName = sorted(
            referenceSets, key=lambda x: x[1].getLocalId())
        expected_names = sorted(
            ["Default", "NCBI37", "example_1", "example_2"])
        for name, (index, rs) in zip(expected_names, referenceSetsByName):
            self.assertEqual(rs.getLocalId(), name)
            self.assertEqual(self._backend.getReferenceSetByIndex(index), rs)
            self.assertEqual(self._backend.getReferenceSet(rs.getId()), rs)
            self.assertEqual(self._backend.getReferenceSetByName(name), rs)


class TestTopLevelObjectGenerator(unittest.TestCase):
    """
    Tests the generator used for top level objects
    """
    def setUp(self):
        class FakeRequest(object):
            pass

        class FakeTopLevelObject(object):
            def toProtocolElement(self):
                return self

        self.request = FakeRequest()
        self.request.pageToken = None
        self.numObjects = 3
        self.objects = [FakeTopLevelObject() for j in range(self.numObjects)]
        self.backend = backend.AbstractBackend()

    def getObjectByIndex(self, index):
        return self.objects[index]

    def testPageToken(self):
        self.request.pageToken = "1"
        self._assertNumItems(2)

    def testPageTokenNone(self):
        self._assertNumItems(3)

    def _assertNumItems(self, numItems):
        iterator = self.backend._topLevelObjectGenerator(
            self.request, self.numObjects, self.getObjectByIndex)
        items = list(iterator)
        self.assertEqual(len(items), numItems)


class TestPrivateBackendMethods(unittest.TestCase):
    """
    keep tests of private backend methods here and not in one of the
    subclasses of TestAbstractBackend, otherwise the tests will needlessly
    be run more than once

    (they could be put in TestAbstractBackend, but I think it's a clearer
    separation to put them in their own test class)
    """
    def testParsePageToken(self):
        goodPageToken = "12:34:567:8:9000"
        parsedToken = backend._parsePageToken(goodPageToken, 5)
        self.assertEqual(parsedToken[2], 567)

    def testParseIntegerArgument(self):
        good = {"one": "1", "minusone": "-1"}
        expected = {"one": 1, "minusone": -1}
        bad = {"string": "A", "float": "0.98"}
        self.assertEqual(backend._parseIntegerArgument({}, "missing", 0), 0)
        for key in good:
            self.assertEqual(
                backend._parseIntegerArgument(good, key, 0), expected[key])
        for key in bad:
            with self.assertRaises(exceptions.BadRequestIntegerException):
                backend._parseIntegerArgument(bad, key, 0)
