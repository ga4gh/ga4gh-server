"""
Tests for the backend objects. We instantiate local copies of
the backends and invoke the entry points for the protocol methods.
We do not set up any server processes or communicate over sockets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.server.exceptions as exceptions
import ga4gh.server.backend as backend
import ga4gh.server.paging as paging
import ga4gh.server.datarepo as datarepo
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.datamodel.references as references

import tests.paths as paths


class TestAbstractBackend(unittest.TestCase):
    """
    Tests for shared functionality between backends.
    """
    def setUp(self):
        self._backend = backend.Backend(datarepo.AbstractDataRepository())
        self._dataRepo = self._backend.getDataRepository()

    def testAddOneDataset(self):
        datasetName = "ds"
        dataset = datasets.Dataset(datasetName)
        self.assertEqual(self._dataRepo.getNumDatasets(), 0)
        self.assertEqual(self._dataRepo.getDatasets(), [])
        self._dataRepo.addDataset(dataset)
        self.assertEqual(self._dataRepo.getNumDatasets(), 1)
        self.assertEqual(self._dataRepo.getDatasets(), [dataset])
        self.assertEqual(self._dataRepo.getDatasetByIndex(0), dataset)
        self.assertEqual(self._dataRepo.getDatasetByName(datasetName), dataset)
        self.assertEqual(self._dataRepo.getDataset(dataset.getId()), dataset)

    def testAddMultipleDatasets(self):
        firstDatasetName = "ds1"
        firstDataset = datasets.Dataset(firstDatasetName)
        secondDatasetName = "ds2"
        secondDataset = datasets.Dataset(secondDatasetName)
        self.assertEqual(self._dataRepo.getNumDatasets(), 0)
        self.assertEqual(self._dataRepo.getDatasets(), [])
        self._dataRepo.addDataset(firstDataset)
        self._dataRepo.addDataset(secondDataset)
        self.assertEqual(self._dataRepo.getNumDatasets(), 2)
        self.assertEqual(self._dataRepo.getDatasets(),
                         [firstDataset, secondDataset])
        self.assertEqual(self._dataRepo.getDatasetByIndex(0), firstDataset)
        self.assertEqual(self._dataRepo.getDatasetByIndex(1), secondDataset)
        self.assertEqual(self._dataRepo.getDatasetByName(firstDatasetName),
                         firstDataset)
        self.assertEqual(self._dataRepo.getDatasetByName(secondDatasetName),
                         secondDataset)
        self.assertEqual(self._dataRepo.getDataset(firstDataset.getId()),
                         firstDataset)
        self.assertEqual(self._dataRepo.getDataset(secondDataset.getId()),
                         secondDataset)

    def testAddOneReferenceSet(self):
        referenceSetLocalId = "id"
        referenceSet = references.AbstractReferenceSet(referenceSetLocalId)
        self.assertEqual(self._dataRepo.getNumReferenceSets(), 0)
        self.assertEqual(self._dataRepo.getReferenceSets(), [])
        self._dataRepo.addReferenceSet(referenceSet)
        self.assertEqual(self._dataRepo.getNumReferenceSets(), 1)
        self.assertEqual(self._dataRepo.getReferenceSets(), [referenceSet])
        self.assertEqual(
            self._dataRepo.getReferenceSetByIndex(0), referenceSet)
        self.assertEqual(
            self._dataRepo.getReferenceSetByName(referenceSet.getLocalId()),
            referenceSet)
        self.assertEqual(self._dataRepo.getReferenceSet(referenceSet.getId()),
                         referenceSet)

    def testAddMultipleReferenceSet(self):
        firstRSLocalId = "id1"
        firstRS = references.AbstractReferenceSet(firstRSLocalId)
        secondRSLocalId = "id2"
        secondRS = references.AbstractReferenceSet(secondRSLocalId)
        self.assertEqual(self._dataRepo.getNumReferenceSets(), 0)
        self.assertEqual(self._dataRepo.getReferenceSets(), [])
        self._dataRepo.addReferenceSet(firstRS)
        self._dataRepo.addReferenceSet(secondRS)
        self.assertEqual(self._dataRepo.getNumReferenceSets(), 2)
        self.assertEqual(self._dataRepo.getReferenceSets(),
                         [firstRS, secondRS])
        self.assertEqual(self._dataRepo.getReferenceSetByIndex(0),
                         firstRS)
        self.assertEqual(self._dataRepo.getReferenceSetByIndex(1),
                         secondRS)
        self.assertEqual(
            self._dataRepo.getReferenceSetByName(firstRS.getLocalId()),
            firstRS)
        self.assertEqual(
            self._dataRepo.getReferenceSetByName(secondRS.getLocalId()),
            secondRS)
        self.assertEqual(self._dataRepo.getReferenceSet(firstRS.getId()),
                         firstRS)
        self.assertEqual(self._dataRepo.getReferenceSet(secondRS.getId()),
                         secondRS)

    def testGetDatasetBadId(self):
        for badId in ["", None, "NO SUCH ID"]:
            self.assertRaises(
                exceptions.DatasetNotFoundException,
                self._dataRepo.getDataset, badId)

    def testGetReferenceSetBadId(self):
        for badId in ["", None, "NO SUCH ID"]:
            self.assertRaises(
                exceptions.ReferenceSetNotFoundException,
                self._dataRepo.getReferenceSet, badId)

    def testGetDatasetBadName(self):
        for badName in ["", None, "NO SUCH NAME"]:
            self.assertRaises(
                exceptions.DatasetNameNotFoundException,
                self._dataRepo.getDatasetByName, badName)

    def testGetReferenceSetBadName(self):
        for badName in ["", None, "NO SUCH NAME"]:
            self.assertRaises(
                exceptions.ReferenceSetNameNotFoundException,
                self._dataRepo.getReferenceSetByName, badName)

    def testGetDatasetByIndexBadIndex(self):
        self.assertRaises(IndexError, self._dataRepo.getDatasetByIndex, 0)
        self.assertRaises(TypeError, self._dataRepo.getDatasetByIndex, None)
        self.assertRaises(TypeError, self._dataRepo.getDatasetByIndex, "")
        datasetName = "ds"
        dataset = datasets.Dataset(datasetName)
        self._dataRepo.addDataset(dataset)
        self.assertRaises(IndexError, self._dataRepo.getDatasetByIndex, 1)

    def testGetReferenceSetByIndexBadIndex(self):
        self.assertRaises(IndexError, self._dataRepo.getReferenceSetByIndex, 0)
        self.assertRaises(TypeError,
                          self._dataRepo.getReferenceSetByIndex, None)
        self.assertRaises(TypeError, self._dataRepo.getReferenceSetByIndex, "")
        referenceSetName = "id"
        referenceSet = references.AbstractReferenceSet(referenceSetName)
        self._dataRepo.addReferenceSet(referenceSet)
        self.assertRaises(IndexError, self._dataRepo.getReferenceSetByIndex, 1)


class TestSqlRepoTestData(unittest.TestCase):
    """
    Tests proper initialization of the SQL repo based on known
    files in the tests/data directory.
    """
    def setUp(self):
        self._dataRepo = datarepo.SqlDataRepository(paths.testDataRepo)
        self._dataRepo.open(datarepo.MODE_READ)

    def testDatasets(self):
        self.assertEqual(self._dataRepo.getNumDatasets(), 1)
        dataset = self._dataRepo.getDatasetByIndex(0)
        self.assertEqual(dataset.getLocalId(), "dataset1")
        self.assertEqual(self._dataRepo.getDatasetByName("dataset1"), dataset)

    def testReferenceSets(self):
        self.assertEqual(self._dataRepo.getNumReferenceSets(), 4)
        referenceSets = enumerate(self._dataRepo.getReferenceSets())
        referenceSetsByName = sorted(
            referenceSets, key=lambda x: x[1].getLocalId())
        expected_names = sorted(
            ["Default", "NCBI37", "example_1", "example_2"])
        for name, (index, rs) in zip(expected_names, referenceSetsByName):
            self.assertEqual(rs.getLocalId(), name)
            self.assertEqual(self._dataRepo.getReferenceSetByIndex(index), rs)
            self.assertEqual(self._dataRepo.getReferenceSet(rs.getId()), rs)
            self.assertEqual(self._dataRepo.getReferenceSetByName(name), rs)


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
        self.request.page_token = ""
        self.num_objects = 3
        self.objects = [FakeTopLevelObject() for j in range(self.num_objects)]
        self.backend = backend.Backend(datarepo.AbstractDataRepository())

    def getObjectByIndex(self, index):
        return self.objects[index]

    def testPageToken(self):
        self.request.page_token = "1"
        self._assertNumItems(2)

    def testPageTokenNone(self):
        self._assertNumItems(3)

    def _assertNumItems(self, numItems):
        iterator = self.backend._topLevelObjectGenerator(
            self.request, self.num_objects, self.getObjectByIndex)
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
        parsedToken = paging._parsePageToken(goodPageToken, 5)
        self.assertEqual(parsedToken[2], 567)

    def testParseIntegerArgument(self):
        good = {"one": "1", "minusone": "-1"}
        expected = {"one": 1, "minusone": -1}
        bad = {"string": "A", "float": "0.98"}
        self.assertEqual(paging._parseIntegerArgument(
            {}, "missing", 0), 0)
        for key in good:
            self.assertEqual(
                paging._parseIntegerArgument(
                    good, key, 0), expected[key])
        for key in bad:
            with self.assertRaises(exceptions.BadRequestIntegerException):
                paging._parseIntegerArgument(bad, key, 0)
