"""
Tests for the backend objects. We instantiate local copies of
the backends and invoke the entry points for the protocol methods.
We do not set up any server processes or communicate over sockets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import unittest

import pysam

import ga4gh.exceptions as exceptions
import ga4gh.backend as backend
import ga4gh.protocol as protocol
import ga4gh.datamodel.references as references
import ga4gh.datamodel.variants as variants


class BackendForTesting(backend.AbstractBackend):
    """
    A backend to test abstract methods
    """


class TestAbstractBackend(unittest.TestCase):
    """
    Provides testing harness for testing methods in AbstractBackend,
    using an instance of the mock SimulatedBackend object.
    """
    def setUp(self):
        self._backend = backend.SimulatedBackend(
            numCalls=100, numVariantSets=10)
        # TODO arbitrary values, pepper to taste

    def getDataset(self):
        return self._backend.getDatasets()[0]

    def resultIterator(
            self, request, pageSize, searchMethod, ResponseClass, listMember):
        """
        Returns an iterator over the list of results from the specified
        request.  All results are returned, and paging is handled
        automatically.
        """
        notDone = True
        request.pageSize = pageSize
        while notDone:
            # TODO validate the response there.
            responseStr = searchMethod(request.toJsonString())
            response = ResponseClass.fromJsonString(responseStr)
            objectList = getattr(response, listMember)
            self.assertLessEqual(len(objectList), pageSize)
            for obj in objectList:
                yield obj
            notDone = response.nextPageToken is not None
            request.pageToken = response.nextPageToken

    def getVariantSets(self, pageSize=100):
        """
        Returns an iterator over the variantSets, abstracting away
        the details of the pageSize.
        """
        request = protocol.SearchVariantSetsRequest()
        request.datasetId = self.getDataset().getId()
        return self.resultIterator(
            request, pageSize, self._backend.runSearchVariantSets,
            protocol.SearchVariantSetsResponse, "variantSets")

    def getVariants(
            self, variantSetIds, referenceName, start=0, end=2 ** 32,
            pageSize=100, callSetIds=None):
        """
        Returns an iterator over the specified list of variants,
        abstracting out paging details.
        """
        request = protocol.SearchVariantsRequest()
        request.variantSetIds = variantSetIds
        request.referenceName = referenceName
        request.start = start
        request.end = end
        request.callSetIds = callSetIds
        return self.resultIterator(
            request, pageSize, self._backend.runSearchVariants,
            protocol.SearchVariantsResponse, "variants")

    def getCallSets(self, variantSetId, pageSize=100):
        """
        Returns an iterator over the callsets in a specified
        variant set.
        """
        request = protocol.SearchCallSetsRequest()
        request.variantSetIds = [variantSetId]
        return self.resultIterator(
            request, pageSize, self._backend.runSearchCallSets,
            protocol.SearchCallSetsResponse, "callSets")

    def testGetVariantSets(self):
        datasetId = self.getDataset().getId()
        sortedVariantSetsFromGetter = sorted(
            self._backend.getDataset(datasetId).getVariantSets())
        sortedVariantSetMapValues = sorted(
            self._backend.getDataset(datasetId)._variantSetIdMap.values())
        self.assertEqual(
            sortedVariantSetMapValues, sortedVariantSetsFromGetter)

    def testRunSearchRequest(self):
        request = protocol.SearchVariantSetsRequest()
        request.datasetId = self.getDataset().getId()
        responseStr = self._backend.runSearchRequest(
            request.toJsonString(), protocol.SearchVariantSetsRequest,
            protocol.SearchVariantSetsResponse,
            self._backend.variantSetsGenerator)
        response = protocol.SearchVariantSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchVariantSetsResponse))

    def testRunGetRequest(self):
        id_ = "anId"
        obj = references.SimulatedReferenceSet(id_)
        idMap = {id_: obj}
        responseStr = self._backend.runGetRequest(idMap, id_)
        class_ = protocol.ReferenceSet
        response = class_.fromJsonString(responseStr)
        self.assertTrue(isinstance(response, class_))

    def testRunListReferenceBases(self):
        id_ = self._backend.getReferenceSets()[0].getReferences()[0].getId()
        self.runListReferenceBases(id_)

    def testSearchVariantSets(self):
        request = protocol.SearchVariantSetsRequest()
        request.datasetId = self.getDataset().getId()
        responseStr = self._backend.runSearchVariantSets(
            request.toJsonString())
        response = protocol.SearchVariantSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchVariantSetsResponse))

    def testSearchVariants(self):
        variantSetIds = [
            variantSet.id for variantSet in self.getVariantSets(pageSize=1)]
        request = protocol.SearchVariantsRequest()
        request.variantSetId = variantSetIds[0]
        responseStr = self._backend.runSearchVariants(request.toJsonString())
        response = protocol.SearchVariantsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchVariantsResponse))

    def testSearchCallSets(self):
        variantSetIds = [
            variantSet.id for variantSet in self.getVariantSets(pageSize=1)]
        request = protocol.SearchCallSetsRequest()
        request.variantSetId = variantSetIds[0]
        responseStr = self._backend.runSearchCallSets(request.toJsonString())
        response = protocol.SearchCallSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchCallSetsResponse))

    def testVariantSetPagination(self):
        results = []
        for pageSize in range(1, 100):
            variantSetIds = [
                variantSet.id for variantSet in self.getVariantSets(
                    pageSize=pageSize)]
            results.append(variantSetIds)
        for result in results[1:]:
            self.assertEqual(result, results[0])

    def runListReferenceBases(self, id_):
        requestArgs = {"start": 5, "end": 10, "pageToken": "0"}
        responseStr = self._backend.runListReferenceBases(id_, requestArgs)
        response = protocol.ListReferenceBasesResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.ListReferenceBasesResponse))


class TestFileSystemBackend(TestAbstractBackend):
    """
    Tests proper initialization of the filesystem backend using indexed
    files in the tests/data directory.
    """
    def setUp(self):
        self._dataDir = os.path.join("tests", "data")
        self._referencesDir = os.path.join(self._dataDir, "references")
        self._datasetDir = os.path.join(self._dataDir, "dataset1")
        self._variantsDir = os.path.join(self._datasetDir, "variants")
        self._vcfs = {}
        self._variants = []
        self._referenceNames = set()
        self._chromFileMap = {}
        for relativePath in os.listdir(self._variantsDir):
            pathToFiles = os.path.join(self._variantsDir, relativePath)
            self._vcfs[relativePath] = []
            for vcfFile in glob.glob(os.path.join(
                    pathToFiles, "*.vcf.gz")):
                self._chromFileMap[relativePath] = {}
                self._vcfs[relativePath].append(vcfFile)
                vcf = pysam.VariantFile(filename=vcfFile)
                for chrom in vcf.index:
                    self._chromFileMap[relativePath][chrom] = vcf
        self._backend = backend.FileSystemBackend(self._dataDir)

    def testVariantSetIds(self):
        variantSets = [variantSet for variantSet in self.getVariantSets()]
        self.assertEqual(len(variantSets), len(self._vcfs))
        ids = set(variantSet.id for variantSet in variantSets)
        dataset = self.getDataset()
        vcfKeys = set()
        for localId in self._vcfs.keys():
            tmpVariantSet = variants.AbstractVariantSet(dataset, localId)
            vcfKeys.add(tmpVariantSet.getId())
        self.assertEqual(ids, vcfKeys)

    def testRunListReferenceBases(self):
        id_ = "example_1:simple"
        self.runListReferenceBases(id_)

    def testDatasetNotFound(self):
        request = protocol.SearchVariantSetsRequest()
        datasetId = 'doesNotExist'
        request.datasetId = datasetId
        with self.assertRaises(exceptions.DatasetNotFoundException):
            self._backend.getDataset(request.datasetId)


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
        self.idMap = {
            "a": FakeTopLevelObject(),
            "b": FakeTopLevelObject(),
            "c": FakeTopLevelObject(),
        }
        self.idList = sorted(self.idMap.keys())
        self.backend = backend.AbstractBackend()

    def testPageToken(self):
        self.request.pageToken = "1"
        self._assertNumItems(2)

    def testPageTokenNone(self):
        self._assertNumItems(3)

    def _assertNumItems(self, numItems):
        iterator = self.backend._topLevelObjectGenerator(
            self.request, self.idMap, self.idList)
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

    def testSafeMapQuery(self):
        # test map
        key = 'a'
        value = 'b'
        idMap = {key: value}
        result = backend._safeMapQuery(idMap, key)
        self.assertEqual(result, value)

        # test array
        arr = [value]
        result = backend._safeMapQuery(arr, 0)
        self.assertEqual(result, value)

        # test exception with custom class
        exceptionClass = exceptions.FileOpenFailedException
        with self.assertRaises(exceptionClass):
            backend._safeMapQuery(idMap, 'notFound', exceptionClass)

        # test exception with custom id string
        with self.assertRaises(exceptions.ObjectWithIdNotFoundException):
            backend._safeMapQuery(idMap, 'notFound', idErrorString='msg')
