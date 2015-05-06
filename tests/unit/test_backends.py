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

import ga4gh.backend as backend
import ga4gh.protocol as protocol


class TestAbstractBackend(unittest.TestCase):
    """
    Provides testing harness for testing methods in AbstractBackend,
    using an instance of the mock SimulatedBackend object.
    """
    def setUp(self):
        self._backend = backend.SimulatedBackend(
            numCalls=100, numVariantSets=10)
        # TODO arbitrary values, pepper to taste

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
        return self.resultIterator(
            request, pageSize, self._backend.searchVariantSets,
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
            request, pageSize, self._backend.searchVariants,
            protocol.SearchVariantsResponse, "variants")

    def getCallSets(self, variantSetId, pageSize=100):
        """
        Returns an iterator over the callsets in a specified
        variant set.
        """
        request = protocol.SearchCallSetsRequest()
        request.variantSetIds = [variantSetId]
        return self.resultIterator(
            request, pageSize, self._backend.searchCallSets,
            protocol.SearchCallSetsResponse, "callSets")

    def testGetVariantSets(self):
        sortedVariantSetsFromGetter = sorted(self._backend.getVariantSets())
        sortedVariantSetMapValues = sorted(
            self._backend._variantSetIdMap.values())
        self.assertEqual(
            sortedVariantSetMapValues, sortedVariantSetsFromGetter)

    def testParsePageToken(self):
        goodPageToken = "12:34:567:8:9000"
        parsedToken = backend._parsePageToken(goodPageToken, 5)
        self.assertEqual(parsedToken[2], 567)

    def testRunSearchRequest(self):
        request = protocol.SearchVariantSetsRequest()
        responseStr = self._backend.runSearchRequest(
            request.toJsonString(), protocol.SearchVariantSetsRequest,
            protocol.SearchVariantSetsResponse,
            self._backend.variantSetsGenerator)
        response = protocol.SearchVariantSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchVariantSetsResponse))

    def testSearchVariantSets(self):
        request = protocol.SearchVariantSetsRequest()
        responseStr = self._backend.searchVariantSets(request.toJsonString())
        response = protocol.SearchVariantSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchVariantSetsResponse))

    def testSearchVariants(self):
        variantSetIds = [
            variantSet.id for variantSet in self.getVariantSets(pageSize=1)]
        request = protocol.SearchVariantsRequest()
        request.variantSetIds = variantSetIds[:1]
        responseStr = self._backend.searchVariants(request.toJsonString())
        response = protocol.SearchVariantsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchVariantsResponse))

    def testSearchCallSets(self):
        variantSetIds = [
            variantSet.id for variantSet in self.getVariantSets(pageSize=1)]
        request = protocol.SearchCallSetsRequest()
        request.variantSetIds = variantSetIds[:1]
        responseStr = self._backend.searchCallSets(request.toJsonString())
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


class TestFileSystemBackend(TestAbstractBackend):
    """
    Tests proper initialization of the filesystem backend using indexed
    files in the tests/data directory.
    """
    def setUp(self):
        self._dataDir = os.path.join("tests", "data")
        self._variantsDir = os.path.join(self._dataDir, "variants")
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
        self.assertEqual(ids, set(self._vcfs.keys()))


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
