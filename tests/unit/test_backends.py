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
            self, request, page_size, searchMethod, ResponseClass, listMember):
        """
        Returns an iterator over the list of results from the specified
        request.  All results are returned, and paging is handled
        automatically.
        """
        notDone = True
        request.page_size = page_size
        while notDone:
            # TODO validate the response there.
            responseStr = searchMethod(request.toJsonString())
            response = ResponseClass.fromJsonString(responseStr)
            objectList = getattr(response, listMember)
            self.assertLessEqual(len(objectList), page_size)
            for obj in objectList:
                yield obj
            notDone = response.next_page_token is not None
            request.page_token = response.next_page_token

    def getVariantSets(self, page_size=100):
        """
        Returns an iterator over the variantSets, abstracting away
        the details of the page_size.
        """
        request = protocol.SearchVariantSetsRequest()
        request.dataset_ids = [self._backend.getDatasetIds()[0]]
        return self.resultIterator(
            request, page_size, self._backend.searchVariantSets,
            protocol.SearchVariantSetsResponse, "variant_sets")

    def getVariants(
            self, variant_set_ids, reference_name, start=0, end=2 ** 32,
            page_size=100, call_set_ids=None):
        """
        Returns an iterator over the specified list of variants,
        abstracting out paging details.
        """
        request = protocol.SearchVariantsRequest()
        request.variant_set_ids = variant_set_ids
        request.reference_name = reference_name
        request.start = start
        request.end = end
        request.call_set_ids = call_set_ids
        return self.resultIterator(
            request, page_size, self._backend.searchVariants,
            protocol.SearchVariantsResponse, "variants")

    def getCallSets(self, variant_set_id, page_size=100):
        """
        Returns an iterator over the callsets in a specified
        variant set.
        """
        request = protocol.SearchCallSetsRequest()
        request.variant_set_ids = [variant_set_id]
        return self.resultIterator(
            request, page_size, self._backend.searchCallSets,
            protocol.SearchCallSetsResponse, "call_sets")

    def testGetVariantSets(self):
        dataset_id = self._backend.getDatasetIds()[0]
        sortedVariantSetsFromGetter = sorted(
            self._backend.getDataset(dataset_id).getVariantSets())
        sortedVariantSetMapValues = sorted(
            self._backend.getDataset(dataset_id)._variantSetIdMap.values())
        self.assertEqual(
            sortedVariantSetMapValues, sortedVariantSetsFromGetter)

    def testParsePageToken(self):
        goodPageToken = "12:34:567:8:9000"
        parsedToken = backend._parsePageToken(goodPageToken, 5)
        self.assertEqual(parsedToken[2], 567)

    def testRunSearchRequest(self):
        request = protocol.SearchVariantSetsRequest()
        request.dataset_ids = [self._backend.getDatasetIds()[0]]
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
        id_ = "referenceSet0:srs0"
        self.runListReferenceBases(id_)

    def testSearchVariantSets(self):
        request = protocol.SearchVariantSetsRequest()
        request.dataset_ids = [self._backend.getDatasetIds()[0]]
        responseStr = self._backend.searchVariantSets(request.toJsonString())
        response = protocol.SearchVariantSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchVariantSetsResponse))

    def testSearchVariants(self):
        variant_set_ids = [
            variantSet.id for variantSet in self.getVariantSets(page_size=1)]
        request = protocol.SearchVariantsRequest()
        request.variant_set_ids = variant_set_ids[:1]
        responseStr = self._backend.searchVariants(request.toJsonString())
        response = protocol.SearchVariantsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchVariantsResponse))

    def testSearchCallSets(self):
        variant_set_ids = [
            variantSet.id for variantSet in self.getVariantSets(page_size=1)]
        request = protocol.SearchCallSetsRequest()
        request.variant_set_ids = variant_set_ids[:1]
        responseStr = self._backend.searchCallSets(request.toJsonString())
        response = protocol.SearchCallSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.SearchCallSetsResponse))

    def testVariantSetPagination(self):
        results = []
        for page_size in range(1, 100):
            variant_set_ids = [
                variantSet.id for variantSet in self.getVariantSets(
                    page_size=page_size)]
            results.append(variant_set_ids)
        for result in results[1:]:
            self.assertEqual(result, results[0])

    def runListReferenceBases(self, id_):
        requestArgs = {"start": 5, "end": 10, "page_token": "0"}
        responseStr = self._backend.listReferenceBases(id_, requestArgs)
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
        variant_sets = [variantSet for variantSet in self.getVariantSets()]
        self.assertEqual(len(variant_sets), len(self._vcfs))
        ids = set(variantSet.id for variantSet in variant_sets)
        self.assertEqual(ids, set(self._vcfs.keys()))

    def testRunListReferenceBases(self):
        id_ = "example_1:simple"
        self.runListReferenceBases(id_)

    def testOneDatasetRestriction(self):
        # no dataset_ids attr
        request = protocol.SearchReadsRequest()
        self._backend._getDatasetFromRequest(request)

        # dataset_ids attr
        request = protocol.SearchVariantSetsRequest()
        with self.assertRaises(exceptions.NotExactlyOneDatasetException):
            self._backend._getDatasetFromRequest(request)
        dataset_id = 'dataset1'
        request.dataset_ids = [dataset_id]
        dataset = self._backend._getDatasetFromRequest(request)
        self.assertEquals(dataset.getId(), dataset_id)
        request.dataset_ids = ['dataset1', 'dataset2']
        with self.assertRaises(exceptions.NotExactlyOneDatasetException):
            self._backend._getDatasetFromRequest(request)


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
        self.request.page_token = None
        self.idMap = {
            "a": FakeTopLevelObject(),
            "b": FakeTopLevelObject(),
            "c": FakeTopLevelObject(),
        }
        self.idList = sorted(self.idMap.keys())
        self.backend = backend.AbstractBackend()

    def testPageToken(self):
        self.request.page_token = "1"
        self._assertNumItems(2)

    def testPageTokenNone(self):
        self._assertNumItems(3)

    def _assertNumItems(self, numItems):
        iterator = self.backend._topLevelObjectGenerator(
            self.request, self.idMap, self.idList)
        items = list(iterator)
        self.assertEqual(len(items), numItems)
