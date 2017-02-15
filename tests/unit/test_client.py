"""
Tests for the client
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import json

import ga4gh.server.backend as backend
import ga4gh.server.datarepo as datarepo

import ga4gh.client.client as client
import ga4gh.common.utils as utils
import ga4gh.schemas.protocol as protocol


class DatamodelObjectWrapper(object):
    """
    Thin wrapper class that allows us to treat data model objects uniformly.
    We should update the data model interface so that objects are always
    returned so that we always call toProtocolElement on the results.
    Variants and Reads are the exceptions here.
    """
    def __init__(self, gaObject):
        self.gaObject = gaObject

    def toProtocolElement(self):
        return self.gaObject


class DummyResponse(object):
    """
    Stand in for requests Response object;
    """
    def __init__(self, text):
        self.text = text
        self.status_code = 200


class DummyRequestsSession(object):
    """
    Takes the place of a requests session so that we can check that all
    values are sent and received correctly.
    """
    def __init__(self, backend, urlPrefix):
        self._backend = backend
        self._urlPrefix = urlPrefix
        self._getMethodMap = {
            "datasets": self._backend.runGetDataset,
            "referencesets": self._backend.runGetReferenceSet,
            "references": self._backend.runGetReference,
            "variantsets": self._backend.runGetVariantSet,
            "variants": self._backend.runGetVariant,
            "readgroupsets": self._backend.runGetReadGroupSet,
            "readgroups": self._backend.runGetReadGroup,
            "rnaquantifications": self._backend.runGetRnaQuantification,
            "rnaquantificationsets": self._backend.runGetRnaQuantificationSet,
            "expressionlevels": self._backend.runGetExpressionLevel,
        }
        self._searchMethodMap = {
            "datasets": self._backend.runSearchDatasets,
            "referencesets": self._backend.runSearchReferenceSets,
            "references": self._backend.runSearchReferences,
            "variantsets": self._backend.runSearchVariantSets,
            "variants": self._backend.runSearchVariants,
            "readgroupsets": self._backend.runSearchReadGroupSets,
            "reads": self._backend.runSearchReads,
            "rnaquantifications": self._backend.runSearchRnaQuantifications,
            "rnaquantificationsets":
                self._backend.runSearchRnaQuantificationSets,
            "expressionlevels": self._backend.runSearchExpressionLevels,
        }
        self.headers = {}

    def checkSessionParameters(self):
        contentType = "Content-type"
        assert contentType in self.headers
        assert self.headers[contentType] == "application/json"

    def get(self, url, params):
        # TODO add some more checks for params to see if Key is set,
        # and we're not sending any extra stuff.
        self.checkSessionParameters()
        assert url.startswith(self._urlPrefix)
        suffix = url[len(self._urlPrefix):]
        splits = suffix.split("/")
        assert len(splits) == 3
        assert splits[0] == ''
        datatype, id_ = splits[1:]
        assert datatype in self._getMethodMap
        method = self._getMethodMap[datatype]
        result = method(id_)
        return DummyResponse(result)

    def post(self, url, params=None, data=None):
        self.checkSessionParameters()
        assert url.startswith(self._urlPrefix)
        suffix = url[len(self._urlPrefix):]
        searchSuffix = "/search"
        if suffix.endswith(searchSuffix):
            datatype = suffix[1:-len(searchSuffix)]
            assert datatype in self._searchMethodMap
            method = self._searchMethodMap[datatype]
            result = method(data)
        else:
            # ListReferenceBases is an oddball and needs to be treated
            # separately.
            data = json.loads(data)
            args = protocol.ListReferenceBasesRequest()
            args.reference_id = data.get('referenceId', "")
            args.start = int(data.get('start', 0))
            args.end = int(data.get('end', 0))
            args.page_token = data.get('pageToken', "")
            result = self._backend.runListReferenceBases(
                protocol.toJson(args))
        return DummyResponse(result)


class DummyHttpClient(client.HttpClient):
    """
    Client in which we intercept calls to the underlying requests connection.
    """
    def __init__(self, backend):
        self._urlPrefix = "http://example.com"
        super(DummyHttpClient, self).__init__(self._urlPrefix)
        self._session = DummyRequestsSession(backend, self._urlPrefix)
        self._setup_http_session()


class ExhaustiveListingsMixin(object):
    """
    Tests exhaustive listings using the high-level API with a Simulated
    backend.
    """
    @classmethod
    def setUpClass(cls):
        cls.backend = backend.Backend(datarepo.SimulatedDataRepository(
            randomSeed=100, numDatasets=3,
            numVariantSets=3, numCalls=3, variantDensity=0.5,
            numReferenceSets=3, numReferencesPerReferenceSet=3,
            numReadGroupSets=3, numReadGroupsPerReadGroupSet=3,
            numAlignments=3, numRnaQuantSets=3))
        cls.dataRepo = cls.backend.getDataRepository()

    def setUp(self):
        self.client = self.getClient()

    def verifyObjectList(self, gaObjects, datamodelObjects, getMethod):
        """
        Verifies that the specified list of protocol objects corresponds
        to the specified list of datamodel objects.
        """
        for gaObject, datamodelObject in utils.zipLists(
                gaObjects, datamodelObjects):
            self.assertEqual(gaObject, datamodelObject.toProtocolElement())
            otherGaObject = getMethod(gaObject.id)
            self.assertEqual(gaObject, otherGaObject)

    def testAllDatasets(self):
        datasets = list(self.client.search_datasets())
        self.verifyObjectList(
            datasets, self.dataRepo.getDatasets(), self.client.get_dataset)

    def testAllReferenceSets(self):
        referenceSets = list(self.client.search_reference_sets())
        self.verifyObjectList(
            referenceSets, self.dataRepo.getReferenceSets(),
            self.client.get_reference_set)

    def testAllReferences(self):
        for referenceSet in self.client.search_reference_sets():
            references = list(self.client.search_references(referenceSet.id))
            datamodelReferences = self.dataRepo.getReferenceSet(
                referenceSet.id).getReferences()
            self.verifyObjectList(
                references, datamodelReferences, self.client.get_reference)
            for datamodelReference in datamodelReferences:
                bases = self.client.list_reference_bases(
                    datamodelReference.getId())
                otherBases = datamodelReference.getBases(
                    0, datamodelReference.getLength())
                self.assertEqual(bases, otherBases)

    def testAllVariantSets(self):
        for dataset in self.client.search_datasets():
            variantSets = list(self.client.search_variant_sets(dataset.id))
            datamodelVariantSets = self.dataRepo.getDataset(
                dataset.id).getVariantSets()
            self.verifyObjectList(
                variantSets, datamodelVariantSets, self.client.get_variant_set)

    def testAllVariants(self):
        for datamodelDataset in self.dataRepo.getDatasets():
            for datamodelVariantSet in datamodelDataset.getVariantSets():
                # TODO the values should be derived from the datamodel
                # variant set object.
                start = 0
                end = 20
                referenceName = "fixme"
                variants = list(self.client.search_variants(
                    datamodelVariantSet.getId(), start=start, end=end,
                    reference_name=referenceName))
                datamodelVariants = [
                    DatamodelObjectWrapper(variant) for variant in
                    datamodelVariantSet.getVariants(
                        referenceName, start, end)]
                self.verifyObjectList(
                    variants, datamodelVariants, self.client.get_variant)

    def testAllReadGroupSets(self):
        for dataset in self.client.search_datasets():
            readGroupSets = list(
                self.client.search_read_group_sets(dataset.id))
            datamodelReadGroupSets = self.dataRepo.getDataset(
                dataset.id).getReadGroupSets()
            self.verifyObjectList(
                readGroupSets, datamodelReadGroupSets,
                self.client.get_read_group_set)
            # Check the readGroups.
            for readGroupSet, datamodelReadGroupSet in zip(
                    readGroupSets, datamodelReadGroupSets):
                datamodelReadGroups = datamodelReadGroupSet.getReadGroups()
                self.verifyObjectList(
                    readGroupSet.read_groups, datamodelReadGroups,
                    self.client.get_read_group)

    def testAllReads(self):
        for dmDataset in self.dataRepo.getDatasets():
            for dmReadGroupSet in dmDataset.getReadGroupSets():
                dmReferenceSet = dmReadGroupSet.getReferenceSet()
                for dmReadGroup in dmReadGroupSet.getReadGroups():
                    for dmReference in dmReferenceSet.getReferences():
                        # TODO fix these coordinates.
                        start = 0
                        end = 10
                        dmReads = list(dmReadGroup.getReadAlignments(
                            dmReference, start, end))
                        reads = list(self.client.search_reads(
                            [dmReadGroup.getId()], dmReference.getId(),
                            start, end))
                        self.assertGreater(len(reads), 0)
                        for dmRead, read in utils.zipLists(dmReads, reads):
                            self.assertEqual(dmRead, read)

    def testAllRnaQuantificationSets(self):
        for dataset in self.client.search_datasets():
            rnaQuantificationSets = \
                list(self.client.search_rna_quantification_sets(dataset.id))
            datamodelRnaQuantificationSets = self.dataRepo.getDataset(
                dataset.id).getRnaQuantificationSets()
            self.verifyObjectList(
                rnaQuantificationSets, datamodelRnaQuantificationSets,
                self.client.get_rna_quantification_set)


class TestExhaustiveListingsHttp(ExhaustiveListingsMixin, unittest.TestCase):
    """
    Tests the exhaustive listings using the HTTP client.
    """

    def getClient(self):
        return DummyHttpClient(self.backend)


class TestExhaustiveListingsLocal(ExhaustiveListingsMixin, unittest.TestCase):
    """
    Tests the exhaustive listings using the local client.
    """

    def getClient(self):
        return client.LocalClient(self.backend)


class PagingMixin(object):
    """
    Tests the paging code using a simulated backend.
    """
    @classmethod
    def setUpClass(cls):
        cls.numReferences = 25
        cls.backend = backend.Backend(datarepo.SimulatedDataRepository(
            randomSeed=100, numDatasets=0,
            numReferenceSets=1,
            numReferencesPerReferenceSet=cls.numReferences))
        cls.dataRepo = cls.backend.getDataRepository()

    def setUp(self):
        self.client = self.getClient()
        self.datamodelReferenceSet = self.dataRepo.getReferenceSetByIndex(0)
        self.datamodelReferences = self.datamodelReferenceSet.getReferences()
        self.references = [
            dmReference.toProtocolElement()
            for dmReference in self.datamodelReferences]
        self.assertEqual(len(self.references), self.numReferences)

    def verifyAllReferences(self):
        """
        Verifies that we correctly return all references.
        """
        references = list(self.client.search_references(
            self.datamodelReferenceSet.getId()))
        self.assertEqual(references, self.references)

    def testDefaultPageSize(self):
        self.verifyAllReferences()

    def verifyPageSize(self, pageSize):
        self.client.set_page_size(pageSize)
        self.assertEqual(pageSize, self.client.get_page_size())
        self.verifyAllReferences()

    def testPageSize1(self):
        self.verifyPageSize(1)

    def testPageSize2(self):
        self.verifyPageSize(2)

    def testPageSize3(self):
        self.verifyPageSize(3)

    def testPageSizeAlmostListLength(self):
        self.verifyPageSize(self.numReferences - 1)

    def testPageSizeListLength(self):
        self.verifyPageSize(self.numReferences)


class TestPagingLocal(PagingMixin, unittest.TestCase):
    """
    Tests paging using the local client.
    """

    def getClient(self):
        return client.LocalClient(self.backend)


class TestPagingHttp(PagingMixin, unittest.TestCase):
    """
    Tests paging using the HTTP client.
    """

    def getClient(self):
        return DummyHttpClient(self.backend)
