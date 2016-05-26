"""
Unit tests for the frontend code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import logging

import tests.paths as paths

import ga4gh.datamodel as datamodel
import ga4gh.frontend as frontend
import ga4gh.protocol as protocol


class TestFrontend(unittest.TestCase):
    """
    Tests the basic routing and HTTP handling for the Flask app.
    """
    exampleUrl = 'www.example.com'

    @classmethod
    def setUpClass(cls):
        config = {
            "DATA_SOURCE": "simulated://",
            "SIMULATED_BACKEND_RANDOM_SEED": 1111,
            "SIMULATED_BACKEND_NUM_CALLS": 1,
            "SIMULATED_BACKEND_VARIANT_DENSITY": 1.0,
            "SIMULATED_BACKEND_NUM_VARIANT_SETS": 1,
            "LANDING_MESSAGE_HTML": paths.landingMessageHtml
            # "DEBUG" : True
        }
        frontend.reset()
        frontend.configure(
            baseConfig="TestConfig", extraConfig=config)
        cls.app = frontend.app.test_client()
        # silence usually unhelpful CORS log
        logging.getLogger('ga4gh.frontend.cors').setLevel(logging.CRITICAL)

        # example test values
        cls.backend = frontend.app.backend
        cls.dataRepo = cls.backend.getDataRepository()
        cls.referenceSet = cls.dataRepo.getReferenceSets()[0]
        cls.referenceSetId = cls.referenceSet.getId()
        cls.reference = cls.referenceSet.getReferences()[0]
        cls.referenceId = cls.reference.getId()
        cls.dataset = cls.backend.getDataRepository().getDatasets()[0]
        cls.datasetId = cls.dataset.getId()
        cls.variantSet = cls.dataset.getVariantSets()[0]
        cls.variantSetId = cls.variantSet.getId()
        gaVariant = cls.variantSet.getVariants("1", 0, 2**32).next()
        cls.variantId = gaVariant.id
        cls.callSet = cls.variantSet.getCallSets()[0]
        cls.callSetId = cls.callSet.getId()
        cls.readGroupSet = cls.dataset.getReadGroupSets()[0]
        cls.readGroupSetId = cls.readGroupSet.getId()
        cls.readGroup = cls.readGroupSet.getReadGroups()[0]
        cls.readGroupId = cls.readGroup.getId()
        cls.readAlignment = cls.readGroup.getReadAlignments().next()
        cls.readAlignmentId = cls.readAlignment.id

    def sendPostRequest(self, path, request):
        """
        Sends the specified GA request object and returns the response.
        """
        headers = {
            'Content-type': 'application/json',
            'Origin': self.exampleUrl,
        }
        return self.app.post(
            path, headers=headers, data=protocol.toJson(request))

    def sendGetRequest(self, path):
        """
        Sends a get request to the specified URL and returns the response.
        """
        headers = {
            'Origin': self.exampleUrl,
        }
        return self.app.get(path, headers=headers)

    def sendVariantsSearch(self):
        response = self.sendVariantSetsSearch()
        variantSets = protocol.fromJson(
            response.data, protocol.SearchVariantSetsResponse).variant_sets
        request = protocol.SearchVariantsRequest()
        request.variant_set_id = variantSets[0].id
        request.reference_name = "1"
        request.start = 0
        request.end = 1
        return self.sendPostRequest('/variants/search', request)

    def sendVariantSetsSearch(self):
        request = protocol.SearchVariantSetsRequest()
        request.dataset_id = self.datasetId
        return self.sendPostRequest('/variantsets/search', request)

    def sendCallSetsSearch(self):
        response = self.sendVariantSetsSearch()
        variantSets = protocol.fromJson(
            response.data, protocol.SearchVariantSetsResponse).variant_sets
        request = protocol.SearchCallSetsRequest()
        request.variant_set_id = variantSets[0].id
        return self.sendPostRequest('/callsets/search', request)

    def sendReadsSearch(self, readGroupIds=None, referenceId=""):
        request = protocol.SearchReadsRequest()
        request.read_group_ids.extend(readGroupIds)
        request.reference_id = referenceId
        return self.sendPostRequest('/reads/search', request)

    def sendDatasetsSearch(self):
        request = protocol.SearchDatasetsRequest()
        return self.sendPostRequest('/datasets/search', request)

    def sendReferencesSearch(self):
        path = "/references/search"
        request = protocol.SearchReferencesRequest()
        response = self.sendPostRequest(path, request)
        return response

    def sendGetVariant(self, id_=None):
        if id_ is None:
            id_ = self.variantId
        path = "/variants/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetVariantSet(self, id_=None):
        if id_ is None:
            id_ = self.variantSetId
        path = "/variantsets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetDataset(self, id_=None):
        if id_ is None:
            id_ = self.datasetId
        path = "/datasets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetReadGroup(self, id_=None):
        if id_ is None:
            id_ = self.readGroupId
        path = "/readgroups/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetReference(self, id_=None):
        if id_ is None:
            id_ = self.referenceId
        path = "/references/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetReadGroupSet(self, id_=None):
        if id_ is None:
            id_ = self.readGroupSetId
        path = "/readgroupsets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetCallSet(self, id_=None):
        if id_ is None:
            id_ = self.callSetId
        path = "/callsets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetReferenceSet(self, id_=None):
        if id_ is None:
            id_ = self.referenceSetId
        path = "/referencesets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendListRequest(self, path, request):
        headers = {
            'Origin': self.exampleUrl,
        }
        data = protocol.toJsonDict(request)
        response = self.app.get(path, data=data, headers=headers)
        return response

    def sendReferenceBasesList(self, id_=None):
        if id_ is None:
            id_ = self.referenceId
        path = "/references/{}/bases".format(id_)
        request = protocol.ListReferenceBasesRequest()
        response = self.sendListRequest(path, request)
        return response

    def test404sReturnJson(self):
        paths = [
            '/doesNotExist',
            '/reads/sea',
            '/variantsets/id/doesNotExist',
        ]
        for path in paths:
            response = self.app.get(path)
            protocol.fromJson(
                response.get_data(), protocol.GAException)
            self.assertEqual(404, response.status_code)

    def testCors(self):
        def assertHeaders(response):
            self.assertEqual(self.exampleUrl,
                             response.headers['Access-Control-Allow-Origin'])
            self.assertTrue('Content-Type' in response.headers)
        # Post-based search methods
        assertHeaders(self.sendVariantsSearch())
        assertHeaders(self.sendVariantSetsSearch())
        assertHeaders(self.sendReadsSearch())
        assertHeaders(self.sendReferencesSearch())
        assertHeaders(self.sendReferenceBasesList())
        assertHeaders(self.sendDatasetsSearch())
        # Get-based accessor methods
        assertHeaders(self.sendGetVariantSet())
        assertHeaders(self.sendGetReference())
        assertHeaders(self.sendGetReferenceSet())
        assertHeaders(self.sendGetReadGroupSet())
        assertHeaders(self.sendGetReadGroup())
        assertHeaders(self.sendGetVariant())
        assertHeaders(self.sendGetDataset())
        # TODO: Test other methods as they are implemented

    def verifySearchRouting(self, path, getDefined=False):
        """
        Verifies that the specified path has the correct routing for a
        search command. If getDefined is False we check to see if it
        returns the correct status code.
        """
        response = self.app.post(path)
        protocol.fromJson(
            response.get_data(), protocol.GAException)
        self.assertEqual(415, response.status_code)
        if not getDefined:
            getResponse = self.app.get(path)
            protocol.fromJson(
                getResponse.get_data(), protocol.GAException)
            self.assertEqual(405, getResponse.status_code)

        # Malformed requests should return 400
        for badJson in ["", None, "JSON", "<xml/>", "{]"]:
            badResponse = self.app.post(
                path, data=badJson,
                headers={'Content-type': 'application/json'})
            self.assertEqual(400, badResponse.status_code)

        # OPTIONS should return success
        self.assertEqual(200, self.app.options(path).status_code)

    def testRouteReferences(self):
        referenceId = self.referenceId
        paths = ['/references/{}', '/references/{}/bases']
        for path in paths:
            path = path.format(referenceId)
            self.assertEqual(200, self.app.get(path).status_code)
        referenceSetId = self.referenceSetId
        paths = ['/referencesets/{}']
        for path in paths:
            path = path.format(referenceSetId)
            self.assertEqual(200, self.app.get(path).status_code)
        self.verifySearchRouting('/referencesets/search', True)
        self.verifySearchRouting('/references/search', True)

    def testRouteCallSets(self):
        path = '/callsets/search'
        self.assertEqual(415, self.app.post(path).status_code)
        self.assertEqual(200, self.app.options(path).status_code)
        self.assertEqual(405, self.app.get(path).status_code)

    def testRouteReads(self):
        paths = ['/reads/search', '/readgroupsets/search']
        for path in paths:
            self.verifySearchRouting(path)

    def testRouteVariants(self):
        self.verifySearchRouting('/variantsets/search', True)
        self.verifySearchRouting('/variants/search', False)

    def testRouteIndex(self):
        path = "/"
        response = self.app.get(path)
        self.assertEqual(200, response.status_code)
        self.assertEqual("text/html", response.mimetype)
        self.assertGreater(len(response.data), 0)

    def testVariantsSearch(self):
        response = self.sendVariantsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(
            response.data, protocol.SearchVariantsResponse)
        self.assertEqual(len(responseData.variants), 1)

    def testVariantSetsSearch(self):
        response = self.sendVariantSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(
            response.data, protocol.SearchVariantSetsResponse)
        self.assertEqual(len(responseData.variant_sets), 1)

    def testGetDataset(self):
        # Test OK: ID found
        response = self.sendDatasetsSearch()
        responseData = protocol.fromJson(
            response.data, protocol.SearchDatasetsResponse)
        datasetId = responseData.datasets[0].id
        response = self.sendGetDataset(datasetId)
        self.assertEqual(200, response.status_code)

        # Test Error: 404, ID not found
        invalidId = datamodel.DatasetCompoundId.getInvalidIdString()
        obfuscated = datamodel.CompoundId.obfuscate(invalidId)
        compoundId = datamodel.DatasetCompoundId.parse(obfuscated)
        response = self.sendGetDataset(str(compoundId))
        self.assertEqual(404, response.status_code)

    def testGetVariantSet(self):
        response = self.sendVariantSetsSearch()
        responseData = protocol.fromJson(
            response.data, protocol.SearchVariantSetsResponse)
        variantSetId = responseData.variant_sets[0].id
        response = self.sendGetVariantSet(variantSetId)
        self.assertEqual(200, response.status_code)
        invalidId = datamodel.VariantSetCompoundId.getInvalidIdString()
        obfuscated = datamodel.CompoundId.obfuscate(invalidId)
        compoundId = datamodel.VariantSetCompoundId.parse(obfuscated)
        response = self.sendGetVariantSet(str(compoundId))
        self.assertEqual(404, response.status_code)

    def testGetReadGroupSet(self):
        response = self.sendGetReadGroupSet()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(response.data, protocol.ReadGroupSet)
        self.assertEqual(
            responseData.id, self.readGroupSetId)

    def testGetReadGroup(self):
        response = self.sendGetReadGroup()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(response.data, protocol.ReadGroup)
        self.assertEqual(
            responseData.id, self.readGroupId)

    def testGetCallSet(self):
        response = self.sendGetCallSet()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(response.data, protocol.CallSet)
        self.assertEqual(
            responseData.id, self.callSetId)

    def testGetVariant(self):
        response = self.sendGetVariant()
        self.assertEqual(200, response.status_code)

    def testCallSetsSearch(self):
        response = self.sendCallSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(
            response.data, protocol.SearchCallSetsResponse)
        self.assertEqual(len(responseData.call_sets), 1)

    def testReadsSearch(self):
        response = self.sendReadsSearch(readGroupIds=[self.readGroupId],
                                        referenceId=self.referenceId)
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(
            response.data, protocol.SearchReadsResponse)
        self.assertEqual(len(responseData.alignments), 2)
        self.assertEqual(
            responseData.alignments[0].id,
            self.readAlignmentId)

    def testDatasetsSearch(self):
        response = self.sendDatasetsSearch()
        responseData = protocol.fromJson(
            response.data, protocol.SearchDatasetsResponse)
        datasets = list(responseData.datasets)
        self.assertEqual(self.datasetId, datasets[0].id)

    def testNoAuthentication(self):
        path = '/oauth2callback'
        self.assertEqual(501, self.app.get(path).status_code)

    def testSearchUnmappedReads(self):
        response = self.sendReadsSearch(readGroupIds=[self.readGroupId],
                                        referenceId="")
        self.assertEqual(501, response.status_code)

    def testSearchReadsMultipleReadGroupSetsSetMismatch(self):
        response = self.sendReadsSearch(
            readGroupIds=[self.readGroupId, "42"],
            referenceId=self.referenceId)
        self.assertEqual(400, response.status_code)
