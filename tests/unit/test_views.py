"""
Unit tests for the frontend code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import logging

import ga4gh.frontend as frontend
import ga4gh.protocol as protocol
import tests.utils as utils


class TestFrontend(unittest.TestCase):
    """
    Tests the basic routing and HTTP handling for the Flask app.
    """
    exampleUrl = 'www.example.com'

    @classmethod
    def setUpClass(cls):
        config = {
            "DATA_SOURCE": "__SIMULATED__",
            "SIMULATED_BACKEND_RANDOM_SEED": 1111,
            "SIMULATED_BACKEND_NUM_CALLS": 1,
            "SIMULATED_BACKEND_VARIANT_DENSITY": 1.0,
            "SIMULATED_BACKEND_NUM_VARIANT_SETS": 1,
            # "DEBUG" : True
        }
        reload(frontend)
        frontend.configure(
            baseConfig="TestConfig", extraConfig=config)
        cls.app = frontend.app.test_client()
        # silence usually unhelpful CORS log
        logging.getLogger('ga4gh.frontend.cors').setLevel(logging.CRITICAL)

    @classmethod
    def tearDownClass(cls):
        cls.app = None

    def sendPostRequest(self, path, request):
        """
        Sends the specified GA request object and returns the response.
        """
        versionedPath = utils.applyVersion(path)
        headers = {
            'Content-type': 'application/json',
            'Origin': self.exampleUrl,
        }
        return self.app.post(
            versionedPath, headers=headers,
            data=request.toJsonString())

    def sendVariantsSearch(self):
        response = self.sendVariantSetsSearch()
        variantSets = protocol.SearchVariantSetsResponse().fromJsonString(
            response.data).variantSets
        request = protocol.SearchVariantsRequest()
        request.variantSetIds = [variantSets[0].id]
        request.referenceName = "1"
        request.start = 0
        request.end = 1
        return self.sendPostRequest('/variants/search', request)

    def sendVariantSetsSearch(self):
        request = protocol.SearchVariantSetsRequest()
        request.datasetIds = ["simulatedDataset1"]
        return self.sendPostRequest('/variantsets/search', request)

    def sendCallSetsSearch(self):
        response = self.sendVariantSetsSearch()
        variantSets = protocol.SearchVariantSetsResponse().fromJsonString(
            response.data).variantSets
        request = protocol.SearchCallSetsRequest()
        request.variantSetIds = [variantSets[0].id]
        return self.sendPostRequest('/callsets/search', request)

    def sendReadsSearch(self, readGroupIds=None):
        if readGroupIds is None:
            readGroupIds = ['aReadGroupSet:one']
        request = protocol.SearchReadsRequest()
        request.readGroupIds = readGroupIds
        return self.sendPostRequest('/reads/search', request)

    def sendDatasetsSearch(self):
        request = protocol.SearchDatasetsRequest()
        return self.sendPostRequest('/datasets/search', request)

    def sendGetRequest(self, path):
        versionedPath = utils.applyVersion(path)
        headers = {
            'Origin': self.exampleUrl,
        }
        response = self.app.get(versionedPath, headers=headers)
        return response

    def sendReferencesGet(self, id_=None):
        if id_ is None:
            id_ = 'simple:simple'
        path = "/references/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendReferenceSetsGet(self, id_=None):
        if id_ is None:
            id_ = 'simple'
        path = "/referencesets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendReferencesSearch(self):
        path = "/references/search"
        request = protocol.SearchReferencesRequest()
        response = self.sendPostRequest(path, request)
        return response

    def sendListRequest(self, path, request):
        versionedPath = utils.applyVersion(path)
        headers = {
            'Origin': self.exampleUrl,
        }
        data = request.toJsonDict()
        response = self.app.get(
            versionedPath, data=data, headers=headers)
        return response

    def sendReferenceBasesList(self, id_=None):
        if id_ is None:
            id_ = 'simple:simple'
        path = "/references/{}/bases".format(id_)
        request = protocol.ListReferenceBasesRequest()
        response = self.sendListRequest(path, request)
        return response

    def test404sReturnJson(self):
        paths = [
            '/doesNotExist',
            utils.applyVersion('/doesNotExist'),
            utils.applyVersion('/reads/sea'),
            utils.applyVersion('/variantsets/id/doesNotExist'),
        ]
        for path in paths:
            response = self.app.get(path)
            protocol.GAException.fromJsonString(response.get_data())
            self.assertEqual(404, response.status_code)

    def testCors(self):
        def assertHeaders(response):
            self.assertEqual(self.exampleUrl,
                             response.headers['Access-Control-Allow-Origin'])
            self.assertTrue('Content-Type' in response.headers)

        assertHeaders(self.sendVariantsSearch())
        assertHeaders(self.sendVariantSetsSearch())
        assertHeaders(self.sendReadsSearch())
        assertHeaders(self.sendReferencesGet())
        assertHeaders(self.sendReferenceSetsGet())
        assertHeaders(self.sendReferencesSearch())
        assertHeaders(self.sendReferenceBasesList())
        assertHeaders(self.sendDatasetsSearch())
        # TODO: Test other methods as they are implemented

    def verifySearchRouting(self, path, getDefined=False):
        """
        Verifies that the specified path has the correct routing for a
        search command. If getDefined is False we check to see if it
        returns the correct status code.
        """
        versionedPath = utils.applyVersion(path)
        response = self.app.post(versionedPath)
        protocol.GAException.fromJsonString(response.get_data())
        self.assertEqual(415, response.status_code)
        if not getDefined:
            getResponse = self.app.get(versionedPath)
            protocol.GAException.fromJsonString(getResponse.get_data())
            self.assertEqual(405, getResponse.status_code)

        # Malformed requests should return 400
        for badJson in ["", None, "JSON", "<xml/>", "{]"]:
            badResponse = self.app.post(
                versionedPath, data=badJson,
                headers={'Content-type': 'application/json'})
            self.assertEqual(400, badResponse.status_code)

        # OPTIONS should return success
        self.assertEqual(200, self.app.options(versionedPath).status_code)

    def testRouteReferences(self):
        referenceId = "referenceSet0:srs0"
        paths = ['/references/{}', '/references/{}/bases']
        for path in paths:
            path = path.format(referenceId)
            versionedPath = utils.applyVersion(path)
            self.assertEqual(200, self.app.get(versionedPath).status_code)
        referenceSetId = "referenceSet0"
        paths = ['/referencesets/{}']
        for path in paths:
            path = path.format(referenceSetId)
            versionedPath = utils.applyVersion(path)
            self.assertEqual(200, self.app.get(versionedPath).status_code)
        self.verifySearchRouting('/referencesets/search', True)
        self.verifySearchRouting('/references/search', True)

    def testRouteCallsets(self):
        path = utils.applyVersion('/callsets/search')
        self.assertEqual(415, self.app.post(path).status_code)
        self.assertEqual(200, self.app.options(path).status_code)
        self.assertEqual(405, self.app.get(path).status_code)

    def testRouteReads(self):
        paths = ['/reads/search', '/readgroupsets/search']
        for path in paths:
            self.verifySearchRouting(path)

    def testRouteVariants(self):
        for path in ['/variantsets/search', '/variants/search']:
            self.verifySearchRouting(path)

    def testRouteIndex(self):
        self._routeIndex("/")

    def testRouteIndexRedirect(self):
        self._routeIndex("/{}".format(protocol.version))

    def _routeIndex(self, path):
        response = self.app.get(path)
        self.assertEqual(200, response.status_code)
        self.assertEqual("text/html", response.mimetype)
        self.assertGreater(len(response.data), 0)

    def testVariantsSearch(self):
        response = self.sendVariantsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchVariantsResponse.fromJsonString(
            response.data)
        self.assertEqual(len(responseData.variants), 1)

    def testVariantSetsSearch(self):
        response = self.sendVariantSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchVariantSetsResponse.fromJsonString(
            response.data)
        self.assertEqual(len(responseData.variantSets), 1)

    def testCallSetsSearch(self):
        response = self.sendCallSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertEqual(len(responseData.callSets), 1)

    def testReadsSearch(self):
        response = self.sendReadsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchReadsResponse.fromJsonString(
            response.data)
        self.assertEqual(len(responseData.alignments), 2)
        self.assertEqual(
            responseData.alignments[0].id,
            "aReadGroupSet:one:simulated0")
        self.assertEqual(
            responseData.alignments[1].id,
            "aReadGroupSet:one:simulated1")

    def testDatasetsSearch(self):
        response = self.sendDatasetsSearch()
        responseData = protocol.SearchDatasetsResponse.fromJsonString(
            response.data)
        datasets = list(responseData.datasets)
        self.assertEqual('simulatedDataset1', datasets[0].id)

    def testWrongVersion(self):
        path = '/v0.1.2/variantsets/search'
        self.assertEqual(404, self.app.options(path).status_code)

    def testCurrentVersion(self):
        path = '/{}/variantsets/search'.format(
            frontend.Version.currentString)
        self.assertEqual(200, self.app.options(path).status_code)

    def testNotImplementedPaths(self):
        pathsNotImplementedPost = [
            '/genotypephenotype/search',
            '/individuals/search',
            '/samples/search',
            '/experiments/search',
            '/individualgroups/search',
            '/analyses/search',
            '/sequences/search',
            '/joins/search',
            '/subgraph/segments',
            '/subgraph/joins',
            '/features/search',
            '/variantsets/<id>/sequences/search',
            '/alleles/search',
        ]
        pathsNotImplementedGet = [
            '/callsets/<id>',
            '/alleles/<id>',
            '/variants/<id>',
            '/variantsets/<id>/sequences/<id>',
            '/variantsets/<id>',
            '/feature/<id>',
            '/sequences/<id>/bases',
            '/mode/<id>',
            '/datasets/<id>',
            '/readgroupsets/<id>',
            '/readgroups/<id>',
        ]

        def runRequest(method, path):
            requestPath = path.replace('<id>', 'someId')
            versionedPath = utils.applyVersion(requestPath)
            response = method(versionedPath)
            protocol.GAException.fromJsonString(response.get_data())
            self.assertEqual(response.status_code, 501)

        for path in pathsNotImplementedGet:
            runRequest(self.app.get, path)
        for path in pathsNotImplementedPost:
            runRequest(self.app.post, path)

    def testNoAuthentication(self):
        path = '/oauth2callback'.format(
            frontend.Version.currentString)
        self.assertEqual(501, self.app.get(path).status_code)
