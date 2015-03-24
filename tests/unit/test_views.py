"""
Unit tests for the frontend code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.frontend as frontend
import ga4gh.protocol as protocol
import tests.utils as utils

_app = None


def setUp(self):
    config = {
        "DATA_SOURCE": "__SIMULATED__",
        "SIMULATED_BACKEND_RANDOM_SEED": 1111,
        "SIMULATED_BACKEND_NUM_CALLS": 0,
        "SIMULATED_BACKEND_VARIANT_DENSITY": 1.0,
        "SIMULATED_BACKEND_NUM_VARIANT_SETS": 1,
        # "DEBUG" : True
    }
    frontend.configure(
        baseConfig="TestConfig", extraConfig=config)
    global _app
    _app = frontend.app.test_client()


def tearDown(self):
    global _app
    _app = None


class TestFrontend(unittest.TestCase):
    """
    Tests the basic routing and HTTP handling for the Flask app.
    """
    exampleUrl = 'www.example.com'

    def setUp(self):
        global _app
        self.app = _app

    def sendRequest(self, path, request):
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

    def sendVariantsSearch(
            self, variantSetIds=[""], referenceName="", start=0, end=0):
        request = protocol.GASearchVariantsRequest()
        request.variantSetIds = variantSetIds
        request.referenceName = referenceName
        request.start = start
        request.end = end
        return self.sendRequest('/variants/search', request)

    def sendVariantSetsSearch(self, datasetIds=[""]):
        request = protocol.GASearchVariantSetsRequest()
        request.datasetIds = datasetIds
        return self.sendRequest('/variantsets/search', request)

    def sendCallSetsSearch(self, variantSetIds=[""]):
        request = protocol.GASearchCallSetsRequest()
        request.variantSetIds = variantSetIds
        return self.sendRequest('/callsets/search', request)

    def test404sReturnJson(self):
        path = utils.applyVersion('/doesNotExist')
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
        # TODO: Test other methods as they are implemented

    def verifySearchRouting(self, path, getDefined=False):
        """
        Verifies that the specified path has the correct routing for a search
        command. If getDefined is False we check to see if it returns the
        correct status code.
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
        paths = ['/references/1', 'references/1/bases', 'referencesets/1']
        for path in paths:
            versionedPath = utils.applyVersion(path)
            self.assertEqual(404, self.app.get(versionedPath).status_code)
        paths = ['/references/search']
        for path in paths:
            versionedPath = utils.applyVersion(path)
            self.assertEqual(404, self.app.get(versionedPath).status_code)
        self.verifySearchRouting('/referencesets/search', True)

    def testRouteCallsets(self):
        path = utils.applyVersion('/callsets/search')
        self.assertEqual(415, self.app.post(path).status_code)
        self.assertEqual(200, self.app.options(path).status_code)
        self.assertEqual(405, self.app.get(path).status_code)

    def testRouteReads(self):
        paths = ['/reads/search']
        for path in paths:
            versionedPath = utils.applyVersion(path)
            self.assertEqual(404, self.app.post(versionedPath).status_code)
        for path in ['/readgroupsets/search']:
            self.verifySearchRouting(path)

    def testRouteVariants(self):
        for path in ['/variantsets/search', '/variants/search']:
            self.verifySearchRouting(path)

    def testRouteIndex(self):
        response = self.app.get("/")
        self.assertEqual(200, response.status_code)
        self.assertEqual("text/html", response.mimetype)
        self.assertGreater(len(response.data), 0)

    def testVariantsSearch(self):
        response = self.sendVariantsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantsResponse.fromJsonString(
            response.data)
        self.assertEqual(responseData.variants, [])

    def testVariantSetsSearch(self):
        response = self.sendVariantSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantSetsResponse.fromJsonString(
            response.data)
        self.assertEqual(len(responseData.variantSets), 1)

    def testCallSetsSearch(self):
        response = self.sendCallSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertEqual(responseData.callSets, [])

    def testWrongVersion(self):
        path = '/v0.1.2/variantsets/search'
        self.assertEqual(404, self.app.options(path).status_code)
