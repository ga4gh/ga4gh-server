"""
Unit tests for the frontend code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.frontend as frontend
import ga4gh.backend as backend
import ga4gh.protocol as protocol
import tests.utils as utils


class TestFrontend(unittest.TestCase):

    def setUp(self):
        frontend.configure("TestWithoutValidationConfig")
        frontend.app.backend = backend.MockBackend()
        self.app = frontend.app.test_client()

    def sendVariantsSearch(self, data):
        path = utils.applyVersion('/variants/search')
        return self.app.post(path,
                             headers={'Content-type': 'application/json'},
                             data=data)

    def sendVariantSetsSearch(self, data):
        path = utils.applyVersion('/variantsets/search')
        return self.app.post(path,
                             headers={'Content-type': 'application/json'},
                             data=data)

    def sendCallSetsSearch(self, data):
        path = utils.applyVersion('/callsets/search')
        return self.app.post(path,
                             headers={'Content-type': 'application/json'},
                             data=data)

    def test404sReturnJson(self):
        path = utils.applyVersion('/doesNotExist')
        response = self.app.get(path)
        protocol.GAException.fromJsonString(response.get_data())
        self.assertEqual(404, response.status_code)

    def testCors(self):
        def assertHeaders(response):
            self.assertEqual('*',
                             response.headers['Access-Control-Allow-Origin'])
            self.assertTrue('Content-Type' in response.headers)

        assertHeaders(self.app.get('/'))
        assertHeaders(self.sendVariantsSearch('{"variantSetIds": [1, 2]}'))
        assertHeaders(self.sendVariantSetsSearch('{"dataSetIds": [1, 2]}'))
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

        # TODO uncomment when error handling correctly implemented
        # badResponse = self.app.post(
        #     versionedPath, headers={'Content-type': 'application/json'})
        # self.assertEqual(400, badResponse.status_code)

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

    def testVariantsSearch(self):
        response = self.sendVariantsSearch('{"variantSetIds": [1, 2]}')
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantsResponse.fromJsonString(
            response.data)
        self.assertEqual(responseData.variants, [])

    def testVariantSetsSearch(self):
        response = self.sendVariantSetsSearch('{"dataSetIds": [1, 2]}')
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantSetsResponse.fromJsonString(
            response.data)
        self.assertEqual(responseData.variantSets, [])

    def testCallSetsSearch(self):
        response = self.sendCallSetsSearch('{"callSetId": "xx"}')
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertEqual(responseData.callSets, [])

    def testWrongVersion(self):
        path = '/v0.1.2/variantsets/search'
        self.assertEqual(404, self.app.options(path).status_code)
