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


class TestFrontend(unittest.TestCase):

    def setUp(self):
        frontend.configure("TestConfig")
        frontend.app.backend = backend.MockBackend()
        self.app = frontend.app.test_client()

    def sendVariantsSearch(self, data):
        return self.app.post('/variants/search',
                             headers={'Content-type': 'application/json'},
                             data=data)

    def sendVariantSetsSearch(self, data):
        return self.app.post('/variantsets/search',
                             headers={'Content-type': 'application/json'},
                             data=data)

    def testServer(self):
        self.assertEqual(404, self.app.get('/doesNotExist').status_code)

    def testCors(self):
        def assertHeaders(response):
            self.assertEqual('*',
                             response.headers['Access-Control-Allow-Origin'])
            self.assertTrue('Content-Type' in response.headers)

        assertHeaders(self.app.get('/'))
        assertHeaders(self.sendVariantsSearch('{"variantSetIds": [1, 2]}'))
        assertHeaders(self.sendVariantSetsSearch('{"dataSetIds": [1, 2]}'))
        # TODO: Test other methods as they are implemented

    def testRouteReferences(self):
        for path in ['/references/1', 'references/1/bases', 'referencesets/1']:
            self.assertEqual(404, self.app.get(path).status_code)

        for path in ['/referencesets/search', '/references/search']:
            self.assertEqual(404, self.app.get(path).status_code)

    def testRouteCallsets(self):
        self.assertEqual(404, self.app.post('/callsets/search').status_code)

    def testRouteReads(self):
        for path in ['/reads/search', '/readgroupsets/search']:
            self.assertEqual(404, self.app.post(path).status_code)

    def testRouteVariants(self):
        for path in ['/variantsets/search', '/variants/search']:
            # POST gets routed to input checking
            self.assertEqual(415, self.app.post(path).status_code)

            # TODO disabling this test until we correctly implement error
            # handling in the frontend.
            # self.assertEqual(400, self.app.post(
            #     path,
            #     headers={'Content-type': 'application/json'}).status_code)

            # OPTIONS should return success
            self.assertEqual(200, self.app.options(path).status_code)

            # GET is not defined
            self.assertEqual(405, self.app.get(path).status_code)

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
