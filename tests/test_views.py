"""
Unit tests for the frontend code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.server as server
import ga4gh.backends as backends
import ga4gh.protocol as protocol


class TestFrontend(unittest.TestCase):

    def setUp(self):
        server.app.config['TESTING'] = True
        server.app.backend = backends.MockBackend()
        self.app = server.app.test_client()

    def testServer(self):
        self.assertEqual(404, self.app.get('/').status_code)

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
        response = self.app.post('/variants/search',
                                 headers={'Content-type': 'application/json'},
                                 data='{"variantSetIds": [1, 2]}')
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantsResponse.fromJSONString(
            response.data)
        self.assertEqual(responseData.variants, [])

    def testVariantSetsSearch(self):
        response = self.app.post('/variantsets/search',
                                 headers={'Content-type': 'application/json'},
                                 data='{"dataSetIds": [1, 2]}')
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantSetsResponse.fromJSONString(
            response.data)
        self.assertEqual(responseData.variantSets, [])
