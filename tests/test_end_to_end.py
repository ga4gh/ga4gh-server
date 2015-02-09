"""
End-to-end tests for each backend configuration. Sets up a server with
the backend, sends some basic queries to that server and verifies results
are as expected.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

from tests.test_backends import WormtableTestFixture
import ga4gh.backend as backend
import ga4gh.datamodel.variants as variants
import ga4gh.frontend as frontend
import ga4gh.protocol as protocol


class EndToEndWormtableTest(unittest.TestCase):
    def setUpServer(self, backend):
        frontend.app.config['TESTING'] = True
        frontend.app.backend = backend
        self.app = frontend.app.test_client()

    def setUp(self):
        self._wtTestFixture = WormtableTestFixture()
        self._wtTestFixture.setUp()
        self.setUpServer(backend.Backend(self._wtTestFixture.dataDir,
                                         variants.WormtableVariantSet))

    def tearDown(self):
        self._wtTestFixture.tearDown()

    def sendJSONPostRequest(self, path, data):
        return self.app.post(path,
                             headers={'Content-type': 'application/json'},
                             data=data)

    def testVariantSetsSearch(self):
        # TODO: Set up the test backend API to return these values rather than
        # hard-coding them here.
        expectedIds = ['example_1', 'example_2', 'example_3', 'example_4']
        request = protocol.GASearchVariantSetsRequest()
        request.pageSize = len(expectedIds)
        response = self.sendJSONPostRequest('/variantsets/search',
                                            request.toJsonString())

        self.assertEqual(200, response.status_code)

        responseData = protocol.GASearchVariantSetsResponse.fromJsonString(
            response.data)
        self.assertTrue(protocol.GASearchVariantSetsResponse.validate(
            responseData.toJsonDict()))

        self.assertIsNone(responseData.nextPageToken)
        self.assertEqual(len(expectedIds), len(responseData.variantSets))
        for variantSet in responseData.variantSets:
            self.assertTrue(variantSet.id in expectedIds)

    def testVariantsSearch(self):
        # TODO: As above, get these from the test backend API
        expectedIds = ['example_1']
        referenceName = '1'

        request = protocol.GASearchVariantsRequest()
        request.referenceName = referenceName
        request.start = 0
        request.end = 0

        request.variantSetIds = expectedIds

        # Request windows is too small, no results
        response = self.sendJSONPostRequest('/variants/search',
                                            request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantsResponse.fromJsonString(
            response.data)
        self.assertIsNone(responseData.nextPageToken)
        self.assertEqual([], responseData.variants)

        # Larger request window, expect results
        request.end = 2 ** 16
        response = self.sendJSONPostRequest('/variants/search',
                                            request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantsResponse.fromJsonString(
            response.data)
        self.assertTrue(protocol.GASearchVariantsResponse.validate(
            responseData.toJsonDict()))
        self.assertNotEqual([], responseData.variants)

        # Verify all results are in the correct range, set and reference
        for variant in responseData.variants:
            self.assertGreaterEqual(variant.start, 0)
            self.assertLessEqual(variant.end, 2 ** 16)
            self.assertTrue(variant.variantSetId in expectedIds)
            self.assertEqual(variant.referenceName, referenceName)

        # TODO: Add more useful test scenarios, including some covering
        # pagination behavior.

        # TODO: Add test cases for other methods when they are implemented.

    # TODO: Add tests for other backends
