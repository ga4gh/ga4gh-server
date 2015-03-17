"""
End-to-end tests for the simulator configuration. Sets up a server with
the backend, sends some basic queries to that server and verifies results
are as expected.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
from nose.tools import nottest  # TODO remove when we fix Calls

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
        "SIMULATED_BACKEND_NUM_VARIANT_SETS": 10,
    }
    frontend.configure(
        baseConfig="TestConfig", extraConfig=config)
    global _app
    _app = frontend.app.test_client()


def tearDown(self):
    global _app
    _app = None


class TestSimulatedStack(unittest.TestCase):
    """
    Tests the full stack for the Simulated backend by using the Flask
    testing client.
    """
    def setUp(self):
        global _app
        self.app = _app
        self.backend = frontend.app.backend
        self.variantSetIds = [
            variantSet.getId() for variantSet in
            self.backend.getVariantSets()]

    def sendJsonPostRequest(self, path, data):
        return self.app.post(
            path, headers={'Content-type': 'application/json'},
            data=data)

    def testVariantSetsSearch(self):
        expectedIds = self.variantSetIds
        request = protocol.GASearchVariantSetsRequest()
        request.pageSize = len(expectedIds)
        path = utils.applyVersion('/variantsets/search')
        response = self.sendJsonPostRequest(
            path, request.toJsonString())

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
        expectedIds = self.variantSetIds[:1]
        referenceName = '1'

        request = protocol.GASearchVariantsRequest()
        request.referenceName = referenceName
        request.start = 0
        request.end = 0
        request.variantSetIds = expectedIds

        # Request windows is too small, no results
        path = utils.applyVersion('/variants/search')
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantsResponse.fromJsonString(
            response.data)
        self.assertIsNone(responseData.nextPageToken)
        self.assertEqual([], responseData.variants)

        # Larger request window, expect results
        request.end = 2 ** 16
        path = utils.applyVersion('/variants/search')
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchVariantsResponse.fromJsonString(
            response.data)
        self.assertTrue(protocol.GASearchVariantsResponse.validate(
            responseData.toJsonDict()))
        self.assertGreater(len(responseData.variants), 0)

        # Verify all results are in the correct range, set and reference
        for variant in responseData.variants:
            self.assertGreaterEqual(variant.start, 0)
            self.assertLessEqual(variant.end, 2 ** 16)
            self.assertTrue(variant.variantSetId in expectedIds)
            self.assertEqual(variant.referenceName, referenceName)

        # TODO: Add more useful test scenarios, including some covering
        # pagination behavior.

        # TODO: Add test cases for other methods when they are implemented.

    @nottest
    def testCallSetsSearch(self):
        # TODO remove the @nottest decorator here once calls have been
        # properly implemented in the simulator.
        request = protocol.GASearchCallSetsRequest()
        request.name = None
        path = utils.applyVersion('/callsets/search')

        # when variantSetIds are wrong, no results
        request.variantSetIds = ["xxxx"]
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertIsNone(responseData.nextPageToken)
        self.assertEqual([], responseData.callSets)

        # if no callset name is given return all callsets
        request.variantSetIds = self.variantSetIds[:1]
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.GASearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertTrue(protocol.GASearchCallSetsResponse.validate(
            responseData.toJsonDict()))
        self.assertNotEqual([], responseData.callSets)
        # TODO test the length of responseData.callSets equal to all callsets

        # Verify all results are of the correct type and range
        for callSet in responseData.callSets:
            self.assertIs(type(callSet.info), dict)
            self.assertIs(type(callSet.variantSetIds), list)
            splits = callSet.id.split(".")
            variantSetId = '.'.join(splits[:2])
            callSetName = splits[-1]
            self.assertIn(variantSetId, callSet.variantSetIds)
            self.assertEqual(callSetName, callSet.name)
            self.assertEqual(callSetName, callSet.sampleId)

        # TODO add tests after name string search schemas is implemented
