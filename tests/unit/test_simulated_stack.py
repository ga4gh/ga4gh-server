"""
End-to-end tests for the simulator configuration. Sets up a server with
the backend, sends some basic queries to that server and verifies results
are as expected.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import hashlib
import logging

import ga4gh.frontend as frontend
import ga4gh.protocol as protocol
import tests.utils as utils


class TestSimulatedStack(unittest.TestCase):
    """
    Tests the full stack for the Simulated backend by using the Flask
    testing client.
    """
    @classmethod
    def setUpClass(cls):
        # silence usually unhelpful CORS log
        logging.getLogger('ga4gh.frontend.cors').setLevel(logging.CRITICAL)

        cls.numReferenceSets = 5
        cls.numReferencesPerReferenceSet = 3
        config = {
            "DATA_SOURCE": "__SIMULATED__",
            "SIMULATED_BACKEND_RANDOM_SEED": 1111,
            "SIMULATED_BACKEND_NUM_CALLS": 0,
            "SIMULATED_BACKEND_VARIANT_DENSITY": 1.0,
            "SIMULATED_BACKEND_NUM_VARIANT_SETS": 10,
            "SIMULATED_BACKEND_NUM_REFERENCE_SETS": cls.numReferenceSets,
            "SIMULATED_BACKEND_NUM_REFERENCES_PER_REFERENCE_SET":
                cls.numReferencesPerReferenceSet,
        }
        frontend.configure(
            baseConfig="TestConfig", extraConfig=config)
        cls.app = frontend.app.test_client()

    @classmethod
    def tearDownClass(cls):
        cls.app = None

    def setUp(self):
        self.backend = frontend.app.backend
        self.dataset_id = self.backend.getDatasetIds()[0]
        self.variant_set_ids = [
            variantSet.getId() for variantSet in
            self.backend.getDataset(self.dataset_id).getVariantSets()]

    def sendJsonPostRequest(self, path, data):
        return self.app.post(
            path, headers={'Content-type': 'application/json'},
            data=data)

    def testVariantSetsSearch(self):
        expectedIds = self.variant_set_ids
        request = protocol.SearchVariantSetsRequest()
        request.page_size = len(expectedIds)
        request.dataset_ids = [self.dataset_id]
        path = utils.applyVersion('/variantsets/search')
        response = self.sendJsonPostRequest(
            path, request.toJsonString())

        self.assertEqual(200, response.status_code)

        responseData = protocol.SearchVariantSetsResponse.fromJsonString(
            response.data)
        self.assertTrue(protocol.SearchVariantSetsResponse.validate(
            responseData.toJsonDict()))

        self.assertIsNone(responseData.next_page_token)
        self.assertEqual(len(expectedIds), len(responseData.variant_sets))
        for variantSet in responseData.variant_sets:
            self.assertTrue(variantSet.id in expectedIds)

    def testVariantsSearch(self):
        expectedIds = self.variant_set_ids[:1]
        reference_name = '1'

        request = protocol.SearchVariantsRequest()
        request.reference_name = reference_name
        request.start = 0
        request.end = 0
        request.variant_set_ids = expectedIds

        # Request windows is too small, no results
        path = utils.applyVersion('/variants/search')
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchVariantsResponse.fromJsonString(
            response.data)
        self.assertIsNone(responseData.next_page_token)
        self.assertEqual([], responseData.variants)

        # Larger request window, expect results
        request.end = 2 ** 16
        path = utils.applyVersion('/variants/search')
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchVariantsResponse.fromJsonString(
            response.data)
        self.assertTrue(protocol.SearchVariantsResponse.validate(
            responseData.toJsonDict()))
        self.assertGreater(len(responseData.variants), 0)

        # Verify all results are in the correct range, set and reference
        for variant in responseData.variants:
            self.assertGreaterEqual(variant.start, 0)
            self.assertLessEqual(variant.end, 2 ** 16)
            self.assertTrue(variant.variant_set_id in expectedIds)
            self.assertEqual(variant.reference_name, reference_name)

        # TODO: Add more useful test scenarios, including some covering
        # pagination behavior.

        # TODO: Add test cases for other methods when they are implemented.

    @unittest.skipIf(True, "")
    def testCallSetsSearch(self):
        # TODO remove the @skipIf decorator here once calls have been
        # properly implemented in the simulator.
        request = protocol.SearchCallSetsRequest()
        request.name = None
        path = utils.applyVersion('/callsets/search')

        # when variant_set_ids are wrong, no results
        request.variant_set_ids = ["xxxx"]
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertIsNone(responseData.next_page_token)
        self.assertEqual([], responseData.call_sets)

        # if no callset name is given return all callsets
        request.variant_set_ids = self.variant_set_ids[:1]
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertTrue(protocol.SearchCallSetsResponse.validate(
            responseData.toJsonDict()))
        self.assertNotEqual([], responseData.call_sets)
        # TODO test the length of responseData.call_sets equal to all callsets

        # Verify all results are of the correct type and range
        for callSet in responseData.call_sets:
            self.assertIs(type(callSet.info), dict)
            self.assertIs(type(callSet.variant_set_ids), list)
            splits = callSet.id.split(".")
            variant_set_id = '.'.join(splits[:2])
            callSetName = splits[-1]
            self.assertIn(variant_set_id, callSet.variant_set_ids)
            self.assertEqual(callSetName, callSet.name)
            self.assertEqual(callSetName, callSet.sampleId)

        # TODO add tests after name string search schemas is implemented

    def testReferences(self):
        # search for reference sets
        path = utils.applyVersion('/referencesets/search')
        request = protocol.SearchReferenceSetsRequest()
        response = self.sendJsonPostRequest(path, request.toJsonString())
        self.assertEqual(response.status_code, 200)
        responseData = protocol.SearchReferenceSetsResponse.fromJsonString(
            response.data)
        reference_sets = responseData.reference_sets
        self.assertEqual(self.numReferenceSets, len(reference_sets))

        # search for references
        path = utils.applyVersion('/references/search')
        request = protocol.SearchReferencesRequest()
        response = self.sendJsonPostRequest(path, request.toJsonString())
        self.assertEqual(response.status_code, 200)
        responseData = protocol.SearchReferencesResponse.fromJsonString(
            response.data)
        references = responseData.references
        self.assertEqual(
            self.numReferenceSets * self.numReferencesPerReferenceSet,
            len(references))

        for referenceSet in reference_sets:
            # fetch the reference set
            path = utils.applyVersion(
                '/referencesets/{}'.format(referenceSet.id))
            response = self.app.get(path)
            self.assertEqual(response.status_code, 200)
            fetchedReferenceSet = protocol.ReferenceSet.fromJsonString(
                response.data)
            self.assertEqual(fetchedReferenceSet, referenceSet)
            self.assertEqual(
                len(fetchedReferenceSet.reference_ids),
                self.numReferencesPerReferenceSet)

            for reference_id in referenceSet.reference_ids:
                # fetch the reference
                path = utils.applyVersion(
                    '/references/{}'.format(reference_id))
                response = self.app.get(path)
                self.assertEqual(response.status_code, 200)
                fetchedReference = protocol.Reference.fromJsonString(
                    response.data)
                self.assertEqual(fetchedReference.id, reference_id)

                # fetch the bases
                path = utils.applyVersion(
                    '/references/{}/bases'.format(reference_id))
                args = protocol.ListReferenceBasesRequest().toJsonDict()
                response = self.app.get(path, data=args)
                self.assertEqual(response.status_code, 200)
                bases = protocol.ListReferenceBasesResponse.fromJsonString(
                    response.data)
                self.assertEqual(len(bases.sequence), 200)
                self.assertEqual(
                    set(bases.sequence),
                    set(['A', 'C', 'T', 'G']))
                calculatedDigest = hashlib.md5(bases.sequence).hexdigest()
                self.assertEqual(
                    calculatedDigest, fetchedReference.md5checksum)
