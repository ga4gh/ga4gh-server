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
        cls.numAlignmentsPerReadGroup = 2
        config = {
            "DATA_SOURCE": "__SIMULATED__",
            "SIMULATED_BACKEND_RANDOM_SEED": 1111,
            "SIMULATED_BACKEND_NUM_CALLS": 0,
            "SIMULATED_BACKEND_VARIANT_DENSITY": 1.0,
            "SIMULATED_BACKEND_NUM_VARIANT_SETS": 10,
            "SIMULATED_BACKEND_NUM_REFERENCE_SETS": cls.numReferenceSets,
            "SIMULATED_BACKEND_NUM_REFERENCES_PER_REFERENCE_SET":
                cls.numReferencesPerReferenceSet,
            "SIMULATED_BACKEND_NUM_ALIGNMENTS_PER_READ_GROUP":
                cls.numAlignmentsPerReadGroup,
        }
        reload(frontend)
        frontend.configure(
            baseConfig="TestConfig", extraConfig=config)
        cls.app = frontend.app.test_client()

    @classmethod
    def tearDownClass(cls):
        cls.app = None

    def setUp(self):
        self.backend = frontend.app.backend
        self.datasetId = self.backend.getDatasetIds()[0]
        self.variantSetIds = [
            variantSet.getId() for variantSet in
            self.backend.getDataset(self.datasetId).getVariantSets()]

    def sendJsonPostRequest(self, path, data):
        """
        Sends a JSON request to the specified path with the specified data
        and returns the response.
        """
        return self.app.post(
            path, headers={'Content-type': 'application/json'},
            data=data)

    def sendObjectGetRequest(self, path, id_):
        """
        Sends a GET request to the specified path for an object with the
        specified ID and returns the response.
        """
        return self.app.get("{}/{}".format(path, id_))

    def testVariantSetsSearch(self):
        expectedIds = self.variantSetIds
        request = protocol.SearchVariantSetsRequest()
        request.pageSize = len(expectedIds)
        request.datasetIds = [self.datasetId]
        path = utils.applyVersion('/variantsets/search')
        response = self.sendJsonPostRequest(
            path, request.toJsonString())

        self.assertEqual(200, response.status_code)

        responseData = protocol.SearchVariantSetsResponse.fromJsonString(
            response.data)
        self.assertTrue(protocol.SearchVariantSetsResponse.validate(
            responseData.toJsonDict()))

        self.assertIsNone(responseData.nextPageToken)
        self.assertEqual(len(expectedIds), len(responseData.variantSets))
        for variantSet in responseData.variantSets:
            self.assertTrue(variantSet.id in expectedIds)

    def testGetVariantSet(self):
        path = utils.applyVersion("/variantsets")
        for variantSetId in self.variantSetIds:
            response = self.sendObjectGetRequest(path, variantSetId)
            self.assertEqual(200, response.status_code)
            responseObject = protocol.VariantSet.fromJsonString(response.data)
            self.assertEqual(responseObject.id, variantSetId)
        for badId in ["", "terribly bad ID value", "x" * 1000]:
            response = self.sendObjectGetRequest(path, badId)
            self.assertEqual(404, response.status_code)

    def testVariantsSearch(self):
        expectedIds = self.variantSetIds[:1]
        referenceName = '1'

        request = protocol.SearchVariantsRequest()
        request.referenceName = referenceName
        request.start = 0
        request.end = 0
        request.variantSetIds = expectedIds

        # Request windows is too small, no results
        path = utils.applyVersion('/variants/search')
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchVariantsResponse.fromJsonString(
            response.data)
        self.assertIsNone(responseData.nextPageToken)
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
            self.assertTrue(variant.variantSetId in expectedIds)
            self.assertEqual(variant.referenceName, referenceName)

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

        # when variantSetIds are wrong, no results
        request.variantSetIds = ["xxxx"]
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertIsNone(responseData.nextPageToken)
        self.assertEqual([], responseData.callSets)

        # if no callset name is given return all callsets
        request.variantSetIds = self.variantSetIds[:1]
        response = self.sendJsonPostRequest(
            path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertTrue(protocol.SearchCallSetsResponse.validate(
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

    def testReferences(self):
        # search for reference sets
        path = utils.applyVersion('/referencesets/search')
        request = protocol.SearchReferenceSetsRequest()
        response = self.sendJsonPostRequest(path, request.toJsonString())
        self.assertEqual(response.status_code, 200)
        responseData = protocol.SearchReferenceSetsResponse.fromJsonString(
            response.data)
        referenceSets = responseData.referenceSets
        self.assertEqual(self.numReferenceSets, len(referenceSets))

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

        for referenceSet in referenceSets:
            # fetch the reference set
            path = utils.applyVersion(
                '/referencesets/{}'.format(referenceSet.id))
            response = self.app.get(path)
            self.assertEqual(response.status_code, 200)
            fetchedReferenceSet = protocol.ReferenceSet.fromJsonString(
                response.data)
            self.assertEqual(fetchedReferenceSet, referenceSet)
            self.assertEqual(
                len(fetchedReferenceSet.referenceIds),
                self.numReferencesPerReferenceSet)

            for referenceId in referenceSet.referenceIds:
                # fetch the reference
                path = utils.applyVersion(
                    '/references/{}'.format(referenceId))
                response = self.app.get(path)
                self.assertEqual(response.status_code, 200)
                fetchedReference = protocol.Reference.fromJsonString(
                    response.data)
                self.assertEqual(fetchedReference.id, referenceId)

                # fetch the bases
                path = utils.applyVersion(
                    '/references/{}/bases'.format(referenceId))
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

    def testReads(self):
        # search read group sets
        path = utils.applyVersion('/readgroupsets/search')
        request = protocol.SearchReadGroupSetsRequest()
        request.datasetIds = ['simulatedDataset1']
        response = self.sendJsonPostRequest(path, request.toJsonString())
        self.assertEqual(response.status_code, 200)
        responseData = protocol.SearchReadGroupSetsResponse.fromJsonString(
            response.data)
        readGroupSets = responseData.readGroupSets
        self.assertEqual(len(readGroupSets), 1)

        # search reads
        path = utils.applyVersion('/reads/search')
        request = protocol.SearchReadsRequest()
        readGroupId = readGroupSets[0].readGroups[0].id
        request.readGroupIds = [readGroupId]
        response = self.sendJsonPostRequest(path, request.toJsonString())
        self.assertEqual(response.status_code, 200)
        responseData = protocol.SearchReadsResponse.fromJsonString(
            response.data)
        alignments = responseData.alignments
        self.assertEqual(len(alignments), self.numAlignmentsPerReadGroup)
        for alignment in alignments:
            self.assertEqual(alignment.readGroupId, readGroupId)
