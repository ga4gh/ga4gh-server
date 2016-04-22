"""
Sequence Annotations testing on the test data
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import logging

import ga4gh.protocol as protocol
import ga4gh.frontend as frontend
import tests.paths as paths


class TestSequenceAnnotations(unittest.TestCase):
    exampleUrl = 'www.example.com'
    datasetId = "YnJjYTE"

    @classmethod
    def setUpClass(cls):
        config = {
            "DATA_SOURCE": paths.testDataRepo,
            "DEBUG": False
        }
        logging.getLogger('ga4gh.frontend.cors').setLevel(logging.CRITICAL)
        frontend.reset()
        frontend.configure(
            baseConfig="TestConfig", extraConfig=config)
        cls.app = frontend.app.test_client()

    @classmethod
    def tearDownClass(cls):
        cls.app = None

    def sendSearchRequest(self, path, request, responseClass):
        """
        Sends the specified protocol request instance as JSON, and
        parses the result into an instance of the specified response.
        """
        response = self.sendJsonPostRequest(path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = responseClass.fromJsonString(response.data)
        self.assertTrue(responseData.validate(responseData.toJsonDict()))
        return responseData

    def getAllDatasets(self):
        path = 'datasets/search'
        request = protocol.SearchDatasetsRequest()
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchDatasetsResponse)
        return responseData.datasets

    def getAllFeatureSets(self):
        datasetId = self.getAllDatasets()[0].id
        path = 'featuresets/search'
        request = protocol.SearchFeatureSetsRequest()
        request.datasetId = datasetId
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchFeatureSetsResponse)
        return responseData.featureSets

    def testSearchFeatures(self):
        datasetId = self.getAllDatasets()[0].id
        featureSets = self.getAllFeatureSets()
        for featureSet in featureSets:
            path = "features/search"
            request = protocol.SearchFeaturesRequest()
            request.datasetId = datasetId
            request.featureSetId = featureSet.id
            request.start = 0
            request.end = 2**16
            request.featureTypes = ["exon"]
            request.referenceName = "chr1"
            request.featureSetId = featureSet.id
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            for feature in responseData.features:
                self.assertIn(
                    feature.featureType.term,
                    request.featureTypes,
                    "Term should be present {} {} \n{}\n{}".format(
                        feature.featureType.term,
                        request.featureTypes,
                        feature, request))

            path = "features/search"
            request = protocol.SearchFeaturesRequest()
            request.datasetId = datasetId
            request.featureSetId = featureSet.id
            request.start = 0
            request.end = 2**16
            request.featureTypes = ["gene", "exon"]
            request.referenceName = "chr1"
            request.featureSetId = featureSet.id
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            for feature in responseData.features:
                self.assertIn(feature.featureType.term, request.featureTypes)

            request = protocol.SearchFeaturesRequest()
            request.datasetId = datasetId
            request.featureSetId = featureSet.id
            request.start = 0
            request.end = 2**16
            request.featureTypes = ["exon"]
            request.referenceName = "chr1"
            request.featureSetId = featureSet.id
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            for feature in responseData.features:
                self.assertIn(feature.featureType.term, request.featureTypes)

    def sendJsonPostRequest(self, path, data):
        """
        Sends a JSON request to the specified path with the specified data
        and returns the response.
        """
        return self.app.post(
            path, headers={'Content-type': 'application/json'},
            data=data)
