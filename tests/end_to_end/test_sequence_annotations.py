"""
Sequence Annotations testing on the test data
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import logging

import ga4gh.server.frontend as frontend
import tests.paths as paths

import ga4gh.schemas.protocol as protocol


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
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(response.data, responseClass)
        self.assertTrue(
            protocol.validate(protocol.toJson(responseData), responseClass))
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
        request.dataset_id = datasetId
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchFeatureSetsResponse)
        return responseData.feature_sets

    def testSearchFeaturesByName(self):
        ran = False
        featureSets = self.getAllFeatureSets()
        for featureSet in featureSets:
            path = "features/search"
            request = protocol.SearchFeaturesRequest()
            request.feature_set_id = featureSet.id
            request.name = "BAD NAME"
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            self.assertEqual(0, len(responseData.features))
            request.name = "exon:ENSTR0000507418.3:5"
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            for feature in responseData.features:
                ran = True
                self.assertEqual(feature.name, request.name)
        self.assertTrue(ran)

    def testSearchFeaturesByGeneSymbol(self):
        ran = False
        featureSets = self.getAllFeatureSets()
        for featureSet in featureSets:
            path = "features/search"
            request = protocol.SearchFeaturesRequest()
            request.feature_set_id = featureSet.id
            request.gene_symbol = "BAD GENE SYMBOL"
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            self.assertEqual(0, len(responseData.features))
            request.gene_symbol = "DDX11L16"
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            for feature in responseData.features:
                ran = True
                self.assertEqual(feature.gene_symbol, request.gene_symbol)
        self.assertTrue(ran)

    def testSearchFeatures(self):
        featureSets = self.getAllFeatureSets()
        for featureSet in featureSets:
            path = "features/search"
            request = protocol.SearchFeaturesRequest()
            request.feature_set_id = featureSet.id
            request.start = 0
            request.end = 2**16
            request.feature_types.extend(["exon"])
            request.reference_name = "chr1"
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            for feature in responseData.features:
                self.assertIn(
                    feature.feature_type.term,
                    request.feature_types,
                    "Term should be present {} {} \n{}\n{}".format(
                        feature.feature_type.term,
                        request.feature_types,
                        feature, request))

            path = "features/search"
            request = protocol.SearchFeaturesRequest()
            request.feature_set_id = featureSet.id
            request.start = 0
            request.end = 2**16
            request.feature_types.extend(["gene", "exon"])
            request.reference_name = "chr1"
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            for feature in responseData.features:
                self.assertIn(feature.feature_type.term, request.feature_types)
            request = protocol.SearchFeaturesRequest()
            request.feature_set_id = featureSet.id
            request.start = 0
            request.end = 2**16
            request.feature_types.extend(["exon"])
            request.reference_name = "chr1"
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchFeaturesResponse)
            for feature in responseData.features:
                self.assertIn(feature.feature_type.term, request.feature_types)

    def sendJsonPostRequest(self, path, data):
        """
        Sends a JSON request to the specified path with the specified data
        and returns the response.
        """
        return self.app.post(
            path, headers={'Content-type': 'application/json'},
            data=data)

    def getAllContinuousSets(self):
        datasetId = self.getAllDatasets()[0].id
        path = 'continuoussets/search'
        request = protocol.SearchContinuousSetsRequest()
        request.dataset_id = datasetId
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchContinuousSetsResponse)
        return responseData.continuous_sets

    def testSearchContinuous(self):
        continuousSets = self.getAllContinuousSets()
        for continuousSet in continuousSets:
            path = "continuous/search"
            request = protocol.SearchContinuousRequest()
            request.continuous_set_id = continuousSet.id
            request.start = 49200000
            request.end = 49308000
            request.reference_name = "chr19"
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchContinuousResponse)
            for continuous in responseData.continuous:
                self.assertGreater(len(continuous.values), 0)
