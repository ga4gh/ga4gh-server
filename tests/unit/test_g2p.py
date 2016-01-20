"""
Unit tests for genotypephenotype objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import logging

import ga4gh.frontend as frontend
import ga4gh.protocol as protocol


class TestGenotypePhenotypeSearch(unittest.TestCase):
    """
    Tests both front end and backend
    genotypephenotype/search POST requests.
    """
    exampleUrl = 'www.example.com'

    @classmethod
    def setUpClass(cls):
        config = {
            "DATA_SOURCE": "simulated://",
            "SIMULATED_BACKEND_RANDOM_SEED": 1111,
            "SIMULATED_BACKEND_NUM_CALLS": 1,
            "SIMULATED_BACKEND_VARIANT_DENSITY": 1.0,
            "SIMULATED_BACKEND_NUM_VARIANT_SETS": 1,
            # "DEBUG" : True
        }
        frontend.reset()
        frontend.configure(
            baseConfig="TestConfig", extraConfig=config)
        cls.app = frontend.app.test_client()
        # silence usually unhelpful CORS log
        logging.getLogger('ga4gh.frontend.cors').setLevel(logging.CRITICAL)
        # example test values
        cls.backend = frontend.app.backend

    def sendPostRequest(self, path, request):
        """
        Sends the specified GA request object and returns the response.
        """
        headers = {
            'Content-type': 'application/json',
            'Origin': self.exampleUrl,
        }
        return self.app.post(
            path, headers=headers, data=request.toJsonString())

    def testGenotypePhenotypeSearchFeature(self):
        """
        Search for evidence on a genomic feature given feature name
        """
        # simple string regexp
        request = protocol.SearchGenotypePhenotypeRequest()

        request.feature = "KIT *wild"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

        request.evidence = "imatinib"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

        request.phenotype = "GIST"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))
        print(response.toJsonString())

        request.phenotype = "FOOBAR"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(0, len(response.associations))

        # identifiers
        request = protocol.SearchGenotypePhenotypeRequest()
        request.feature = protocol.ExternalIdentifierQuery()
        id = protocol.ExternalIdentifier()
        id.database = "http://ohsu.edu/cgd/"
        id.identifier = "4841bf74"
        id.version = "*"
        request.feature.ids = [id]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

        request.phenotype = protocol.ExternalIdentifierQuery()
        id = protocol.ExternalIdentifier()
        id.database = "http://ohsu.edu/cgd/"
        id.identifier = "37da8697"
        id.version = "*"
        request.phenotype.ids = [id]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

        request.evidence = protocol.ExternalIdentifierQuery()
        id = protocol.ExternalIdentifier()
        id.database = "http://www.drugbank.ca/drugs/"
        id.identifier = "DB00398"
        id.version = "*"
        request.evidence.ids = [id]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertNotEqual(0, len(response.associations))
        self.assertEqual(1, len(response.associations[0].features))

        request.evidence = protocol.ExternalIdentifierQuery()
        id = protocol.ExternalIdentifier()
        id.database = "FOO"
        id.identifier = "DB00619"
        id.version = "*"
        request.evidence.ids = [id]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(0, len(response.associations))
        self.testGenotypePhenotypeSearchFeaturePagingOne()
        self.testGenotypePhenotypeSearchFeaturePagingMore()

    def testGenotypePhenotypeSearchFeaturePagingOne(self):
        """
        If page size is set to 1 only one association should be returned
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.pageSize = 1
        request.feature = "KIT *wild"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)

        self.assertEqual(1, len(response.associations))
        self.assertEqual(1, len(response.associations[0].features))
        self.assertIsNotNone(response.nextPageToken)

    def testGenotypePhenotypeSearchFeaturePagingMore(self):
        """
        If page size is not set to more than one association should be returned
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.feature = "KIT *wild"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)

        self.assertGreater(len(response.associations), 1)
        self.assertIsNone(response.nextPageToken)

    def testGenotypePhenotypeSearchFeaturePagingAll(self):
        """
        Loop through all pages
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.pageSize = 1
        request.feature = "KIT *wild"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations))
        self.assertIsNotNone(response.nextPageToken)

        for i in range(3):
            previous_id = response.associations[0].id
            request = protocol.SearchGenotypePhenotypeRequest()
            request.pageToken = response.nextPageToken
            request.pageSize = 1
            request.feature = "KIT *wild"
            response = self.sendPostRequest(
                         '/genotypephenotype/search', request)
            self.assertEqual(200, response.status_code)
            response = protocol.SearchGenotypePhenotypeResponse().\
                fromJsonString(response.data)
            self.assertEqual(1, len(response.associations))
            self.assertNotEqual(previous_id, response.associations[0].id)
            if i != 2:
                self.assertIsNotNone(response.nextPageToken)
        # from IPython.core.debugger import Pdb ;        Pdb().set_trace()

    def testGenotypePheontypeSearchEnsureEvidenceLevel(self):
        """
        Ensure evidence level is serialized in responses
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.feature = "KIT *wild"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertTrue(hasattr(response.associations[0],
                                'evidence'))
        sample_evidence = response.associations[0].toJsonDict()['evidence'][0]
        sample_evidence_type = sample_evidence['evidenceType']
        self.assertIn('ontologySource', sample_evidence_type)
        self.assertEqual(sample_evidence_type['ontologySource'],
                         'http://ohsu.edu/cgd/')
        self.assertIn('id', sample_evidence_type)
        self.assertEqual(sample_evidence_type['id'], 'c703f7ab')
        self.assertIn('name', sample_evidence_type)
        self.assertEqual(sample_evidence_type['name'], 'early trials')
