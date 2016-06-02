"""
G2P testing on the test data
"""

import unittest

import ga4gh.protocol as protocol
import ga4gh.frontend as frontend
import tests.paths as paths


class TestG2P(unittest.TestCase):
    exampleUrl = 'www.example.com'
    phenotypeAssociationSetId = ""

    @classmethod
    def setUpClass(cls):
        config = {
            "DATA_SOURCE": paths.testDataRepo,
            # "DEBUG": True
        }
        frontend.reset()
        frontend.configure(
            baseConfig="DevelopmentConfig", extraConfig=config)
        cls.app = frontend.app.test_client()

    def getPhenotypeAssociationSetId(self):
        request = protocol.SearchDatasetsRequest()
        response = self.sendPostRequest("datasets/search", request)
        response = protocol.SearchDatasetsResponse().fromJsonString(
             response.data)
        datasetId = response.datasets[0].id
        request = protocol.SearchPhenotypeAssociationSetsRequest()
        request.datasetId = datasetId
        response = self.sendPostRequest("phenotypeassociationsets/search",
                                        request)
        response = protocol.SearchPhenotypeAssociationSetsResponse(
            ).fromJsonString(response.data)
        return response.phenotypeAssociationSets[0].id

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

    def testPhenotypeAssociationSetSearch(self):
        request = protocol.SearchDatasetsRequest()
        response = self.sendPostRequest("datasets/search", request)
        response = protocol.SearchDatasetsResponse().fromJsonString(
             response.data)
        datasetId = response.datasets[0].id
        request = protocol.SearchPhenotypeAssociationSetsRequest()
        request.datasetId = datasetId
        response = self.sendPostRequest("phenotypeassociationsets/search",
                                        request)
        response = protocol.SearchPhenotypeAssociationSetsResponse(
            ).fromJsonString(response.data)
        self.assertIsNotNone(response.phenotypeAssociationSets)

    def testGenotypePhenotypeSearchFeature(self):
        """
        Search for evidence on a genomic feature given feature name
        """
        # simple string regexp
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
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

        # TODO search with protocol.PhenotypeQuery
        pass

    def testNoFind(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.phenotype = "FOOBAR"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(0, len(response.associations))

    def testFindEvidenceExternalIdentifier(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.evidence = protocol.ExternalIdentifierQuery()
        id = protocol.ExternalIdentifier()
        id.database = "http://www.drugbank.ca/drugs/"
        id.identifier = "DB00619"
        id.version = "*"
        request.evidence.ids = [id]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertNotEqual(0, len(response.associations))
        self.assertEqual(1, len(response.associations[0].features))

    def testFindFeatureExternalIdentifier(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.feature = protocol.ExternalIdentifierQuery()
        id = protocol.ExternalIdentifier()
        id.database = "http://ohsu.edu/cgd/"
        id.identifier = "055b872c"
        id.version = "*"
        request.feature.ids = [id]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

    def testFindPhenotypeExternalIdentifier(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.phenotype = protocol.ExternalIdentifierQuery()
        id = protocol.ExternalIdentifier()
        id.database = "http://ohsu.edu/cgd/"
        id.identifier = "032c97e8"
        id.version = "*"
        request.phenotype.ids = [id]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

    def testGenotypePhenotypeSearchFeaturePagingOne(self):
        """
        If page size is set to 1 only one association should be returned
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
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
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
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
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.pageSize = 1
        request.feature = "KIT *wild"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations))
        self.assertIsNotNone(response.nextPageToken)
        pageCount = 1
        while response.nextPageToken:
            previous_id = response.associations[0].id
            request = protocol.SearchGenotypePhenotypeRequest()
            request.phenotypeAssociationSetId =\
                self.getPhenotypeAssociationSetId()
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
            pageCount += 1
        self.assertEqual(3, pageCount)

    def testGenotypePhenotypeSearchEnsureEvidenceLevel(self):
        """
        Ensure evidence level is serialized in responses
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.feature = "KIT *wild"
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertTrue(hasattr(response.associations[0],
                                'evidence'))
        sample_evidence = response.associations[0].toJsonDict()['evidence'][0]
        sample_evidence_type = sample_evidence['evidenceType']
        self.assertIn('sourceName', sample_evidence_type)
        self.assertEqual(sample_evidence_type['sourceName'],
                         'CGD')
        self.assertIn('id', sample_evidence_type)
        self.assertEqual(sample_evidence_type['id'],
                         'http://ohsu.edu/cgd/87752f6c')
        self.assertIn('term', sample_evidence_type)
        self.assertEqual(sample_evidence_type['term'],
                         'http://purl.obolibrary.org/obo/ECO_0000033')
        self.assertEqual(sample_evidence['description'],
                         'decreased_sensitivity')

    def testGenotypePheontypeSearchOntologyTermPrefixTerm(self):
        """
        ensure we can read ontology terms usng term
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()

        term = protocol.OntologyTerm()
        term.term = "DrugBank:DB01268"
        request.evidence = protocol.OntologyTermQuery()
        request.evidence.terms = [term]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertTrue(hasattr(response.associations[0],
                                'evidence'))

    def testGenotypePheontypeSearchEnsureOntologyTermFeature(self):
        """
        Ensure evidence level is serialized in responses
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        term = protocol.OntologyTerm()
        term.term = "CGD:27d2169c"
        request.feature = protocol.OntologyTermQuery()
        request.feature.terms = [term]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertTrue(hasattr(response.associations[0],
                                'evidence'))
