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

    def testFeaturesSearch(self):
        request = protocol.SearchFeaturesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup the external identifiers query
        idquery = protocol.ExternalIdentifierQuery()
        extid = protocol.ExternalIdentifier()
        extid.identifier = "rs6920220"
        extid.version = "*"
        extid.database = "dbSNP"
        idquery.ids = [extid]
        # setup the term query
        termquery = protocol.TermQuery()
        termquery.term = idquery
        request.termQueries = [termquery]
        response = self.sendPostRequest('features/search', request)
        self.assertEqual(200, response.status_code)

    def testPhenotypesSearch(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        phenotypeQuery = protocol.PhenotypeQuery()
        phenotypeQuery.id = "p12345"
        request.phenotype = phenotypeQuery
        postUrl = 'associations/%s/phenotypes/search' % \
                  request.phenotypeAssociationSetId
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(0, len(response.phenotypes))
        self.assertEqual("p12345", response.phenotypes[0].id)

    def testPhenotypesSearchOntologyTerm(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        phenotypeQuery = protocol.PhenotypeQuery()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://www.ebi.ac.uk/efo/EFO_0003767"
        phenotypeQuery.type = ontologyterm
        postUrl = 'associations/%s/phenotypes/search' % \
                  request.phenotypeAssociationSetId
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(0, response.phenotypes)

    def testPhenotypeSearchQualifiers(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        phenotypeQuery = protocol.PhenotypeQuery()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://purl.obolibrary.org/obo/PATO_0001899"
        phenotypeQuery.qualifiers = [ontologyterm]
        postUrl = 'associations/%s/phenotypes/search' % \
                  request.phenotypeAssociationSetId
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(0, response.phenotypes)

    def testPhenotypeSearchMultipleQualifiers(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        phenotypeQuery = protocol.PhenotypeQuery()
        ontologyterm1 = protocol.OntologyTerm()
        ontologyterm1.id = "http://purl.obolibrary.org/obo/PATO_0000396"
        ontologyterm2 = protocol.OntologyTerm()
        ontologyterm2.id = "http://purl.obolibrary.org/obo/PATO_0000460"
        phenotypeQuery.qualifiers = [ontologyterm1, ontologyterm2]
        postUrl = 'associations/%s/phenotypes/search' % \
                  request.phenotypeAssociationSetId
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(0, response.phenotypes)

    def testPhenotypesSearchDescription(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        phenotypeQuery = protocol.PhenotypeQuery()
        phenotypeQuery.description = "inflammatory bowel disease"
        postUrl = 'associations/%s/phenotypes/search' % \
                  request.phenotypeAssociationSetId
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(0, response.phenotypes)

    def testPhenotypesSearchMultipleTerms(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        phenotypeQuery = protocol.PhenotypeQuery()
        phenotypeQuery.description = "AML"
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://purl.obolibrary.org/obo/HP_0003581"
        phenotypeQuery.ageOfOnset = ontologyterm
        postUrl = 'associations/%s/phenotypes/search' % \
                  request.phenotypeAssociationSetId
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(0, response.phenotypes)

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
        request = protocol.SearchFeaturesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup the external identifiers query
        idquery = protocol.ExternalIdentifierQuery()
        extid = protocol.ExternalIdentifier()
        extid.identifier = "rs6920220"
        extid.version = "*"
        extid.database = "dbSNP"
        idquery.ids = [extid]
        # setup the term query
        termquery = protocol.TermQuery()
        termquery.term = idquery
        request.termQueries = [termquery]
        response = self.sendPostRequest('features/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchFeaturesResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.features))

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
