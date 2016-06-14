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
            "DEBUG": False
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
        # there should be an array
        self.assertIsNotNone(response.phenotypeAssociationSets)
        # there should be at least one entry
        self.assertGreater(len(response.phenotypeAssociationSets), 0)

    def testGenotypesSearchByExternalIdentifier(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup the external identifiers query
        extid = protocol.ExternalIdentifier()
        # http://www.ncbi.nlm.nih.gov/SNP/121908585
        extid.identifier = "121908585"
        extid.version = "*"
        extid.database = "dbSNP"
        request.externalIdentifiers = [extid]
        response = self.sendPostRequest('/genotypes/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertEqual(1, len(response.genotypes))

    def testFindFeatureExternalIdentifier(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup the external identifiers query
        extid = protocol.ExternalIdentifier()
        # http://www.ncbi.nlm.nih.gov/SNP/121908585
        extid.identifier = "121908585"
        extid.version = "*"
        extid.database = "dbSNP"
        request.externalIdentifiers = [extid]
        response = self.sendPostRequest('/genotypes/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertEqual(1, len(response.genotypes))
        genotypeId = response.genotypes[0].id

        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.genotypeIds = [genotypeId]
        response = self.sendPostRequest('/genotypephenotypes/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

    def testGenotypesSearchById(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.id = \
            "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=965"
        postUrl = '/genotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertEqual(request.id, response.genotypes[0].id)

    def testGenotypesSearchByReferenceName(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.referenceName = "RET M918T missense mutation"
        postUrl = '/genotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertEqual(
            "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=965",
            response.genotypes[0].id)

    def testGenotypesSearchByReferenceNameKIT(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.id = "http://ohsu.edu/cgd/30ebfd1a"
        postUrl = '/phenotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        print(response.data)
        response = protocol.SearchPhenotypeResponse() \
                           .fromJsonString(response.data)
        self.assertEqual(
            "http://ohsu.edu/cgd/27d2169c",
            response.genotypes[0].id)

    def testPhenotypesSearchById(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.id = "http://ohsu.edu/cgd/30ebfd1a"
        postUrl = '/phenotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse()\
                           .fromJsonString(response.data)
        self.assertEqual(request.id, response.phenotypes[0].id)

    def testPhenotypesSearchOntologyTerm(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        phenotypeQuery = protocol.PhenotypeQuery()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://ohsu.edu/cgd/5c895709"
        request.type = ontologyterm
        postUrl = '/phenotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypeSearchQualifiersSensitivity(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        phenotypeQuery = protocol.PhenotypeQuery()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://ohsu.edu/cgd/sensitivity"
        request.qualifiers = [ontologyterm]
        postUrl = '/phenotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypeSearchQualifiersSensitivityPATO_0000396(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://purl.obolibrary.org/obo/PATO_0000396"
        request.qualifiers = [ontologyterm]
        postUrl = '/phenotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypeSearchMultipleQualifiers(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://purl.obolibrary.org/obo/PATO_0000396"
        ontologyterm2 = protocol.OntologyTerm()
        ontologyterm2.id = "http://purl.obolibrary.org/obo/PATO_0000396"
        request.qualifiers = [ontologyterm, ontologyterm2]
        postUrl = '/phenotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(len(response.phenotypes), 0)

    @unittest.skip
    def testPhenotypesSearchDescription(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.description = \
                "Papillary thyroid carcinoma with sensitivity to therapy"  # noqa
        postUrl = '/phenotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypesSearchDescriptionWildcard(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.description = ".*sensitivity.*"
        postUrl = '/phenotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertEquals(7, len(response.phenotypes))

    def testPhenotypesSearchMultipleTerms(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.description = "Melanoma, NOS with response to therapy"
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://purl.obolibrary.org/obo/HP_0003581"
        phenotypeQuery.ageOfOnset = ontologyterm
        postUrl = 'associations/%s/phenotypes/search' % \
                  request.phenotypeAssociationSetId
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchPhenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(len(response.phenotypes), 0)

    def testGenotypePhenotypeSearchFeature(self):
        """
        Search for associations given a feature
        """
        # simple string regexp
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.genotypeIds = ["http://ohsu.edu/cgd/27d2169c"]
        response = self.sendPostRequest('/genotypephenotypes/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

    def testGenotypePhenotypeSearchEvidence(self):
        """
        Search for associations given an evidence
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        eq = protocol.EvidenceQuery()
        eq.description = "imatinib"
        request.evidence = [eq]
        response = self.sendPostRequest('/genotypephenotypes/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

    def testGenotypePhenotypeSearchPhenotype(self):
        """
        Search for associations given a phenotype
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.phenotypeIds = ["http://ohsu.edu/cgd/25abbb09"]
        response = self.sendPostRequest('/genotypephenotypes/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

    def testNoFind(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.genotypeIds = ["FOOBAR"]
        response = self.sendPostRequest('/genotypephenotypes/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(0, len(response.associations))

    def testGenotypePhenotypeSearchEnsureEvidence(self):
        """
        Ensure evidence level is serialized in responses
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.genotypeIds = ["http://ohsu.edu/cgd/27d2169c"]
        response = self.sendPostRequest('/genotypephenotypes/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)

        self.assertEqual(1, len(response.associations[0].evidence))
        evidence = response.associations[0].evidence[0]
        self.assertEqual('decreased_sensitivity', evidence.description)

    def testGenotypePhenotypeSearchEnsureEnvironment(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        request.genotypeIds = ["http://ohsu.edu/cgd/27d2169c"]
        eq = protocol.EvidenceQuery()
        eq.description = "imatinib"
        request.evidence = [eq]
        response = self.sendPostRequest('/genotypephenotypes/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].environmentalContexts))
        environmentalContext = response.associations[0].environmentalContexts[0]
        self.assertEqual('imatinib', environmentalContext.description)

    def testGenotypeSearchFeaturePagingOne(self):
        """
        If page size is set to 1 only one association should be returned
        """
        request = protocol.SearchGenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.pageSize = 1
        request.referenceName = \
            "KIT *wild"
        postUrl = '/genotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertEqual(1, len(response.genotypes))
        self.assertIsNotNone(response.nextPageToken)

    def testGenotypeSearchFeaturePagingMore(self):
        """
        If page size is not set to more than one association should be returned
        """
        request = protocol.SearchGenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.referenceName = \
            "KIT *wild"
        postUrl = '/genotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertGreater(len(response.genotypes), 1)
        self.assertIsNone(response.nextPageToken)

    def testGenotypeSearchFeaturePagingAll(self):
        """
        Loop through all pages
        """
        request = protocol.SearchGenotypesRequest()
        request.phenotypeAssociationSetId = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.pageSize = 1
        request.referenceName = \
            "KIT *wild"
        postUrl = '/genotypes/search'
        response = self.sendPostRequest(postUrl, request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypesResponse() \
                           .fromJsonString(response.data)
        self.assertEqual(1, len(response.genotypes))
        self.assertIsNotNone(response.nextPageToken)
        pageCount = 1
        while response.nextPageToken:
            previous_id = response.genotypes[0].id
            request = protocol.SearchGenotypesRequest()
            request.phenotypeAssociationSetId =\
                self.getPhenotypeAssociationSetId()
            request.pageToken = response.nextPageToken
            request.pageSize = 1
            request.referenceName = "KIT *wild"
            response = self.sendPostRequest(postUrl, request)
            self.assertEqual(200, response.status_code)
            response = protocol.SearchGenotypesResponse().\
                fromJsonString(response.data)
            self.assertEqual(1, len(response.genotypes))
            self.assertNotEqual(previous_id, response.genotypes[0].id)
            pageCount += 1
        self.assertEqual(3, pageCount)
