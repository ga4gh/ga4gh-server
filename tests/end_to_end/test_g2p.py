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

    def getPhenotypeAssociationSetId(self):
        request = protocol.SearchDatasetsRequest()
        response = self.sendSearchRequest(
            "datasets/search",
            request,
            protocol.SearchDatasetsResponse)
        datasetId = response.datasets[0].id
        request = protocol.SearchPhenotypeAssociationSetsRequest()
        request.dataset_id = datasetId
        response = self.sendPostRequest("phenotypeassociationsets/search",
                                        request)
        response = protocol.fromJson(
            response.data,protocol.SearchPhenotypeAssociationSetsResponse)
        return response.phenotype_association_sets[0].id

    def sendPostRequest(self, path, request):
        """
        Sends the specified GA request object and returns the response.
        """
        headers = {
            'Content-type': 'application/json',
            'Origin': self.exampleUrl,
        }
        return self.app.post(
            path, headers=headers, data=protocol.toJson(request))

    def sendJsonPostRequest(self, path, data):
        """
        Sends a JSON request to the specified path with the specified data
        and returns the response.
        """
        return self.app.post(
            path, headers={'Content-type': 'application/json'},
            data=data)

    def testPhenotypeAssociationSetSearch(self):
        request = protocol.SearchDatasetsRequest()
        response = self.sendSearchRequest(
            "datasets/search",
            request,
            protocol.SearchDatasetsResponse)
        datasetId = response.datasets[0].id
        request = protocol.SearchPhenotypeAssociationSetsRequest()
        request.dataset_id = datasetId
        response = self.sendSearchRequest(
            "phenotypeassociationsets/search",
            request,
            protocol.SearchPhenotypeAssociationSetsResponse)
        # there should be an array
        self.assertIsNotNone(response.phenotype_association_sets)
        # there should be at least one entry
        self.assertGreater(len(response.phenotype_association_sets), 0)

    def testGenotypesSearchByExternalIdentifier(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        # setup the external identifiers query
        extid = protocol.ExternalIdentifier()
        # http://www.ncbi.nlm.nih.gov/SNP/121908585
        extid.identifier = "121908585"
        extid.version = "*"
        extid.database = "dbSNP"
        request.external_identifiers.extend([extid])
        response = self.sendSearchRequest(
            '/genotypes/search',
            request,
            protocol.SearchGenotypesResponse)
        self.assertEqual(1, len(response.genotypes))

    def testFindFeatureExternalIdentifier(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        # setup the external identifiers query
        extid = protocol.ExternalIdentifier()
        # http://www.ncbi.nlm.nih.gov/SNP/121908585
        extid.identifier = "121908585"
        extid.version = "*"
        extid.database = "dbSNP"
        request.external_identifiers.extend([extid])
        response = self.sendSearchRequest(
            '/genotypes/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations))
        genotypeId = response.associations[0].id

        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.genotype_ids = [genotypeId]
        response = self.sendSearchRequest(
            '/genotypephenotypes/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations[0].features))

    def testGenotypesSearchById(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.id = \
            "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=965"
        postUrl = '/genotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchGenotypesResponse)
        self.assertEqual(request.id, response.genotypes[0].id)

    def testGenotypesSearchByReferenceName(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.reference_name = "RET M918T missense mutation"
        postUrl = '/genotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchGenotypesResponse)
        self.assertEqual(
            "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=965",
            response.genotypes[0].id)

    def testGenotypesSearchByReferenceNameKIT(self):
        request = protocol.SearchGenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.reference_name = \
            "KIT *wild"
        postUrl = '/genotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchGenotypesResponse)
        self.assertEqual(
            "http://ohsu.edu/cgd/27d2169c",
            response.genotypes[0].id)

    def testPhenotypesSearchById(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.id = "http://ohsu.edu/cgd/30ebfd1a"
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertEqual(request.id, response.phenotypes[0].id)

    def testPhenotypesSearchOntologyTerm(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.type.id = "http://ohsu.edu/cgd/5c895709"
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypeSearchQualifiersSensitivity(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://ohsu.edu/cgd/sensitivity"
        request.qualifiers.extend([ontologyterm])
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypeSearchQualifiersSensitivityPATO_0000396(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://purl.obolibrary.org/obo/PATO_0000396"
        request.qualifiers.extend([ontologyterm])
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypeSearchMultipleQualifiers(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.id = "http://purl.obolibrary.org/obo/PATO_0000396"
        ontologyterm2 = protocol.OntologyTerm()
        ontologyterm2.id = "http://purl.obolibrary.org/obo/PATO_0000396"
        request.qualifiers.extend([ontologyterm, ontologyterm2])
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertGreater(len(response.phenotypes), 0)

    @unittest.skip
    def testPhenotypesSearchDescription(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.description = \
                "Papillary thyroid carcinoma with sensitivity to therapy"  # noqa
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypesSearchDescriptionWildcard(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.description = ".*sensitivity.*"
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertEquals(7, len(response.phenotypes))

    def testPhenotypesSearchMultipleTerms(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.description = "Melanoma, NOS with response to therapy"
        ontologyterm = protocol.OntologyTerm()
        request.age_of_on_set.id = "http://purl.obolibrary.org/obo/HP_0003581"
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertGreater(len(response.phenotypes), 0)

    def testGenotypePhenotypeSearchFeature(self):
        """
        Search for associations given a feature
        """
        # simple string regexp
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.genotype_ids.extend(["http://ohsu.edu/cgd/27d2169c"])
        response = self.sendSearchRequest(
            '/genotypephenotypes/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations[0].features))

    def testGenotypePhenotypeSearchEvidence(self):
        """
        Search for associations given an evidence
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        eq = protocol.EvidenceQuery()
        eq.description = "imatinib"
        request.evidence.extend([eq])
        response = self.sendSearchRequest(
            '/genotypephenotypes/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations[0].features))

    def testGenotypePhenotypeSearchPhenotype(self):
        """
        Search for associations given a phenotype
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.phenotype_ids.extend(["http://ohsu.edu/cgd/25abbb09"])
        response = self.sendSearchRequest(
            '/genotypephenotypes/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations[0].features))

    def testNoFind(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.feature_ids.extend(["FOOBAR"])
        response = self.sendSearchRequest(
            '/genotypephenotypes/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(0, len(response.associations))

    def testGenotypePhenotypeSearchEnsureEvidence(self):
        """
        Ensure evidence level is serialized in responses
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.feature_ids.extend(["http://ohsu.edu/cgd/27d2169c"])
        response = self.sendSearchRequest(
            '/genotypephenotypes/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations[0].evidence))
        evidence = response.associations[0].evidence[0]
        self.assertEqual('decreased_sensitivity', evidence.description)

    def testGenotypePhenotypeSearchEnsureEnvironment(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        request.feature_ids.extend(["http://ohsu.edu/cgd/27d2169c"])
        eq = protocol.EvidenceQuery()
        eq.description = "imatinib"
        request.evidence.extend([eq])
        response = self.sendSearchRequest(
            '/genotypephenotypes/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(
            1, len(response.associations[0].environmental_contexts))
        environmentalContext = response.associations[0].\
                               environmental_contexts[0]
        self.assertEqual('imatinib', environmentalContext.description)

    def testGenotypeSearchFeaturePagingOne(self):
        """
        If page size is set to 1 only one association should be returned
        """
        request = protocol.SearchGenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.page_size = 1
        request.reference_name = \
            "KIT *wild"
        postUrl = '/genotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchGenotypesResponse)
        self.assertEqual(1, len(response.genotypes))
        self.assertIsNotNone(response.nextPageToken)

    def testGenotypeSearchFeaturePagingMore(self):
        """
        If page size is not set to more than one association should be returned
        """
        request = protocol.SearchGenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.reference_name = \
            "KIT *wild"
        postUrl = '/genotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchGenotypesResponse)
        self.assertGreater(len(response.genotypes), 1)
        self.assertIsNone(response.nextPageToken)

    def testGenotypeSearchFeaturePagingAll(self):
        """
        Loop through all pages
        """
        request = protocol.SearchGenotypesRequest()
        request.phenotype_association_set_id = self.getPhenotypeAssociationSetId()
        # setup phenotype query
        request.page_size = 1
        request.reference_name = \
            "KIT *wild"
        postUrl = '/genotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchGenotypesResponse)
        self.assertEqual(1, len(response.genotypes))
        self.assertIsNotNone(response.next_page_token)
        pageCount = 1
        while response.nextPageToken:
            previous_id = response.genotypes[0].id
            request = protocol.SearchGenotypesRequest()
            request.phenotype_association_set_id =\
                self.getPhenotypeAssociationSetId()
            request.page_token = response.next_page_token
            request.page_size = 1
            request.reference_name = "KIT *wild"
            response = self.sendSearchRequest(
                postUrl,
                request,
                protocol.SearchGenotypesResponse)
            self.assertEqual(1, len(response.genotypes))
            self.assertNotEqual(previous_id, response.genotypes[0].id)
            pageCount += 1
        self.assertEqual(3, pageCount)
