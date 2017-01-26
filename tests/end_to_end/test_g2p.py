"""
G2P testing on the test data
"""
import unittest

import ga4gh.server.datamodel as datamodel
import ga4gh.server.frontend as frontend
import tests.paths as paths

import ga4gh.schemas.protocol as protocol


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

    def sendGetRequest(self, path):
        """
        Sends a get request to the specified URL and returns the response.
        """
        return self.app.get(path)

    def getPhenotypeAssociationSetId(self):
        """
        Gets the dataset phenotype association set ID
        """
        request = protocol.SearchDatasetsRequest()
        response = self.sendSearchRequest(
            "datasets/search",
            request,
            protocol.SearchDatasetsResponse)
        datasetId = response.datasets[0].id
        request = protocol.SearchPhenotypeAssociationSetsRequest()
        request.dataset_id = datasetId
        response = self.sendPostRequest(
            "phenotypeassociationsets/search", request)
        response = protocol.fromJson(
            response.data, protocol.SearchPhenotypeAssociationSetsResponse)
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

    def getAllDatasets(self):
        """
        Gets all datasets available
        """
        path = 'datasets/search'
        request = protocol.SearchDatasetsRequest()
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchDatasetsResponse)
        return responseData.datasets

    def getAllFeatureSets(self):
        """
        Gets all feature sets available
        """
        datasetId = self.getAllDatasets()[0].id
        datasetName = self.getAllDatasets()[0].name
        path = 'featuresets/search'
        request = protocol.SearchFeatureSetsRequest()
        request.dataset_id = datasetId
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchFeatureSetsResponse)
        return datasetName, responseData.feature_sets

    def getCGDDataSetFeatureSet(self):
        """
        Gets CGD data feature set
        """
        datasetName, featureSets = self.getAllFeatureSets()
        for featureSet in featureSets:
            if featureSet.name == 'cgd':
                return datasetName, featureSet

    def getObfuscatedFeatureCompoundId(
            self, dataSetName, featureSetName, featureId):
        """
        Gets the obfuscated feature compound Id
        """
        splits = [
            dataSetName,
            featureSetName,
            featureId]
        joined = datamodel.FeatureSetCompoundId.join(splits)
        obfuscated = datamodel.FeatureCompoundId.obfuscate(joined)
        return obfuscated

    def testEnsureCGDFeatureSet(self):
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        self.assertIsNotNone(datasetName)
        self.assertIsNotNone(featureSet)
        self.assertIsNotNone(featureSet.name)

    def testEnsureCGDFeatureId(self):
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        featureId = \
            "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=736"
        obfuscated = self.getObfuscatedFeatureCompoundId(
            datasetName, featureSet.name, featureId)
        compoundId = datamodel.FeatureCompoundId.parse(obfuscated)
        self.assertEqual(featureId, compoundId.featureId)

    def testCompoundFeatureSearch(self):
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        featureId = \
            "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=736"
        obfuscated = self.getObfuscatedFeatureCompoundId(
            datasetName, featureSet.name, featureId)
        request = protocol.GetFeatureRequest
        request.feature_id = obfuscated
        response = self.sendGetRequest('/features/{}'.format(obfuscated))

        feature = protocol.fromJson(response.data, protocol.Feature)

        self.assertIsNotNone(feature)
        featureId = feature.id

        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        request.feature_ids.append(featureId)
        response = self.sendSearchRequest(
            '/featurephenotypeassociations/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations))
        self.assertEqual(1, len(response.associations[0].feature_ids))

    def testFeaturesSearchById(self):
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        featureId = \
            "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=965"
        obfuscated = self.getObfuscatedFeatureCompoundId(
            datasetName, featureSet.name, featureId)
        request = protocol.GetFeatureRequest
        request.feature_id = obfuscated
        response = self.sendGetRequest(
            '/features/{}'.format(obfuscated))

        feature = protocol.fromJson(response.data, protocol.Feature)
        self.assertIsNotNone(feature)
        self.assertEqual(request.feature_id, feature.id)
        self.assertIsNotNone(feature.feature_type)
        self.assertIsNotNone(feature.feature_type.term_id)
        self.assertEqual(feature.reference_name,  "chr10")
        self.assertEqual(feature.start,  43617416)
        self.assertEqual(feature.end,  43617416)

    def testGenotypesSearchByName(self):
        # setup phenotype query
        request = protocol.SearchFeaturesRequest()
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        request.feature_set_id = featureSet.id
        request.name = "RET M918T missense mutation"

        postUrl = "features/search"
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchFeaturesResponse)
        self.assertEqual(1, len(response.features))
        self.assertEqual(
            "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=965",
            datamodel.FeatureCompoundId
            .parse(response.features[0].id)
            .featureId
            )
        self.assertEqual(
            request.name,
            response.features[0].name
            )

    def testGenotypesSearchByNameKIT(self):
        request = protocol.SearchFeaturesRequest()
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        request.feature_set_id = featureSet.id
        request.name = \
            "KIT *wild"
        postUrl = "features/search"
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchFeaturesResponse)
        self.assertEqual(3, len(response.features))

    def testPhenotypesSearchById(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
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
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        request.type.term_id = "http://ohsu.edu/cgd/5c895709"
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypeSearchQualifiersSensitivity(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.term_id = "http://ohsu.edu/cgd/sensitivity"
        request.qualifiers.extend([ontologyterm])
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypeSearchQualifiersSensitivityPATO_0000396(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.term_id = "http://purl.obolibrary.org/obo/PATO_0000396"
        request.qualifiers.extend([ontologyterm])
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertGreater(len(response.phenotypes), 0)

    def testPhenotypeSearchMultipleQualifiers(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        ontologyterm = protocol.OntologyTerm()
        ontologyterm.term_id = "http://purl.obolibrary.org/obo/PATO_0000396"
        ontologyterm2 = protocol.OntologyTerm()
        ontologyterm2.term_id = "http://purl.obolibrary.org/obo/PATO_0000460"
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
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
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
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        request.description = ".*sensitivity.*"
        postUrl = '/phenotypes/search'
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchPhenotypesResponse)
        self.assertEquals(7, len(response.phenotypes))

    def testPhenotypesSearchMultipleTerms(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        request.description = "Melanoma, NOS with response to therapy"
        request.age_of_onset.term_id = \
            "http://purl.obolibrary.org/obo/HP_0003581"
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
        # simulate user interacting with sequenceAnnotations
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        featureId = \
            "http://ohsu.edu/cgd/27d2169c"
        obfuscated = self.getObfuscatedFeatureCompoundId(
            datasetName, featureSet.name, featureId)
        # use the feature to look up associations
        request.feature_ids.extend([obfuscated])
        response = self.sendSearchRequest(
            '/featurephenotypeassociations/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations[0].feature_ids))

    def testGenotypePhenotypeSearchEvidence(self):
        """
        Search for associations given an evidence
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        eq = protocol.EvidenceQuery()
        eq.description = "imatinib"
        request.evidence.extend([eq])
        response = self.sendSearchRequest(
            '/featurephenotypeassociations/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations[0].feature_ids))

    def testGenotypePhenotypeSearchPhenotype(self):
        """
        Search for associations given a phenotype
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        request.phenotype_ids.extend(["http://ohsu.edu/cgd/25abbb09"])
        response = self.sendSearchRequest(
            '/featurephenotypeassociations/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations[0].feature_ids))

    def testNoFind(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        request.feature_ids.extend(["FOOBAR"])
        response = self.sendSearchRequest(
            '/featurephenotypeassociations/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(0, len(response.associations))

    def testGenotypePhenotypeSearchEnsureEvidence(self):
        """
        Ensure evidence level is serialized in responses
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        featureId = \
            "http://ohsu.edu/cgd/27d2169c"
        obfuscated = self.getObfuscatedFeatureCompoundId(datasetName,
                                                         featureSet.name,
                                                         featureId)
        request.feature_ids.extend([obfuscated])
        response = self.sendSearchRequest(
            '/featurephenotypeassociations/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(1, len(response.associations[0].evidence))
        evidence = response.associations[0].evidence[0]
        self.assertEqual('decreased_sensitivity', evidence.description)

    def testGenotypePhenotypeSearchEnsureEnvironment(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = \
            self.getPhenotypeAssociationSetId()
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        featureId = \
            "http://ohsu.edu/cgd/27d2169c"
        obfuscated = self.getObfuscatedFeatureCompoundId(
            datasetName, featureSet.name, featureId)
        request.feature_ids.extend([obfuscated])
        eq = protocol.EvidenceQuery()
        eq.description = "imatinib"
        request.evidence.extend([eq])
        response = self.sendSearchRequest(
            '/featurephenotypeassociations/search',
            request,
            protocol.SearchGenotypePhenotypeResponse)
        self.assertEqual(
            1, len(response.associations[0].environmental_contexts))
        environmentalContext = response.associations[0] \
                                       .environmental_contexts[0]
        self.assertEqual('imatinib', environmentalContext.description)

    def _createPagingRequest(self):
        request = protocol.SearchFeaturesRequest()
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        request.feature_set_id = featureSet.id
        request.name = "KIT *wild"
        return request

    def testGenotypeSearchFeaturePagingOne(self):
        """
        If page size is set to 1 only one association should be returned
        """
        request = self._createPagingRequest()
        request.page_size = 1
        postUrl = "features/search"
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchFeaturesResponse)
        self.assertEqual(1, len(response.features))
        self.assertIsNotNone(response.next_page_token)

    def testGenotypeSearchFeaturePagingMore(self):
        """
        If page size is not set to more than one association should be returned
        """
        request = self._createPagingRequest()
        postUrl = "features/search"
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchFeaturesResponse)
        self.assertGreater(len(response.features), 1)
        self.assertEqual(response.next_page_token, '')

    def testGenotypeSearchFeaturePagingAll(self):
        """
        Loop through all pages
        """
        request = self._createPagingRequest()
        request.page_size = 1
        feature_set_id = request.feature_set_id
        postUrl = "features/search"
        response = self.sendSearchRequest(
            postUrl,
            request,
            protocol.SearchFeaturesResponse)

        self.assertEqual(1, len(response.features))
        self.assertIsNotNone(response.next_page_token)
        pageCount = 1

        while response.next_page_token:
            previous_id = response.features[0].id
            request = protocol.SearchFeaturesRequest()
            request.feature_set_id = feature_set_id
            request.page_size = 1
            request.page_token = response.next_page_token
            request.name = "KIT *wild"
            response = self.sendSearchRequest(
                postUrl,
                request,
                protocol.SearchFeaturesResponse)
            self.assertEqual(1, len(response.features))
            self.assertNotEqual(previous_id, response.features[0].id)
            pageCount += 1
        self.assertEqual(3, pageCount)

    def testGenotypesSearchByNameError(self):
        """
        Search for feature by name with a malformed regular expression.
        """
        # setup phenotype query
        request = protocol.SearchFeaturesRequest()
        datasetName, featureSet = self.getCGDDataSetFeatureSet()
        request.feature_set_id = featureSet.id
        request.name = "*"  # invalid regular expression

        postUrl = "features/search"
        response = self.sendJsonPostRequest(postUrl, protocol.toJson(request))
        self.assertEqual(400, response.status_code)
