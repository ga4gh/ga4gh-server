"""
Unit tests for the frontend code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import logging

import ga4gh.datamodel as datamodel
import ga4gh.frontend as frontend
import ga4gh.protocol as protocol


class TestFrontend(unittest.TestCase):
    """
    Tests the basic routing and HTTP handling for the Flask app.
    """
    exampleUrl = 'www.example.com'

    @classmethod
    def setUpClass(cls):
        config = {
            "DATA_SOURCE": "__SIMULATED__",
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
        cls.referenceSet = cls.backend.getReferenceSets()[0]
        cls.referenceSetId = cls.referenceSet.getId()
        cls.reference = cls.referenceSet.getReferences()[0]
        cls.referenceId = cls.reference.getId()
        cls.dataset = cls.backend.getDatasets()[0]
        cls.datasetId = cls.dataset.getId()
        cls.variantSet = cls.dataset.getVariantSets()[0]
        cls.variantSetId = cls.variantSet.getId()
        gaVariant = cls.variantSet.getVariants("1", 0, 2**32).next()
        cls.variantId = gaVariant.id
        cls.callSet = cls.variantSet.getCallSets()[0]
        cls.callSetId = cls.callSet.getId()
        cls.readGroupSet = cls.dataset.getReadGroupSets()[0]
        cls.readGroupSetId = cls.readGroupSet.getId()
        cls.readGroup = cls.readGroupSet.getReadGroups()[0]
        cls.readGroupId = cls.readGroup.getId()
        cls.readAlignment = cls.readGroup.getReadAlignments().next()
        cls.readAlignmentId = cls.readAlignment.id

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

    def sendGetRequest(self, path):
        """
        Sends a get request to the specified URL and returns the response.
        """
        headers = {
            'Origin': self.exampleUrl,
        }
        return self.app.get(path, headers=headers)

    def sendVariantsSearch(self):
        response = self.sendVariantSetsSearch()
        variantSets = protocol.SearchVariantSetsResponse().fromJsonString(
            response.data).variantSets
        request = protocol.SearchVariantsRequest()
        request.variantSetId = variantSets[0].id
        request.referenceName = "1"
        request.start = 0
        request.end = 1
        return self.sendPostRequest('/variants/search', request)

    def sendVariantSetsSearch(self):
        request = protocol.SearchVariantSetsRequest()
        request.datasetId = self.datasetId
        return self.sendPostRequest('/variantsets/search', request)

    def sendCallSetsSearch(self):
        response = self.sendVariantSetsSearch()
        variantSets = protocol.SearchVariantSetsResponse().fromJsonString(
            response.data).variantSets
        request = protocol.SearchCallSetsRequest()
        request.variantSetId = variantSets[0].id
        return self.sendPostRequest('/callsets/search', request)

    def sendReadsSearch(self, readGroupIds=None):
        if readGroupIds is None:
            readGroupIds = [self.readGroupId]
        request = protocol.SearchReadsRequest()
        request.readGroupIds = readGroupIds
        request.referenceId = self.referenceId
        return self.sendPostRequest('/reads/search', request)

    def sendDatasetsSearch(self):
        request = protocol.SearchDatasetsRequest()
        return self.sendPostRequest('/datasets/search', request)

    def sendReferencesSearch(self):
        path = "/references/search"
        request = protocol.SearchReferencesRequest()
        response = self.sendPostRequest(path, request)
        return response

    def sendGenotypePhenotypeSearch(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        response = self.sendPostRequest('/genotypephenotype/search', request)
        return response

    def sendGetVariant(self, id_=None):
        if id_ is None:
            id_ = self.variantId
        path = "/variants/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetVariantSet(self, id_=None):
        if id_ is None:
            id_ = self.variantSetId
        path = "/variantsets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetDataset(self, id_=None):
        if id_ is None:
            id_ = self.datasetId
        path = "/datasets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetReadGroup(self, id_=None):
        if id_ is None:
            id_ = self.readGroupId
        path = "/readgroups/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetReference(self, id_=None):
        if id_ is None:
            id_ = self.referenceId
        path = "/references/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetReadGroupSet(self, id_=None):
        if id_ is None:
            id_ = self.readGroupSetId
        path = "/readgroupsets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetCallset(self, id_=None):
        if id_ is None:
            id_ = self.callSetId
        path = "/callsets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendGetReferenceSet(self, id_=None):
        if id_ is None:
            id_ = self.referenceSetId
        path = "/referencesets/{}".format(id_)
        response = self.sendGetRequest(path)
        return response

    def sendListRequest(self, path, request):
        headers = {
            'Origin': self.exampleUrl,
        }
        data = request.toJsonDict()
        response = self.app.get(path, data=data, headers=headers)
        return response

    def sendReferenceBasesList(self, id_=None):
        if id_ is None:
            id_ = self.referenceId
        path = "/references/{}/bases".format(id_)
        request = protocol.ListReferenceBasesRequest()
        response = self.sendListRequest(path, request)
        return response

    def test404sReturnJson(self):
        paths = [
            '/doesNotExist',
            '/reads/sea',
            '/variantsets/id/doesNotExist',
        ]
        for path in paths:
            response = self.app.get(path)
            protocol.GAException.fromJsonString(response.get_data())
            self.assertEqual(404, response.status_code)

    def testCors(self):
        def assertHeaders(response):
            self.assertEqual(self.exampleUrl,
                             response.headers['Access-Control-Allow-Origin'])
            self.assertTrue('Content-Type' in response.headers)
        # Post-based search methods
        assertHeaders(self.sendVariantsSearch())
        assertHeaders(self.sendVariantSetsSearch())
        assertHeaders(self.sendReadsSearch())
        assertHeaders(self.sendReferencesSearch())
        assertHeaders(self.sendReferenceBasesList())
        assertHeaders(self.sendDatasetsSearch())
        assertHeaders(self.sendGenotypePhenotypeSearch())
        # Get-based accessor methods
        assertHeaders(self.sendGetVariantSet())
        assertHeaders(self.sendGetReference())
        assertHeaders(self.sendGetReferenceSet())
        assertHeaders(self.sendGetReadGroupSet())
        assertHeaders(self.sendGetReadGroup())
        assertHeaders(self.sendGetVariant())
        assertHeaders(self.sendGetDataset())
        # TODO: Test other methods as they are implemented

    def verifySearchRouting(self, path, getDefined=False):
        """
        Verifies that the specified path has the correct routing for a
        search command. If getDefined is False we check to see if it
        returns the correct status code.
        """
        response = self.app.post(path)
        protocol.GAException.fromJsonString(response.get_data())
        self.assertEqual(415, response.status_code)
        if not getDefined:
            getResponse = self.app.get(path)
            protocol.GAException.fromJsonString(getResponse.get_data())
            self.assertEqual(405, getResponse.status_code)

        # Malformed requests should return 400
        for badJson in ["", None, "JSON", "<xml/>", "{]"]:
            badResponse = self.app.post(
                path, data=badJson,
                headers={'Content-type': 'application/json'})
            self.assertEqual(400, badResponse.status_code)

        # OPTIONS should return success
        self.assertEqual(200, self.app.options(path).status_code)

    def testRouteReferences(self):
        referenceId = self.referenceId
        paths = ['/references/{}', '/references/{}/bases']
        for path in paths:
            path = path.format(referenceId)
            self.assertEqual(200, self.app.get(path).status_code)
        referenceSetId = self.referenceSetId
        paths = ['/referencesets/{}']
        for path in paths:
            path = path.format(referenceSetId)
            self.assertEqual(200, self.app.get(path).status_code)
        self.verifySearchRouting('/referencesets/search', True)
        self.verifySearchRouting('/references/search', True)

    def testRouteCallsets(self):
        path = '/callsets/search'
        self.assertEqual(415, self.app.post(path).status_code)
        self.assertEqual(200, self.app.options(path).status_code)
        self.assertEqual(405, self.app.get(path).status_code)

    def testRouteReads(self):
        paths = ['/reads/search', '/readgroupsets/search']
        for path in paths:
            self.verifySearchRouting(path)

    def testRouteVariants(self):
        self.verifySearchRouting('/variantsets/search', True)
        self.verifySearchRouting('/variants/search', False)

    def testRouteIndex(self):
        path = "/"
        response = self.app.get(path)
        self.assertEqual(200, response.status_code)
        self.assertEqual("text/html", response.mimetype)
        self.assertGreater(len(response.data), 0)

    def testVariantsSearch(self):
        response = self.sendVariantsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchVariantsResponse.fromJsonString(
            response.data)
        self.assertEqual(len(responseData.variants), 1)

    def testVariantSetsSearch(self):
        response = self.sendVariantSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchVariantSetsResponse.fromJsonString(
            response.data)
        self.assertEqual(len(responseData.variantSets), 1)

    def testGetDataset(self):
        # Test OK: ID found
        response = self.sendDatasetsSearch()
        responseData = protocol.SearchDatasetsResponse.fromJsonString(
            response.data)
        datasetId = responseData.datasets[0].id
        response = self.sendGetDataset(datasetId)
        self.assertEqual(200, response.status_code)

        # Test Error: 404, ID not found
        obfuscated = datamodel.CompoundId.obfuscate("notValid")
        compoundId = datamodel.DatasetCompoundId.parse(obfuscated)
        response = self.sendGetDataset(str(compoundId))
        self.assertEqual(404, response.status_code)

    def testGetVariantSet(self):
        response = self.sendVariantSetsSearch()
        responseData = protocol.SearchVariantSetsResponse.fromJsonString(
            response.data)
        variantSetId = responseData.variantSets[0].id
        response = self.sendGetVariantSet(variantSetId)
        self.assertEqual(200, response.status_code)
        obfuscated = datamodel.CompoundId.obfuscate("notValid:notValid")
        compoundId = datamodel.VariantSetCompoundId.parse(obfuscated)
        response = self.sendGetVariantSet(str(compoundId))
        self.assertEqual(404, response.status_code)

    def testGetReadGroupSet(self):
        response = self.sendGetReadGroupSet()
        self.assertEqual(200, response.status_code)
        responseData = protocol.ReadGroupSet.fromJsonString(
            response.data)
        self.assertEqual(
            responseData.id, self.readGroupSetId)

    def testGetReadGroup(self):
        response = self.sendGetReadGroup()
        self.assertEqual(200, response.status_code)
        responseData = protocol.ReadGroup.fromJsonString(
            response.data)
        self.assertEqual(
            responseData.id, self.readGroupId)

    def testGetCallset(self):
        response = self.sendGetCallset()
        self.assertEqual(200, response.status_code)
        responseData = protocol.CallSet.fromJsonString(
            response.data)
        self.assertEqual(
            responseData.id, self.callSetId)

    def testGetVariant(self):
        response = self.sendGetVariant()
        self.assertEqual(200, response.status_code)

    def testCallSetsSearch(self):
        response = self.sendCallSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchCallSetsResponse.fromJsonString(
            response.data)
        self.assertEqual(len(responseData.callSets), 1)

    def testReadsSearch(self):
        response = self.sendReadsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.SearchReadsResponse.fromJsonString(
            response.data)
        self.assertEqual(len(responseData.alignments), 2)
        self.assertEqual(
            responseData.alignments[0].id,
            self.readAlignmentId)

    def testDatasetsSearch(self):
        response = self.sendDatasetsSearch()
        responseData = protocol.SearchDatasetsResponse.fromJsonString(
            response.data)
        datasets = list(responseData.datasets)
        self.assertEqual(self.datasetId, datasets[0].id)

    def testGenotypePhenotypeSearchBadRequest(self):
        """
        If nothing is provided to search on,
        the system should return `bad request`
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(501, response.status_code)

    def testNoAuthentication(self):
        path = '/oauth2callback'
        self.assertEqual(501, self.app.get(path).status_code)

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
        id.database = "http://www.monarchinitiative.org/_"
        id.identifier = "CGD:d8c2d551UniProtKB:P10721#P10721-1Region"
        id.version = "*"
        request.feature.ids = [id]
        response = self.sendPostRequest('/genotypephenotype/search', request)
        self.assertEqual(200, response.status_code)
        response = protocol.SearchGenotypePhenotypeResponse().fromJsonString(
            response.data)
        self.assertEqual(1, len(response.associations[0].features))

        request.phenotype = protocol.ExternalIdentifierQuery()
        id = protocol.ExternalIdentifier()
        id.database = "http://purl.obolibrary.org/obo/OMIM_"
        id.identifier = "606764"
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
        id.identifier = "DB00619"
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
            self.assertIsNotNone(response.nextPageToken)
        # from IPython.core.debugger import Pdb ;        Pdb().set_trace()
