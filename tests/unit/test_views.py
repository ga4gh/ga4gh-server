"""
Unit tests for the frontend code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import logging

import tests.paths as paths

import ga4gh.server.datamodel as datamodel
import ga4gh.server.frontend as frontend
import ga4gh.schemas.protocol as protocol
import ga4gh.server.exceptions as exceptions


class TestFrontend(unittest.TestCase):
    """
    Tests the basic routing and HTTP handling for the Flask app.
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
            "LANDING_MESSAGE_HTML": paths.landingMessageHtml
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
        cls.dataRepo = cls.backend.getDataRepository()
        cls.referenceSet = cls.dataRepo.getReferenceSets()[0]
        cls.referenceSetId = cls.referenceSet.getId()
        cls.reference = cls.referenceSet.getReferences()[0]
        cls.referenceId = cls.reference.getId()
        cls.dataset = cls.backend.getDataRepository().getDatasets()[0]
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
        cls.phenotypeAssociationSet = \
            cls.dataset.getPhenotypeAssociationSets()[0]
        cls.phenotypeAssociationSetId = cls.phenotypeAssociationSet.getId()
        cls.association = cls.phenotypeAssociationSet.getAssociations()[0]
        cls.phenotype = cls.association.phenotype
        cls.phenotypeId = cls.phenotype.id
        cls.featureSets = cls.dataset.getFeatureSets()
        cls.genotypePhenotype = cls.phenotypeAssociationSet.getAssociations(
            request=None, featureSets=cls.featureSets)[0]
        cls.genotypePhenotypeId = cls.genotypePhenotype.id
        cls.rnaQuantificationSet = cls.dataset.getRnaQuantificationSets()[0]
        cls.rnaQuantificationSetId = cls.rnaQuantificationSet.getId()
        cls.rnaQuantification = cls.rnaQuantificationSet.getRnaQuantifications(
            )[0]
        cls.rnaQuantificationId = cls.rnaQuantification.getId()
        cls.expressionLevel = cls.rnaQuantification.getExpressionLevels(
            1, 2)[0]
        cls.expressionLevelId = cls.expressionLevel.getId()

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
        variantSets = protocol.fromJson(
            response.data, protocol.SearchVariantSetsResponse).variant_sets
        request = protocol.SearchVariantsRequest()
        request.variant_set_id = variantSets[0].id
        request.reference_name = "1"
        request.start = 0
        request.end = 1
        return self.sendPostRequest('/variants/search', request)

    def sendVariantSetsSearch(self):
        request = protocol.SearchVariantSetsRequest()
        request.dataset_id = self.datasetId
        return self.sendPostRequest('/variantsets/search', request)

    def sendCallSetsSearch(self):
        response = self.sendVariantSetsSearch()
        variantSets = protocol.fromJson(
            response.data, protocol.SearchVariantSetsResponse).variant_sets
        request = protocol.SearchCallSetsRequest()
        request.variant_set_id = variantSets[0].id
        return self.sendPostRequest('/callsets/search', request)

    def sendReadsSearch(self, readGroupIds=None, referenceId=""):
        request = protocol.SearchReadsRequest()
        request.read_group_ids.extend(readGroupIds)
        request.reference_id = referenceId
        return self.sendPostRequest('/reads/search', request)

    def sendDatasetsSearch(self):
        request = protocol.SearchDatasetsRequest()
        return self.sendPostRequest('/datasets/search', request)

    def sendPhenotypesSearch(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = self.phenotypeAssociationSetId
        return self.sendPostRequest('/phenotypes/search', request)

    def sendGenotypePhenotypesSearch(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = self.phenotypeAssociationSetId
        return self.sendPostRequest(
            '/featurephenotypeassociations/search', request)

    def sendPhenotypeAssociationSetsSearch(self):
        request = protocol.SearchPhenotypeAssociationSetsRequest()
        request.dataset_id = self.datasetId
        return self.sendPostRequest(
            '/phenotypeassociationsets/search', request)

    def sendRnaQuantificationSetsSearch(self):
        request = protocol.SearchRnaQuantificationSetsRequest()
        request.dataset_id = self.datasetId
        return self.sendPostRequest(
            '/rnaquantificationsets/search', request)

    def sendRnaQuantificationsSearch(self):
        request = protocol.SearchRnaQuantificationsRequest()
        request.rna_quantification_set_id = self.rnaQuantificationSetId
        return self.sendPostRequest('/rnaquantifications/search', request)

    def sendExpressionLevelsSearch(self):
        request = protocol.SearchExpressionLevelsRequest()
        request.rna_quantification_id = self.rnaQuantificationId
        return self.sendPostRequest('/expressionlevels/search', request)

    def sendReferencesSearch(self):
        path = "/references/search"
        request = protocol.SearchReferencesRequest()
        response = self.sendPostRequest(path, request)
        return response

    def sendGetObject(self, id_, defaultId, path):
        if id_ is None:
            id_ = defaultId
        getPath = path.format(id_)
        response = self.sendGetRequest(getPath)
        return response

    def sendGetVariant(self, id_=None):
        return self.sendGetObject(id_, self.variantId, "/variants/{}")

    def sendGetVariantSet(self, id_=None):
        return self.sendGetObject(id_, self.variantSetId, "/variantsets/{}")

    def sendGetDataset(self, id_=None):
        return self.sendGetObject(id_, self.datasetId, "/datasets/{}")

    def sendGetReadGroup(self, id_=None):
        return self.sendGetObject(id_, self.readGroupId, "/readgroups/{}")

    def sendGetReference(self, id_=None):
        return self.sendGetObject(id_, self.referenceId, "/references/{}")

    def sendGetReadGroupSet(self, id_=None):
        return self.sendGetObject(
            id_, self.readGroupSetId, "/readgroupsets/{}")

    def sendGetCallSet(self, id_=None):
        return self.sendGetObject(id_, self.callSetId, "/callsets/{}")

    def sendGetReferenceSet(self, id_=None):
        return self.sendGetObject(
            id_, self.referenceSetId, "/referencesets/{}")

    def sendGetRnaQuantificationSet(self, id_=None):
        return self.sendGetObject(
            id_, self.rnaQuantificationSetId, "/rnaquantificationsets/{}")

    def sendGetRnaQuantification(self, id_=None):
        return self.sendGetObject(
            id_, self.rnaQuantificationId, "/rnaquantifications/{}")

    def sendGetExpressionLevel(self, id_=None):
        return self.sendGetObject(
            id_, self.expressionLevelId, "/expressionlevels/{}")

    def sendListRequest(self, path, request):
        headers = {
            'Origin': self.exampleUrl,
        }
        data = protocol.toJsonDict(request)
        response = self.app.post(path, data=data, headers=headers)
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
            protocol.fromJson(
                response.get_data(), protocol.GAException)
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
        assertHeaders(self.sendPhenotypesSearch())
        assertHeaders(self.sendGenotypePhenotypesSearch())
        assertHeaders(self.sendPhenotypeAssociationSetsSearch())
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
        if not getDefined:
            getResponse = self.app.get(path)
            protocol.fromJson(
                getResponse.get_data(), protocol.GAException)
            self.assertEqual(405, getResponse.status_code)

        # Malformed requests should return 400
        for badJson in ["JSON", "<xml/>", "{]"]:
            badResponse = self.app.post(
                path, data=badJson,
                headers={'Content-type': 'application/json'})
            self.assertEqual(400, badResponse.status_code)

        # OPTIONS should return success
        self.assertEqual(200, self.app.options(path).status_code)

    def testRouteReferences(self):
        referenceId = self.referenceId
        path = '/references/{}'
        path = path.format(referenceId)
        self.assertEqual(200, self.app.get(path).status_code)
        path = '/listreferencebases'
        self.assertEqual(404, self.app.post(path).status_code)
        referenceSetId = self.referenceSetId
        path = '/referencesets/{}'
        path = path.format(referenceSetId)
        self.assertEqual(200, self.app.get(path).status_code)
        path = 'references/{}'
        self.verifySearchRouting('/referencesets/search', True)
        self.verifySearchRouting('/references/search', True)

    def testRouteCallSets(self):
        path = '/callsets/search'
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
        responseData = protocol.fromJson(
            response.data, protocol.SearchVariantsResponse)
        self.assertEqual(len(responseData.variants), 1)

    def testVariantSetsSearch(self):
        response = self.sendVariantSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(
            response.data, protocol.SearchVariantSetsResponse)
        self.assertEqual(len(responseData.variant_sets), 1)

    def testGetDataset(self):
        # Test OK: ID found
        response = self.sendDatasetsSearch()
        responseData = protocol.fromJson(
            response.data, protocol.SearchDatasetsResponse)
        datasetId = responseData.datasets[0].id
        response = self.sendGetDataset(datasetId)
        self.assertEqual(200, response.status_code)

        # Test Error: 404, ID not found
        invalidId = datamodel.DatasetCompoundId.getInvalidIdString()
        obfuscated = datamodel.CompoundId.obfuscate(invalidId)
        compoundId = datamodel.DatasetCompoundId.parse(obfuscated)
        response = self.sendGetDataset(str(compoundId))
        self.assertEqual(404, response.status_code)

    def testGetVariantSet(self):
        response = self.sendVariantSetsSearch()
        responseData = protocol.fromJson(
            response.data, protocol.SearchVariantSetsResponse)
        variantSetId = responseData.variant_sets[0].id
        response = self.sendGetVariantSet(variantSetId)
        self.assertEqual(200, response.status_code)
        invalidId = datamodel.VariantSetCompoundId.getInvalidIdString()
        obfuscated = datamodel.CompoundId.obfuscate(invalidId)
        compoundId = datamodel.VariantSetCompoundId.parse(obfuscated)
        response = self.sendGetVariantSet(str(compoundId))
        self.assertEqual(404, response.status_code)

    def testCallSetsSearch(self):
        response = self.sendCallSetsSearch()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(
            response.data, protocol.SearchCallSetsResponse)
        self.assertEqual(len(responseData.call_sets), 1)

    def testReadsSearch(self):
        response = self.sendReadsSearch(readGroupIds=[self.readGroupId],
                                        referenceId=self.referenceId)
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(
            response.data, protocol.SearchReadsResponse)
        self.assertEqual(len(responseData.alignments), 2)
        self.assertEqual(
            responseData.alignments[0].id,
            self.readAlignmentId)

    def testPhenotypeAssociationSetsSearch(self):
        response = self.sendPhenotypeAssociationSetsSearch()
        responseData = protocol.fromJson(
            response.data, protocol.SearchPhenotypeAssociationSetsResponse)
        pasets = list(responseData.phenotype_association_sets)
        foundPASet = False
        for paset in pasets:
            if self.phenotypeAssociationSetId == paset.id:
                foundPASet = True
        self.assertTrue(foundPASet)

    def testNoAuthentication(self):
        path = '/oauth2callback'
        self.assertEqual(501, self.app.get(path).status_code)

    def testSearchUnmappedReads(self):
        response = self.sendReadsSearch(
            readGroupIds=[self.readGroupId], referenceId="")
        self.assertEqual(501, response.status_code)

    def testSearchReadsMultipleReadGroupSetsSetMismatch(self):
        response = self.sendReadsSearch(
            readGroupIds=[self.readGroupId, "42"],
            referenceId=self.referenceId)
        self.assertEqual(400, response.status_code)

    def getObjectTest(self, getMethod, protocolClass, objectId):
        response = getMethod()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(response.data, protocolClass)
        self.assertEqual(responseData.id, objectId)

    def testGetReadGroupSet(self):
        self.getObjectTest(
            self.sendGetReadGroupSet,
            protocol.ReadGroupSet,
            self.readGroupSetId)

    def testGetReadGroup(self):
        self.getObjectTest(
            self.sendGetReadGroup,
            protocol.ReadGroup,
            self.readGroupId)

    def testGetCallSet(self):
        self.getObjectTest(
            self.sendGetCallSet,
            protocol.CallSet,
            self.callSetId)

    def testGetVariant(self):
        self.getObjectTest(
            self.sendGetVariant,
            protocol.Variant,
            self.variantId)

    def testGetExpressionLevel(self):
        self.getObjectTest(
            self.sendGetExpressionLevel,
            protocol.ExpressionLevel,
            self.expressionLevelId)

    def testGetRnaQuantification(self):
        self.getObjectTest(
            self.sendGetRnaQuantification,
            protocol.RnaQuantification,
            self.rnaQuantificationId)

    def testGetRnaQuantificationSet(self):
        self.getObjectTest(
            self.sendGetRnaQuantificationSet,
            protocol.RnaQuantificationSet,
            self.rnaQuantificationSetId)

    def searchObjectTest(
            self, responseMethod, responseClass, attributeName, objectId):
        response = responseMethod()
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(response.data, responseClass)
        responseList = getattr(responseData, attributeName)
        objectList = list(responseList)
        self.assertEqual(objectId, objectList[0].id)

    def testDatasetsSearch(self):
        self.searchObjectTest(
            self.sendDatasetsSearch,
            protocol.SearchDatasetsResponse,
            "datasets",
            self.datasetId)

    def testPhenotypesSearch(self):
        self.searchObjectTest(
            self.sendPhenotypesSearch,
            protocol.SearchPhenotypesResponse,
            "phenotypes",
            self.phenotypeId)

    def testGenotypePhenotypesSearch(self):
        self.searchObjectTest(
            self.sendGenotypePhenotypesSearch,
            protocol.SearchGenotypePhenotypeResponse,
            "associations",
            self.genotypePhenotypeId)

    def testExpressionLevelsSearch(self):
        self.searchObjectTest(
            self.sendExpressionLevelsSearch,
            protocol.SearchExpressionLevelsResponse,
            "expression_levels",
            self.expressionLevelId)

    def testRnaQuantificationsSearch(self):
        self.searchObjectTest(
            self.sendRnaQuantificationsSearch,
            protocol.SearchRnaQuantificationsResponse,
            "rna_quantifications",
            self.rnaQuantificationId)

    def testRnaQuantificationSetsSearch(self):
        self.searchObjectTest(
            self.sendRnaQuantificationSetsSearch,
            protocol.SearchRnaQuantificationSetsResponse,
            "rna_quantification_sets",
            self.rnaQuantificationSetId)

    def testNoCallback(self):
        response = self.sendGetRequest("callback")
        self.assertEqual(
            response.status_code,
            404, "Ensure that when Auth0 is turned off the callback"
                 "URL returns a 404 but got {}".format(response.status_code))

    def testNoLogin(self):
        response = self.sendGetRequest("login")
        self.assertEqual(
            response.status_code,
            404, "Ensure that when Auth0 is turned off the callback"
                 "URL returns a 404 but got {}".format(response.status_code))

    def testSimplePost(self):
        path = "/datasets/search"
        response = protocol.fromJson(self.app.post(
            path, headers={}).get_data(), protocol.SearchDatasetsResponse)
        self.assertIsNotNone(
            response.datasets,
            "When an empty JSON document "
            "without a mimetype is sent we can still"
            "get datasets.")

    def testHandleHttpPost(self):

        class Mock(object):
            pass
        request = Mock()
        request.mimetype = "garbage"
        # A bad mimetype should throw an exception
        with self.assertRaises(exceptions.UnsupportedMediaTypeException):
            response = frontend.handleHttpPost(request, lambda x: x)

        # An empty mimetype should work OK
        request = Mock()
        request.mimetype = None
        request.get_data = lambda: "data"
        response = frontend.handleHttpPost(request, lambda x: x)
        self.assertEquals(response.get_data(), "data")
