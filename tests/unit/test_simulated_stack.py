"""
End-to-end tests for the simulator configuration. Sets up a server with
the backend, sends some basic queries to that server and verifies results
are as expected.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import logging
import os
import tempfile
import unittest

import ga4gh.frontend as frontend
import ga4gh.protocol as protocol
import tests.utils as utils


class TestSimulatedStack(unittest.TestCase):
    """
    Tests the full stack for the Simulated backend by using the Flask
    testing client.
    """
    @classmethod
    def setUpClass(cls):
        # silence usually unhelpful CORS log
        logging.getLogger('ga4gh.frontend.cors').setLevel(logging.CRITICAL)
        fd, cls.db_file = tempfile.mkstemp(prefix="ga4gh_test_simulated_stack")
        os.close(fd)
        db_url = "sqlite:///" + cls.db_file
        registry_db = utils.create_simulated_registry_db(
            db_url=db_url, random_seed=1111, num_calls=5, num_variant_sets=4,
            num_variant_annotation_sets=3, num_reference_sets=3,
            num_references_per_reference_set=4, num_read_group_sets=3,
            num_read_groups_per_read_group_set=2, num_feature_sets=5)
        registry_db.close()
        config = {
            "DATA_SOURCE": db_url,
            # "DEBUG": True
        }
        frontend.reset()
        frontend.configure(
            baseConfig="TestConfig", extraConfig=config)
        cls.app = frontend.app.test_client()

    @classmethod
    def tearDownClass(cls):
        cls.app = None
        os.unlink(cls.db_file)

    def setUp(self):
        self.backend = frontend.app.backend
        self.registry_db = self.backend.get_registry_db()
        self.dataset = self.registry_db.get_datasets()[0]
        self.readGroupSet = self.dataset.read_group_sets[0]
        self.readGroup = self.readGroupSet.read_groups[0]
        referenceSet = self.readGroupSet.reference_set
        self.reference = referenceSet.references[0]
        self.variantSet = self.dataset.variant_sets[0]
        self.variantAnnotationSet = self.variantSet.variant_annotation_sets[0]
        self.backend.setMaxResponseLength(10000)

    def getBadIds(self):
        """
        Returns a list of IDs that should not exist in the server and should
        raise a 404 error.
        """
        return ["1234:", "x"*100, ":", ":xx", "::", ":::", "::::"]

    def sendJsonPostRequest(self, path, data):
        """
        Sends a JSON request to the specified path with the specified data
        and returns the response.
        """
        return self.app.post(
            path, headers={'Content-type': 'application/json'},
            data=data)

    def sendSearchRequest(self, path, request, responseClass):
        """
        Sends the specified protocol request instance as JSON, and
        parses the result into an instance of the specified response.
        """
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        self.assertEqual(200, response.status_code)
        responseData = protocol.fromJson(response.data, responseClass)
        self.assertTrue(
            protocol.validate(
                protocol.toJson(responseData),
                type(responseData)))
        return responseData

    def sendObjectGetRequest(self, path, id_):
        """
        Sends a GET request to the specified path for an object with the
        specified ID and returns the response.
        """
        return self.app.get("{}/{}".format(path, id_))

    def sendGetObject(self, path, id_, responseClass):
        """
        Sends a get request and parses the value into an instance of the
        specified class.
        """
        response = self.sendObjectGetRequest(path, id_)
        self.assertEqual(200, response.status_code)
        obj = protocol.fromJson(response.data, responseClass)
        self.assertIsInstance(obj, responseClass)
        return obj

    def sendListReferenceBasesRequest(self, id_, request):
        """
        Sends a ListReferenceBasesRequest and parses the result into a
        ListReferenceBasesResponse.
        """
        path = '/references/{}/bases'.format(id_)
        response = self.app.get(
            path, query_string=protocol.toJsonDict(request))
        self.assertEqual(response.status_code, 200)
        obj = protocol.fromJson(
            response.data, protocol.ListReferenceBasesResponse)
        self.assertIsInstance(obj, protocol.ListReferenceBasesResponse)
        return obj

    def verifyVariantSetsEqual(self, gaVariantSet, variantSet):
        dataset = variantSet.dataset
        self.assertEqual(gaVariantSet.id, str(variantSet.id))
        self.assertEqual(gaVariantSet.dataset_id, str(dataset.id))
        self.assertEqual(gaVariantSet.name, variantSet.name)
        metadataList = sorted(
            variantSet.variant_set_metadata, key=lambda x: str(x.id))
        gaMetadataList = sorted(gaVariantSet.metadata, key=lambda x: x.id)
        self.assertEqual(len(metadataList), len(gaMetadataList))
        for metadata, gaMetadata in zip(metadataList, gaMetadataList):
            self.assertEqual(str(metadata.id), gaMetadata.id)
            self.assertEqual(metadata.key, gaMetadata.key)
            self.assertEqual(metadata.value, gaMetadata.value)
            self.assertEqual(metadata.type, gaMetadata.type)
            self.assertEqual(metadata.number, gaMetadata.number)
            self.assertEqual(metadata.description, gaMetadata.description)

    def verifyVariantAnnotationSetsEqual(self, ga_vas, vas):
        self.assertEqual(ga_vas.id, str(vas.id))
        self.assertEqual(ga_vas.name, str(vas.name))
        self.assertEqual(ga_vas.variantSet_id, str(vas.variant_set_id))
        ga_analysis = ga_vas.analysis
        analysis = vas.analysis
        self.assertEqual(ga_analysis.id, str(analysis.id))
        self.assertEqual(ga_analysis.name, analysis.name)
        self.assertEqual(ga_analysis.description, analysis.description)
        self.assertEqual(ga_analysis.software[0], analysis.software)
        self.assertEqual(
            ga_analysis.created,
            protocol.datetime_to_iso8601(analysis.creation_timestamp))
        self.assertEqual(
            ga_analysis.updated,
            protocol.datetime_to_iso8601(analysis.update_timestamp))

    def verifyCallSetsEqual(self, gaCallSet, callSet):
        self.assertEqual(gaCallSet.id, str(callSet.id))
        self.assertEqual(gaCallSet.name, callSet.name)
        self.assertEqual(
            gaCallSet.variant_set_ids,
            [str(vs.id) for vs in callSet.variant_sets])
        for key, value in gaCallSet.info.items():
            self.assertEqual(value[0], callSet.getInfo()[key])

    def verifyReadGroupSetsEqual(self, gaReadGroupSet, readGroupSet):
        dataset = readGroupSet.dataset
        self.assertEqual(gaReadGroupSet.id, str(readGroupSet.id))
        self.assertEqual(gaReadGroupSet.dataset_id, str(dataset.id))
        self.assertEqual(gaReadGroupSet.name, readGroupSet.name)
        self.assertEqual(
            len(gaReadGroupSet.read_groups), len(readGroupSet.read_groups))
        for gaReadGroup, readGroup in zip(
                gaReadGroupSet.read_groups, readGroupSet.read_groups):
            self.verifyReadGroupsEqual(gaReadGroup, readGroup)

    def verifyReadGroupsEqual(self, gaReadGroup, readGroup):
        self.assertEqual(gaReadGroup.id, str(readGroup.id))
        self.assertEqual(gaReadGroup.name, readGroup.name)
        self.assertEqual(gaReadGroup.description, readGroup.description)
        self.assertEqual(
            gaReadGroup.bio_sample_id, str(readGroup.bio_sample_id))
        self.assertEqual(
            gaReadGroup.dataset_id, str(readGroup.read_group_set.dataset.id))
        self.assertEqual(
            gaReadGroup.reference_set_id,
            str(readGroup.read_group_set.reference_set.id))
        self.assertEqual(
            gaReadGroup.predicted_insert_size, readGroup.predicted_insert_size)
        # TODO stats, programs, experiment.

    def verifyDatasetsEqual(self, gaDataset, dataset):
        self.assertEqual(gaDataset.id, str(dataset.id))
        self.assertEqual(gaDataset.name, dataset.name)
        self.assertEqual(gaDataset.description, dataset.description)

    def verifyReferenceSetsEqual(self, gaReferenceSet, referenceSet):
        self.assertEqual(gaReferenceSet.id, str(referenceSet.id))
        self.assertEqual(
            gaReferenceSet.md5checksum, referenceSet.md5checksum)
        self.assertEqual(
            gaReferenceSet.ncbi_taxon_id, referenceSet.ncbi_taxon_id)
        self.assertEqual(
            gaReferenceSet.assembly_id, referenceSet.assembly_id)
        self.assertEqual(
            gaReferenceSet.source_uri, referenceSet.source_uri)
        self.assertEqual(
            gaReferenceSet.source_accessions,
            [acc.name for acc in referenceSet.source_accessions])
        self.assertEqual(gaReferenceSet.is_derived, referenceSet.is_derived)
        self.assertEqual(gaReferenceSet.name, referenceSet.name)

    def verifyFeatureSetsEqual(self, ga_feature_set, feature_set):
        dataset = feature_set.dataset
        self.assertEqual(ga_feature_set.id, str(feature_set.id))
        self.assertEqual(ga_feature_set.dataset_id, str(dataset.id))
        self.assertEqual(ga_feature_set.name, feature_set.name)
        self.assertEqual(ga_feature_set.source_uri, feature_set.source_uri)
        # TODO compare info values.

    def verifyFeaturesEquivalent(self, f1, f2):
        # at least modulo featureId. They can obviously have different
        # start/end/etc parameters if randomly generated from search vs get.
        self.assertEqual(f1.id, f2.id)
        self.assertEqual(f1.parent_id, f2.parent_id)
        self.assertEqual(f1.feature_set_id, f2.feature_set_id)

    def verifyReferencesEqual(self, gaReference, reference):
        self.assertEqual(gaReference.id, str(reference.id))
        self.assertEqual(gaReference.name, reference.name)
        self.assertEqual(gaReference.length, reference.length)
        self.assertEqual(gaReference.md5checksum, reference.md5checksum)
        self.assertEqual(gaReference.ncbi_taxon_id, reference.ncbi_taxon_id)
        self.assertEqual(gaReference.source_uri, reference.source_uri)
        self.assertEqual(
            gaReference.source_accessions,
            [acc.name for acc in reference.source_accessions])
        self.assertEqual(gaReference.is_derived, reference.is_derived)
        self.assertEqual(
            gaReference.source_divergence, reference.source_divergence)

    def verifyIndividualsEqual(self, gaIndividual, individual):
        self.assertEqual(gaIndividual.id, str(individual.id))
        self.assertEqual(gaIndividual.name, str(individual.name))
        self.assertEqual(
            gaIndividual.dataset_id, str(individual.dataset_id))
        # TODO time stamps and other fields.

    def verifyBioSamplesEqual(self, gaBioSample, bio_sample):
        self.assertEqual(gaBioSample.id, str(bio_sample.id))
        self.assertEqual(gaBioSample.name, str(bio_sample.name))
        self.assertEqual(
            gaBioSample.dataset_id, str(bio_sample.dataset_id))
        self.assertEqual(
            gaBioSample.individual_id, str(bio_sample.individual_id))
        # TODO time stamps and other fields.

    def verifySearchMethod(
            self, request, path, responseClass, objects, objectVerifier):
        """
        Verifies that the specified search request operates correctly
        and returns all the speficied objects. The specified verifier
        function checks that all the returned objects are equivalent
        to their datamodel counterparts.
        """
        request.page_size = len(objects)
        self.assertGreater(request.page_size, 0)
        responseData = self.sendSearchRequest(path, request, responseClass)
        self.assertEqual(responseData.next_page_token, "")
        responseList = getattr(
            responseData, protocol.getValueListName(responseClass))
        self.assertEqual(len(objects), len(responseList))
        # Sort both lists by ID
        objects = sorted(objects, key=lambda x: str(x.id))
        responseList = sorted(responseList, key=lambda x: x.id)
        for gaObject, datamodelObject in zip(responseList, objects):
            objectVerifier(gaObject, datamodelObject)

    def verifySearchResultsEmpty(self, request, path, responseClass):
        """
        Verifies that we get a successful response with an empty list of
        results.
        """
        responseData = self.sendSearchRequest(path, request, responseClass)
        self.assertEqual("", responseData.next_page_token)
        responseList = getattr(
            responseData, protocol.getValueListName(responseClass))
        self.assertEqual(0, len(responseList))

    def assertObjectNotFound(self, response):
        """
        Checks that the specified response contains a search failure.
        """
        self.assertEqual(404, response.status_code)
        error = protocol.fromJson(response.data, protocol.GAException)
        self.assertTrue(protocol.validate(protocol.toJson(error), type(error)))
        self.assertGreater(error.error_code, 0)
        self.assertGreater(len(error.message), 0)

    def assertObjectNotSupported(self, response):
        """
        Checks that the specified response returns a not supported 501 status
        """
        self.assertEqual(501, response.status_code)
        error = protocol.fromJson(response.data, protocol.GAException)
        self.assertTrue(protocol.validate(protocol.toJson(error), type(error)))
        self.assertGreater(error.error_code, 0)
        self.assertGreater(len(error.message), 0)

    def verifySearchMethodFails(self, request, path):
        """
        Verify that the specified search request fails with a 404.
        """
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        self.assertObjectNotFound(response)

    def verifySearchMethodNotSupported(self, request, path):
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        self.assertObjectNotSupported(response)

    def verifyGetMethodFails(self, path, id_):
        """
        Verifies the specified GET request failes with a 404.
        """
        response = self.sendObjectGetRequest(path, id_)
        self.assertObjectNotFound(response)

    def testGetDataset(self):
        path = "/datasets"
        for dataset in self.registry_db.get_datasets():
            responseObject = self.sendGetObject(
                path, dataset.id, protocol.Dataset)
            self.verifyDatasetsEqual(responseObject, dataset)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testDatasetsSearch(self):
        request = protocol.SearchDatasetsRequest()
        datasets = self.registry_db.get_datasets()
        path = '/datasets/search'
        self.verifySearchMethod(
            request, path, protocol.SearchDatasetsResponse, datasets,
            self.verifyDatasetsEqual)

    def testVariantSetsSearch(self):
        path = '/variantsets/search'
        for dataset in self.registry_db.get_datasets():
            variantSets = list(dataset.variant_sets)
            request = protocol.SearchVariantSetsRequest()
            request.dataset_id = str(dataset.id)
            self.verifySearchMethod(
                request, path, protocol.SearchVariantSetsResponse, variantSets,
                self.verifyVariantSetsEqual)
        for badId in self.getBadIds():
            request = protocol.SearchVariantSetsRequest()
            request.dataset_id = badId
            self.verifySearchMethodFails(request, path)

    @unittest.skip("TODO fix CallSet info field.")
    def testCallSetsSearch(self):
        path = '/callsets/search'
        for dataset in self.registry_db.get_datasets():
            for variantSet in dataset.variant_sets:
                callSets = list(variantSet.call_sets)
                self.assertGreater(len(callSets), 0)
                request = protocol.SearchCallSetsRequest()
                request.variant_set_id = str(variantSet.id)
                self.verifySearchMethod(
                    request, path, protocol.SearchCallSetsResponse, callSets,
                    self.verifyCallSetsEqual)
                # Check if we can search for the callset with a good name.
                for callSet in callSets:
                    request = protocol.SearchCallSetsRequest()
                    request.variant_set_id = str(variantSet.id)
                    request.name = callSet.name
                    self.verifySearchMethod(
                        request, path, protocol.SearchCallSetsResponse,
                        [callSet], self.verifyCallSetsEqual)
                # Check if we can search for the callset with a bad name.
                for badId in self.getBadIds():
                    request = protocol.SearchCallSetsRequest()
                    request.variant_set_id = str(variantSet.id)
                    request.name = badId
                    self.verifySearchResultsEmpty(
                        request, path, protocol.SearchCallSetsResponse)
        # Check for searches within missing variantSets.
        for badId in self.getBadIds():
            request = protocol.SearchCallSetsRequest()
            request.variant_set_id = badId
            self.verifySearchMethodFails(request, path)

    def testReadGroupSetsSearch(self):
        path = '/readgroupsets/search'
        for dataset in self.registry_db.get_datasets():
            readGroupSets = list(dataset.read_group_sets)
            request = protocol.SearchReadGroupSetsRequest()
            request.dataset_id = str(dataset.id)
            self.verifySearchMethod(
                request, path, protocol.SearchReadGroupSetsResponse,
                readGroupSets, self.verifyReadGroupSetsEqual)
            # Check if we can search for the readGroupSet with a good name.
            for readGroupSet in readGroupSets:
                request = protocol.SearchReadGroupSetsRequest()
                request.dataset_id = str(dataset.id)
                request.name = readGroupSet.name
                self.verifySearchMethod(
                    request, path, protocol.SearchReadGroupSetsResponse,
                    [readGroupSet], self.verifyReadGroupSetsEqual)
            # Check if we can search for the readGroupSet with a bad name.
            for badId in self.getBadIds():
                request = protocol.SearchReadGroupSetsRequest()
                request.dataset_id = str(dataset.id)
                request.name = badId
                self.verifySearchResultsEmpty(
                    request, path, protocol.SearchReadGroupSetsResponse)
        for badId in self.getBadIds():
            request = protocol.SearchReadGroupSetsRequest()
            request.dataset_id = badId
            self.verifySearchMethodFails(request, path)

    def testReferenceSetsSearch(self):
        request = protocol.SearchReferenceSetsRequest()
        referenceSets = self.registry_db.get_reference_sets()
        path = '/referencesets/search'
        self.verifySearchMethod(
            request, path, protocol.SearchReferenceSetsResponse, referenceSets,
            self.verifyReferenceSetsEqual)

    def testReferencesSearch(self):
        path = '/references/search'
        for referenceSet in self.registry_db.get_reference_sets():
            references = referenceSet.references
            request = protocol.SearchReferencesRequest()
            request.reference_set_id = str(referenceSet.id)
            self.verifySearchMethod(
                request, path, protocol.SearchReferencesResponse, references,
                self.verifyReferencesEqual)
        for badId in self.getBadIds():
            request = protocol.SearchReferencesRequest()
            request.reference_set_id = badId
            self.verifySearchMethodFails(request, path)

    def verifyReferenceSearchFilters(
            self, objectList, hasAssemblyId, path, requestFactory,
            responseClass, objectVerifier):
        """
        Verifies the filtering functionality for the specified list of
        reference-like objects.
        """
        self.assertGreater(len(objectList), 2)
        for obj in objectList[1:]:
            request = requestFactory()
            # First, check the simple cases; 1 filter set, others null.
            request.md5checksum = obj.md5checksum
            self.verifySearchMethod(
                request, path, responseClass, [obj], objectVerifier)
            request.md5checksum = ""
            request.accession = obj.source_accessions[0].name
            self.verifySearchMethod(
                request, path, responseClass, [obj], objectVerifier)
            request.accession = ""
            if hasAssemblyId:
                request.assembly_id = obj.assembly_id
                self.verifySearchMethod(
                    request, path, responseClass, [obj], objectVerifier)
                request.assembly_id = ""
            # Now check one good value and some bad values.
            request.md5checksum = obj.md5checksum
            badAccessions = [
                "no such accession", objectList[0].source_accessions[0].name]
            for accession in badAccessions:
                request.accession = accession
                self.verifySearchResultsEmpty(request, path, responseClass)
            request.accession = ""
            if hasAssemblyId:
                badAssemblyIds = [
                    "no such asssembly", objectList[0].assembly_id]
                for assemblyId in badAssemblyIds:
                    request.assembly_id = assemblyId
                    self.verifySearchResultsEmpty(request, path, responseClass)
                request.assembly_id = ""

    def testReferencesSearchFilters(self):
        path = '/references/search'
        for referenceSet in self.registry_db.get_reference_sets():

            def requestFactory():
                request = protocol.SearchReferencesRequest()
                request.reference_set_id = str(referenceSet.id)
                return request
            self.verifyReferenceSearchFilters(
                referenceSet.references, False, path, requestFactory,
                protocol.SearchReferencesResponse, self.verifyReferencesEqual)

    def testReferenceSetsSearchFilters(self):
        path = '/referencesets/search'

        def requestFactory():
            return protocol.SearchReferenceSetsRequest()
        self.verifyReferenceSearchFilters(
            self.registry_db.get_reference_sets(), True, path, requestFactory,
            protocol.SearchReferenceSetsResponse,
            self.verifyReferenceSetsEqual)

    def testGetVariantSet(self):
        path = "/variantsets"
        for dataset in self.registry_db.get_datasets():
            for variantSet in dataset.variant_sets:
                responseObject = self.sendGetObject(
                    path, variantSet.id, protocol.VariantSet)
                self.verifyVariantSetsEqual(responseObject, variantSet)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testGetVariantAnnotationSet(self):
        path = "/variantannotationsets"
        for dataset in self.registry_db.get_datasets():
            for variantSet in dataset.variant_sets:
                for vas in variantSet.variant_annotation_sets:
                    responseObject = self.sendGetObject(
                        path, vas.id, protocol.VariantAnnotationSet)
                    self.verifyVariantAnnotationSetsEqual(responseObject, vas)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    @unittest.skip("TODO fix simulated variant compound ID.")
    def testGetVariant(self):
        # get a variant from the search method
        referenceName = '1'
        start = 0
        request = protocol.SearchVariantsRequest()
        request.variant_set_id = str(self.variantSet.id)
        request.reference_name = referenceName
        request.start = start
        request.end = 2**16
        path = '/variants/search'
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchVariantsResponse)
        variants = responseData.variants[:10]

        # get 'the same' variant using the get method
        for variant in variants:
            print(variant.id)
            path = '/variants'
            responseObject = self.sendGetObject(
                path, variant.id, protocol.Variant)
            self.assertEqual(responseObject, variant)

    def testGetReferenceSet(self):
        path = "/referencesets"
        for referenceSet in self.registry_db.get_reference_sets():
            responseObject = self.sendGetObject(
                path, referenceSet.id, protocol.ReferenceSet)
            self.verifyReferenceSetsEqual(responseObject, referenceSet)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testGetReference(self):
        path = "/references"
        for referenceSet in self.registry_db.get_reference_sets():
            for reference in referenceSet.references:
                responseObject = self.sendGetObject(
                    path, reference.id, protocol.Reference)
                self.verifyReferencesEqual(responseObject, reference)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testGetCallSet(self):
        path = "/callsets"
        for dataset in self.registry_db.get_datasets():
            for variantSet in dataset.variant_sets:
                for callSet in variantSet.call_sets:
                    responseObject = self.sendGetObject(
                        path, callSet.id, protocol.CallSet)
                    self.verifyCallSetsEqual(responseObject, callSet)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testGetReadGroup(self):
        path = "/readgroups"
        for dataset in self.registry_db.get_datasets():
            for readGroupSet in dataset.read_group_sets:
                for readGroup in readGroupSet.read_groups:
                    responseObject = self.sendGetObject(
                        path, readGroup.id, protocol.ReadGroup)
                    self.verifyReadGroupsEqual(responseObject, readGroup)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testVariantsSearch(self):
        referenceName = '1'

        request = protocol.SearchVariantsRequest()
        request.reference_name = referenceName
        request.start = 0
        request.end = 0
        request.variant_set_id = str(self.variantSet.id)

        # Request windows is too small, no results
        path = '/variants/search'
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchVariantsResponse)
        self.assertEqual("", responseData.next_page_token)
        self.assertEqual(0, len(responseData.variants))

        # Larger request window, expect results
        request.end = 2 ** 16
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchVariantsResponse)
        self.assertTrue(protocol.validate(
            protocol.toJson(responseData), protocol.SearchVariantsResponse))
        self.assertGreater(len(responseData.variants), 0)

        # Verify all results are in the correct range, set and reference
        for variant in responseData.variants:
            self.assertGreaterEqual(variant.start, 0)
            self.assertLessEqual(variant.end, 2 ** 16)
            self.assertEqual(variant.variant_set_id, str(self.variantSet.id))
            self.assertEqual(variant.reference_name, referenceName)

        # TODO: Add more useful test scenarios, including some covering
        # pagination behavior.

    def testVariantAnnotationSetsSearch(self):
        path = '/variantannotationsets/search'
        for variant_set in self.registry_db.get_variant_sets():
            request = protocol.SearchVariantAnnotationSetsRequest()
            request.variant_set_id = str(variant_set.id)
            self.verifySearchMethod(
                request, path, protocol.SearchVariantAnnotationSetsResponse,
                variant_set.variant_annotation_sets,
                self.verifyVariantAnnotationSetsEqual)
        for badId in self.getBadIds():
            request = protocol.SearchVariantAnnotationSetsRequest()
            request.variant_set_id = badId
            self.verifySearchMethodFails(request, path)

    def testVariantAnnotationsSearch(self):
        self.assertIsNotNone(self.variantAnnotationSet)

        request = protocol.SearchVariantAnnotationsRequest()
        # TODO split these into separate tests, and factor out the duplicated
        # code.

        path = '/variantannotations/search'
        request.start = 0
        request.end = 10000
        request.page_size = 1
        request.reference_name = "1"
        request.variant_annotation_set_id = str(self.variantAnnotationSet.id)
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        responseData = protocol.fromJson(response.data, protocol.
                                         SearchVariantAnnotationsResponse)
        self.assertGreater(len(responseData.variant_annotations), 0)
        self.assertIsNotNone(
            responseData.next_page_token,
            "Expected more than one page of results")

        request = protocol.SearchVariantAnnotationsRequest()
        request.variant_annotation_set_id = str(self.variantAnnotationSet.id)
        request.start = 0
        request.end = 10
        request.reference_name = "1"

        request.effects.add().id = "ThisIsNotAnEffect"

        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        responseData = protocol.fromJson(response.data, protocol.
                                         SearchVariantAnnotationsResponse)
        self.assertEquals(
            len(responseData.variant_annotations), 0,
            "There should be no results for a nonsense effect")

        request = protocol.SearchVariantAnnotationsRequest()
        request.variant_annotation_set_id = str(self.variantAnnotationSet.id)
        request.start = 0
        request.end = 10
        request.reference_name = "1"
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        responseData = protocol.fromJson(response.data, protocol.
                                         SearchVariantAnnotationsResponse)
        self.assertGreater(len(responseData.variant_annotations), 0)
        for ann in responseData.variant_annotations:
            self.assertGreater(
                len(ann.transcript_effects), 0,
                ("When no effects are requested ensure "
                    "some transcript effects are still present"))

        request = protocol.SearchVariantAnnotationsRequest()
        request.variant_annotation_set_id = str(self.variantAnnotationSet.id)
        request.start = 0
        request.end = 5
        request.reference_name = "1"
        request.effects.add().id = "SO:0001627"
        request.effects.add().id = "B4DID"
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        responseData = protocol.fromJson(response.data, protocol.
                                         SearchVariantAnnotationsResponse)
        responseLength = len(responseData.variant_annotations)
        self.assertGreater(
            responseLength, 0,
            "There should be some results for a known effect")
        for ann in responseData.variant_annotations:
            effectPresent = False
            for effect in ann.transcript_effects:
                for featureType in effect.effects:
                    if featureType.id in map(
                            lambda e: e.id, request.effects):
                        effectPresent = True
            self.assertEquals(
                True, effectPresent,
                "The ontology term should appear at least once")

        request = protocol.SearchVariantAnnotationsRequest()
        request.variant_annotation_set_id = str(self.variantAnnotationSet.id)
        request.start = 0
        request.end = 5
        request.reference_name = "1"
        request.effects.add().id = "B4DID"
        request.effects.add().id = "SO:0001627"
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        responseData = protocol.fromJson(response.data, protocol.
                                         SearchVariantAnnotationsResponse)
        self.assertEqual(
            len(responseData.variant_annotations),
            responseLength,
            "Order shall not affect results")
        for ann in responseData.variant_annotations:
            effectPresent = False
            for effect in ann.transcript_effects:
                for featureType in effect.effects:
                    if featureType.id in map(
                            lambda e: e.id, request.effects):
                        effectPresent = True
            self.assertEquals(
                True,
                effectPresent,
                "The ontology term should appear at least once")

        request = protocol.SearchVariantAnnotationsRequest()
        request.variant_annotation_set_id = str(self.variantAnnotationSet.id)
        request.start = 0
        request.end = 5
        request.reference_name = "1"
        request.effects.add().id = "SO:0001627"
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        responseData = protocol.fromJson(response.data, protocol.
                                         SearchVariantAnnotationsResponse)
        self.assertGreater(len(responseData.variant_annotations), 0,
                           "There should be some results for a good effect ID")
        for ann in responseData.variant_annotations:
            effectPresent = False
            for effect in ann.transcript_effects:
                for featureType in effect.effects:
                    if featureType.id in map(
                            lambda e: e.id, request.effects):
                        effectPresent = True
            self.assertEquals(True, effectPresent,
                              "The ontology term should appear at least once")

        request = protocol.SearchVariantAnnotationsRequest()
        request.variant_annotation_set_id = str(self.variantAnnotationSet.id)
        request.start = 0
        request.end = 10
        request.reference_name = "1"
        request.effects.add().id = "SO:0001627"
        request.effects.add().id = "SO:0001791"
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        responseData = protocol.fromJson(response.data, protocol.
                                         SearchVariantAnnotationsResponse)
        self.assertGreater(len(responseData.variant_annotations), 0)

    def testGetFeatureSet(self):
        path = "/featuresets"
        for dataset in self.registry_db.get_datasets():
            for featureSet in dataset.feature_sets:
                responseObject = self.sendGetObject(
                    path, featureSet.id, protocol.FeatureSet)
                self.verifyFeatureSetsEqual(responseObject, featureSet)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testFeatureSetsSearch(self):
        path = '/featuresets/search'
        for dataset in self.registry_db.get_datasets():
            featureSets = list(dataset.feature_sets)
            request = protocol.SearchFeatureSetsRequest()
            request.dataset_id = str(dataset.id)
            self.verifySearchMethod(
                request, path, protocol.SearchFeatureSetsResponse, featureSets,
                self.verifyFeatureSetsEqual)
        for badId in self.getBadIds():
            request = protocol.SearchFeatureSetsRequest()
            request.dataset_id = badId
            self.verifySearchMethodFails(request, path)

    @unittest.skip("TODO simulatoed feature Compound IDs")
    def testGetFeature(self):
        dataset = self.registry_db.get_datasets()[0]
        featureSet = dataset.feature_sets[0]
        request = protocol.SearchFeaturesRequest()
        request.feature_set_id = str(featureSet.id)
        request.reference_name = "chr1"
        request.start = 0
        request.end = 10
        path = '/features/search'
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchFeaturesResponse)
        features = responseData.features
        self.assertGreater(len(features), 0)

        # get 'the same' feature using the get method
        for feature in features:
            path = '/features'
            responseObject = self.sendGetObject(
                path, feature.id, protocol.Feature)
            self.verifyFeaturesEquivalent(responseObject, feature)

    def testFeaturesSearch(self):
        dataset = self.registry_db.get_datasets()[0]
        featureSet = dataset.feature_sets[0]
        referenceName = featureSet.reference_set.references[0].name

        request = protocol.SearchFeaturesRequest()
        request.reference_name = referenceName
        request.feature_set_id = str(featureSet.id)

        path = "/features/search"
        request.start = 0
        request.end = 10
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchFeaturesResponse)
        self.assertEqual(len(responseData.features), 10)

        # Verify all results are in the correct range, set and reference
        for feature in responseData.features:
            self.assertGreaterEqual(feature.start, 0)
            self.assertGreaterEqual(feature.end, feature.start)
            self.assertEqual(feature.feature_set_id, str(featureSet.id))
            self.assertEqual(feature.reference_name, referenceName)

    def testListReferenceBases(self):
        for referenceSet in self.registry_db.get_reference_sets():
            for reference in referenceSet.references:
                id_ = reference.id
                length = reference.length
                sequence = reference.run_get_bases(0, length)
                # fetch the bases
                args = protocol.ListReferenceBasesRequest()
                response = self.sendListReferenceBasesRequest(id_, args)
                self.assertEqual(response.sequence, sequence)
                # Try some simple slices.
                ranges = [(0, length), (0, 1), (length - 1, length)]
                for start, end in ranges:
                    args = protocol.ListReferenceBasesRequest()
                    args.start, args.end = start, end
                    response = self.sendListReferenceBasesRequest(id_, args)
                    self.assertEqual(response.sequence, sequence[start:end])
                    self.assertEqual("", response.next_page_token)
                    self.assertEqual(response.offset, start)

    def testListReferenceBasesErrors(self):
        for badId in self.getBadIds():
            path = '/references/{}/bases'.format(badId)
            response = self.app.get(path)
            self.assertEqual(response.status_code, 404)
        path = '/references/{}/bases'.format(self.reference.id)
        length = self.reference.length
        badRanges = [(-1, 0), (-1, -1), (length, 1), (0, length + 1)]
        for start, end in badRanges:
            args = protocol.ListReferenceBasesRequest()
            args.start, args.end = start, end
            response = self.app.get(
                path, query_string=protocol.toJsonDict(args))
            self.assertEqual(response.status_code, 416)

    def testListReferenceBasesPaging(self):
        id_ = self.reference.id
        length = self.reference.length
        completeSequence = self.reference.run_get_bases(0, length)
        for start, end in [(0, length), (5, 10), (length // 2, length)]:
            sequence = completeSequence[start: end]
            for pageSize in [1, 2, length - 1]:
                self.backend.setMaxResponseLength(pageSize)
                args = protocol.ListReferenceBasesRequest()
                args.start, args.end = start, end
                response = self.sendListReferenceBasesRequest(id_, args)
                self.assertEqual(response.sequence, sequence[:pageSize])
                self.assertEqual(response.offset, start)
                sequenceFragments = [response.sequence]
                while response.next_page_token is not "":
                    args = protocol.ListReferenceBasesRequest()
                    args.page_token = response.next_page_token
                    args.start, args.end = start, end
                    response = self.sendListReferenceBasesRequest(id_, args)
                    self.assertGreater(len(response.sequence), 0)
                    sequenceFragments.append(response.sequence)
                    offset = response.offset
                    self.assertEqual(
                        response.sequence,
                        completeSequence[
                            offset: offset + len(response.sequence)])

                self.assertEqual("".join(sequenceFragments), sequence)

    def testReads(self):
        path = '/reads/search'
        for dataset in self.registry_db.get_datasets():
            for readGroupSet in dataset.read_group_sets:
                referenceSet = readGroupSet.reference_set
                for reference in referenceSet.references:
                    for readGroup in readGroupSet.read_groups:
                        # search reads
                        request = protocol.SearchReadsRequest()
                        request.read_group_ids.append(str(readGroup.id))
                        request.reference_id = str(reference.id)
                        request.start = 0
                        request.end = 3
                        responseData = self.sendSearchRequest(
                            path, request, protocol.SearchReadsResponse)
                        alignments = responseData.alignments
                        self.assertGreater(len(alignments), 0)
                        for alignment in alignments:
                            # TODO more tests here: this is very weak.
                            self.assertEqual(
                                alignment.read_group_id, str(readGroup.id))

    def testUnsupportedReadOperations(self):
        path = '/reads/search'

        # unmapped Reads
        request = protocol.SearchReadsRequest()
        request.read_group_ids.append(str(self.readGroup.id))
        request.reference_id = ""
        self.verifySearchMethodNotSupported(request, path)

        # multiple ReadGroupSets set mismatch
        request = protocol.SearchReadsRequest()
        readGroupSets = list(self.dataset.read_group_sets)
        self.assertGreater(len(readGroupSets), 1)
        for readGroupSet in readGroupSets:
            request.read_group_ids.extend(
                str(rg.id) for rg in readGroupSet.read_groups)
        request.reference_id = str(self.reference.id)
        response = self.sendJsonPostRequest(path, protocol.toJson(request))
        self.assertEqual(400, response.status_code)

    def testReadsMultipleReadGroupSets(self):
        path = '/reads/search'
        readGroupIds = [
            str(readGroup.id) for readGroup in
            self.readGroupSet.read_groups]
        referenceId = str(self.reference.id)
        request = protocol.SearchReadsRequest()
        request.read_group_ids.extend(readGroupIds)
        request.reference_id = referenceId
        request.start = 0
        request.end = 10
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchReadsResponse)
        alignments = responseData.alignments
        self.assertGreater(len(alignments), 0)
        for alignment in alignments:
            self.assertIn(alignment.read_group_id, readGroupIds)

    @unittest.skip("Skipping until ReadGroupSet bioSamplesId semantics fixed")
    def testBioSamplesFromReadGroupSets(self):
        path = 'readgroupsets/search'
        dataset = self.registry_db.get_datasets()[0]
        # get all the read group sets
        request = protocol.SearchReadGroupSetsRequest()
        request.dataset_id = str(dataset.id)
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchReadGroupSetsResponse)
        # go through each read group
        bioSamplesRgs = []
        ran = False
        for rgs in responseData.read_group_sets:
            for rg in rgs.read_groups:
                ran = True
                # request biosample record
                if rg.bio_sample_id:
                    bioSample = self.sendGetObject(
                        'biosamples',
                        rg.bio_sample_id,
                        protocol.BioSample)
                    bioSamplesRgs.append((bioSample.id, rgs.id))
                    self.assertNotEqual(
                        None, bioSample,
                        "A BioSample should exist for reach read")
        self.assertTrue(ran)
        # search reads by biosample
        ran = False
        for bsId, rgsId in bioSamplesRgs:
            request = protocol.SearchReadGroupSetsRequest()
            request.dataset_id = str(dataset.id)
            request.bio_sample_id = bsId
            request.name = "A BAD NAME"
            request.page_size = 1
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchReadGroupSetsResponse)
            self.assertEquals(
                len(responseData.read_group_sets), 0,
                "A good biosample ID and bad name should return 0")
            request = protocol.SearchReadGroupSetsRequest()
            request.dataset_id = str(dataset.id)
            request.bio_sample_id = bsId
            responseData = self.sendSearchRequest(
                path, request, protocol.SearchReadGroupSetsResponse)
            for rgs in responseData.read_group_sets:
                for rg in rgs.read_groups:
                    ran = True
                    self.assertEqual(
                        rg.bio_sample_id, bsId,
                        "Only read groups matching the BioSample ID")
        self.assertTrue(ran)

    def testBioSamplesFromCallSets(self):
        path = 'callsets/search'
        dataset = self.registry_db.get_datasets()[0]
        variantSet = dataset.variant_sets[0]
        callSet = variantSet.call_sets[0]
        request = protocol.SearchCallSetsRequest()
        request.variant_set_id = str(variantSet.id)
        request.bio_sample_id = "A BAD ID"
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchCallSetsResponse)
        self.assertEqual(len(responseData.call_sets), 0)

        request = protocol.SearchCallSetsRequest()
        request.variant_set_id = str(variantSet.id)
        request.bio_sample_id = str(callSet.bio_sample_id)
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchCallSetsResponse)
        self.assertGreater(len(responseData.call_sets), 0)
        for cs in responseData.call_sets:
            self.assertEqual(cs.bio_sample_id, request.bio_sample_id)

        request = protocol.SearchCallSetsRequest()
        request.variant_set_id = str(variantSet.id)
        request.bio_sample_id = str(callSet.bio_sample_id)
        request.name = "A BAD NAME"
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchCallSetsResponse)
        self.assertEqual(
            len(responseData.call_sets), 0,
            "None should be returned")

    def testBioSamplesSearch(self):
        path = 'biosamples/search'
        dataset = self.registry_db.get_datasets()[0]
        bio_sample = dataset.bio_samples[0]
        request = protocol.SearchBioSamplesRequest()
        request.name = "BAD NAME"
        request.dataset_id = str(dataset.id)
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchBioSamplesResponse)
        self.assertEqual(
            len(responseData.biosamples), 0,
            "A bad name should return none")
        request = protocol.SearchBioSamplesRequest()
        request.name = bio_sample.name
        request.dataset_id = str(dataset.id)
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchBioSamplesResponse)
        # Currently always returns a singleton
        self.assertGreater(
            len(responseData.biosamples), 0,
            "A good name should return some")

        request = protocol.SearchBioSamplesRequest()
        request.individual_id = "BAD ID"
        request.dataset_id = str(dataset.id)
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchBioSamplesResponse)
        self.assertEqual(
            len(responseData.biosamples), 0,
            "A bad individual ID should return none")

        request = protocol.SearchIndividualsRequest()
        request.dataset_id = str(dataset.id)
        responseData = self.sendSearchRequest(
            "individuals/search", request, protocol.SearchIndividualsResponse)
        self.assertGreater(
            len(responseData.individuals), 0,
            "Some individuals should be returned")

        individualId = responseData.individuals[0].id

        request = protocol.SearchBioSamplesRequest()
        request.individual_id = individualId
        request.dataset_id = str(dataset.id)
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchBioSamplesResponse)
        self.assertGreater(
            len(responseData.biosamples), 0,
            "A good individual ID should return some")

        request = protocol.SearchBioSamplesRequest()
        request.individual_id = individualId
        request.page_size = 1
        request.dataset_id = str(dataset.id)
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchBioSamplesResponse)
        self.assertIsNotNone(
            responseData.next_page_token,
            "More than one page should be returned")
        request.page_token = responseData.next_page_token
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchBioSamplesResponse)
        self.assertEqual(
            responseData.biosamples[0].individual_id,
            individualId,
            "Results on the second page should match")

    def testSearchIndividuals(self):
        path = 'individuals/search'
        dataset = self.registry_db.get_datasets()[0]
        individual = dataset.individuals[0]
        request = protocol.SearchIndividualsRequest()
        request.name = "BAD NAME"
        request.dataset_id = str(dataset.id)
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchIndividualsResponse)
        self.assertEqual(
            len(responseData.individuals), 0,
            "A bad individual name should return none")
        request = protocol.SearchIndividualsRequest()
        request.name = individual.name
        request.dataset_id = str(dataset.id)
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchIndividualsResponse)
        self.assertGreater(
            len(responseData.individuals), 0,
            "A good individual name should return some")

    def testGetIndividual(self):
        path = "/individuals"
        for dataset in self.registry_db.get_datasets():
            for individual in dataset.individuals:
                responseObject = self.sendGetObject(
                    path, individual.id, protocol.Individual)
                self.verifyIndividualsEqual(responseObject, individual)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testGetBioSample(self):
        path = "/biosamples"
        for dataset in self.registry_db.get_datasets():
            for bioSample in dataset.bio_samples:
                responseObject = self.sendGetObject(
                    path, bioSample.id, protocol.BioSample)
                self.verifyBioSamplesEqual(responseObject, bioSample)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)
