"""
End-to-end tests for the simulator configuration. Sets up a server with
the backend, sends some basic queries to that server and verifies results
are as expected.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import logging
import random
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.references as references
import ga4gh.datamodel.variants as variants
import ga4gh.frontend as frontend
import ga4gh.protocol as protocol


class TestSimulatedStack(unittest.TestCase):
    """
    Tests the full stack for the Simulated backend by using the Flask
    testing client.
    """
    @classmethod
    def setUpClass(cls):
        # silence usually unhelpful CORS log
        logging.getLogger('ga4gh.frontend.cors').setLevel(logging.CRITICAL)
        # Set the random seed to make tests reproducible.
        random.seed(1)
        config = {
            "DATA_SOURCE": "simulated://",
            "SIMULATED_BACKEND_RANDOM_SEED": 1111,
            "SIMULATED_BACKEND_NUM_CALLS": 5,
            "SIMULATED_BACKEND_VARIANT_DENSITY": 1.0,
            "SIMULATED_BACKEND_NUM_VARIANT_SETS": 4,
            "SIMULATED_BACKEND_NUM_REFERENCE_SETS": 3,
            "SIMULATED_BACKEND_NUM_REFERENCES_PER_REFERENCE_SET": 4,
            "SIMULATED_BACKEND_NUM_ALIGNMENTS_PER_READ_GROUP": 5
        }
        frontend.reset()
        frontend.configure(
            baseConfig="TestConfig", extraConfig=config)
        cls.app = frontend.app.test_client()

    @classmethod
    def tearDownClass(cls):
        cls.app = None

    def setUp(self):
        self.backend = frontend.app.backend
        self.backend.setMaxResponseLength(10000)

    def getBadIds(self):
        """
        Returns a list of IDs that should not exist in the server and should
        raise a 404 error.
        """
        return ["", "1234:", "x"*100, ":", ":xx", "::", ":::", "::::"]

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
        response = self.sendJsonPostRequest(path, request.toJsonString())
        self.assertEqual(200, response.status_code)
        responseData = responseClass.fromJsonString(response.data)
        self.assertTrue(responseData.validate(responseData.toJsonDict()))
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
        obj = responseClass.fromJsonString(response.data)
        self.assertTrue(responseClass.validate(obj.toJsonDict()))
        return obj

    def sendListReferenceBasesRequest(self, id_, request):
        """
        Sends a ListReferenceBasesRequest and parses the result into a
        ListReferenceBasesResponse.
        """
        path = '/references/{}/bases'.format(id_)
        response = self.app.get(path, query_string=request.toJsonDict())
        self.assertEqual(response.status_code, 200)
        obj = protocol.ListReferenceBasesResponse.fromJsonString(
            response.data)
        self.assertTrue(obj.validate(obj.toJsonDict()))
        return obj

    def verifyVariantSetsEqual(self, gaVariantSet, variantSet):
        dataset = variantSet.getParentContainer()
        self.assertEqual(gaVariantSet.id, variantSet.getId())
        self.assertEqual(gaVariantSet.datasetId, dataset.getId())
        self.assertEqual(gaVariantSet.name, variantSet.getLocalId())
        # TODO verify the metadata and other attributes.

    def verifyCallSetsEqual(self, gaCallSet, callSet):
        variantSet = callSet.getParentContainer()
        self.assertEqual(gaCallSet.id, callSet.getId())
        self.assertEqual(gaCallSet.name, callSet.getLocalId())
        self.assertEqual(gaCallSet.variantSetIds, [variantSet.getId()])
        # TODO add some simulated info and check

    def verifyReadGroupSetsEqual(self, gaReadGroupSet, readGroupSet):
        dataset = readGroupSet.getParentContainer()
        self.assertEqual(gaReadGroupSet.id, readGroupSet.getId())
        self.assertEqual(gaReadGroupSet.datasetId, dataset.getId())
        self.assertEqual(gaReadGroupSet.name, readGroupSet.getLocalId())
        self.assertEqual(
            len(gaReadGroupSet.readGroups), len(readGroupSet.getReadGroups()))
        for gaReadGroup, readGroup in zip(
                gaReadGroupSet.readGroups, readGroupSet.getReadGroups()):
            self.verifyReadGroupsEqual(gaReadGroup, readGroup)

    def verifyReadGroupsEqual(self, gaReadGroup, readGroup):
        self.assertEqual(gaReadGroup.id, readGroup.getId())

    def verifyDatasetsEqual(self, gaDataset, dataset):
        self.assertEqual(gaDataset.id, dataset.getId())
        self.assertEqual(gaDataset.name, dataset.getLocalId())
        self.assertEqual(gaDataset.description, dataset.getDescription())

    def verifyReferenceSetsEqual(self, gaReferenceSet, referenceSet):
        self.assertEqual(gaReferenceSet.id, referenceSet.getId())
        self.assertEqual(
            gaReferenceSet.md5checksum, referenceSet.getMd5Checksum())
        self.assertEqual(
            gaReferenceSet.ncbiTaxonId, referenceSet.getNcbiTaxonId())
        self.assertEqual(
            gaReferenceSet.assemblyId, referenceSet.getAssemblyId())
        self.assertEqual(
            gaReferenceSet.sourceURI, referenceSet.getSourceUri())
        self.assertEqual(
            gaReferenceSet.sourceAccessions,
            referenceSet.getSourceAccessions())
        self.assertEqual(
            gaReferenceSet.isDerived, referenceSet.getIsDerived())
        self.assertEqual(
            gaReferenceSet.name, referenceSet.getLocalId())

    def verifyReferencesEqual(self, gaReference, reference):
        self.assertEqual(gaReference.id, reference.getId())
        self.assertEqual(gaReference.name, reference.getName())
        self.assertEqual(gaReference.length, reference.getLength())
        self.assertEqual(gaReference.md5checksum, reference.getMd5Checksum())
        self.assertEqual(gaReference.ncbiTaxonId, reference.getNcbiTaxonId())
        self.assertEqual(gaReference.sourceURI, reference.getSourceUri())
        self.assertEqual(
            gaReference.sourceAccessions, reference.getSourceAccessions())
        self.assertEqual(gaReference.isDerived, reference.getIsDerived())
        self.assertEqual(
            gaReference.sourceDivergence, reference.getSourceDivergence())

    def verifySearchMethod(
            self, request, path, responseClass, objects, objectVerifier):
        """
        Verifies that the specified search request operates correctly
        and returns all the speficied objects. The specified verifier
        function checks that all the returned objects are equivalent
        to their datamodel counterparts.
        """
        request.pageSize = len(objects)
        self.assertGreater(request.pageSize, 0)
        responseData = self.sendSearchRequest(path, request, responseClass)
        self.assertIsNone(responseData.nextPageToken)
        responseList = getattr(responseData, responseClass.getValueListName())
        self.assertEqual(len(objects), len(responseList))
        for gaObject, datamodelObject in zip(responseList, objects):
            objectVerifier(gaObject, datamodelObject)

    def verifySearchResultsEmpty(self, request, path, responseClass):
        """
        Verifies that we get a successful response with an empty list of
        results.
        """
        responseData = self.sendSearchRequest(path, request, responseClass)
        self.assertIsNone(responseData.nextPageToken)
        responseList = getattr(responseData, responseClass.getValueListName())
        self.assertEqual(0, len(responseList))

    def assertObjectNotFound(self, response):
        """
        Checks that the specified response contains a search failure.
        """
        self.assertEqual(404, response.status_code)
        error = protocol.GAException.fromJsonString(response.data)
        self.assertTrue(error.validate(error.toJsonDict()))
        self.assertGreater(error.errorCode, 0)
        self.assertGreater(len(error.message), 0)

    def verifySearchMethodFails(self, request, path):
        """
        Verify that the specified search request fails with a 404.
        """
        response = self.sendJsonPostRequest(path, request.toJsonString())
        self.assertObjectNotFound(response)

    def verifyGetMethodFails(self, path, id_):
        """
        Verifies the specified GET request failes with a 404.
        """
        response = self.sendObjectGetRequest(path, id_)
        self.assertObjectNotFound(response)

    def testGetDataset(self):
        path = "/datasets"
        for dataset in self.backend.getDatasets():
            responseObject = self.sendGetObject(
                path, dataset.getId(), protocol.Dataset)
            self.verifyDatasetsEqual(responseObject, dataset)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testDatasetsSearch(self):
        request = protocol.SearchDatasetsRequest()
        datasets = self.backend.getDatasets()
        path = '/datasets/search'
        self.verifySearchMethod(
            request, path, protocol.SearchDatasetsResponse, datasets,
            self.verifyDatasetsEqual)

    def testVariantSetsSearch(self):
        path = '/variantsets/search'
        for dataset in self.backend.getDatasets():
            variantSets = dataset.getVariantSets()
            request = protocol.SearchVariantSetsRequest()
            request.datasetId = dataset.getId()
            self.verifySearchMethod(
                request, path, protocol.SearchVariantSetsResponse, variantSets,
                self.verifyVariantSetsEqual)
        for badId in self.getBadIds():
            request = protocol.SearchVariantSetsRequest()
            request.datasetId = badId
            self.verifySearchMethodFails(request, path)

    def testCallSetsSearch(self):
        path = '/callsets/search'
        for dataset in self.backend.getDatasets():
            for variantSet in dataset.getVariantSets():
                callSets = variantSet.getCallSets()
                self.assertGreater(len(callSets), 0)
                request = protocol.SearchCallSetsRequest()
                request.variantSetId = variantSet.getId()
                self.verifySearchMethod(
                    request, path, protocol.SearchCallSetsResponse, callSets,
                    self.verifyCallSetsEqual)
                # Check if we can search for the callset with a good name.
                for callSet in callSets:
                    request = protocol.SearchCallSetsRequest()
                    request.variantSetId = variantSet.getId()
                    request.name = callSet.getLocalId()
                    self.verifySearchMethod(
                        request, path, protocol.SearchCallSetsResponse,
                        [callSet], self.verifyCallSetsEqual)
                # Check if we can search for the callset with a bad name.
                for badId in self.getBadIds():
                    request = protocol.SearchCallSetsRequest()
                    request.variantSetId = variantSet.getId()
                    request.name = badId
                    self.verifySearchResultsEmpty(
                        request, path, protocol.SearchCallSetsResponse)
        # Check for searches within missing variantSets.
        for badId in self.getBadIds():
            request = protocol.SearchCallSetsRequest()
            request.variantSetId = badId
            self.verifySearchMethodFails(request, path)

    def testReadGroupSetsSearch(self):
        path = '/readgroupsets/search'
        for dataset in self.backend.getDatasets():
            readGroupSets = dataset.getReadGroupSets()
            request = protocol.SearchReadGroupSetsRequest()
            request.datasetId = dataset.getId()
            self.verifySearchMethod(
                request, path, protocol.SearchReadGroupSetsResponse,
                readGroupSets, self.verifyReadGroupSetsEqual)
            # Check if we can search for the readGroupSet with a good name.
            for readGroupSet in readGroupSets:
                request = protocol.SearchReadGroupSetsRequest()
                request.datasetId = dataset.getId()
                request.name = readGroupSet.getLocalId()
                self.verifySearchMethod(
                    request, path, protocol.SearchReadGroupSetsResponse,
                    [readGroupSet], self.verifyReadGroupSetsEqual)
            # Check if we can search for the readGroupSet with a bad name.
            for badId in self.getBadIds():
                request = protocol.SearchReadGroupSetsRequest()
                request.datasetId = dataset.getId()
                request.name = badId
                self.verifySearchResultsEmpty(
                    request, path, protocol.SearchReadGroupSetsResponse)
        for badId in self.getBadIds():
            request = protocol.SearchReadGroupSetsRequest()
            request.datasetId = badId
            self.verifySearchMethodFails(request, path)

    def testReferenceSetsSearch(self):
        request = protocol.SearchReferenceSetsRequest()
        referenceSets = self.backend.getReferenceSets()
        path = '/referencesets/search'
        self.verifySearchMethod(
            request, path, protocol.SearchReferenceSetsResponse, referenceSets,
            self.verifyReferenceSetsEqual)

    def testReferencesSearch(self):
        path = '/references/search'
        for referenceSet in self.backend.getReferenceSets():
            references = referenceSet.getReferences()
            request = protocol.SearchReferencesRequest()
            request.referenceSetId = referenceSet.getId()
            self.verifySearchMethod(
                request, path, protocol.SearchReferencesResponse, references,
                self.verifyReferencesEqual)
        for badId in self.getBadIds():
            request = protocol.SearchReferencesRequest()
            request.referenceSetId = badId
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
            request.md5checksum = obj.getMd5Checksum()
            self.verifySearchMethod(
                request, path, responseClass, [obj], objectVerifier)
            request.md5checksum = None
            request.accession = obj.getSourceAccessions()[0]
            self.verifySearchMethod(
                request, path, responseClass, [obj], objectVerifier)
            request.accession = None
            if hasAssemblyId:
                request.assemblyId = obj.getAssemblyId()
                self.verifySearchMethod(
                    request, path, responseClass, [obj], objectVerifier)
                request.assemblyId = None
            # Now check one good value and some bad values.
            request.md5checksum = obj.getMd5Checksum()
            badAccessions = [
                "no such accession", objectList[0].getSourceAccessions()[0]]
            for accession in badAccessions:
                request.accession = accession
                self.verifySearchResultsEmpty(request, path, responseClass)
            request.accession = None
            if hasAssemblyId:
                badAssemblyIds = [
                    "no such asssembly", objectList[0].getAssemblyId()]
                for assemblyId in badAssemblyIds:
                    request.assemblyId = assemblyId
                    self.verifySearchResultsEmpty(request, path, responseClass)
                request.assemblyId = None

    def testReferencesSearchFilters(self):
        path = '/references/search'
        for referenceSet in self.backend.getReferenceSets():

            def requestFactory():
                request = protocol.SearchReferencesRequest()
                request.referenceSetId = referenceSet.getId()
                return request
            self.verifyReferenceSearchFilters(
                referenceSet.getReferences(), False, path, requestFactory,
                protocol.SearchReferencesResponse, self.verifyReferencesEqual)

    def testReferenceSetsSearchFilters(self):
        path = '/referencesets/search'

        def requestFactory():
            return protocol.SearchReferenceSetsRequest()
        self.verifyReferenceSearchFilters(
            self.backend.getReferenceSets(), True, path, requestFactory,
            protocol.SearchReferenceSetsResponse,
            self.verifyReferenceSetsEqual)

    def testGetVariantSet(self):
        path = "/variantsets"
        for dataset in self.backend.getDatasets():
            for variantSet in dataset.getVariantSets():
                responseObject = self.sendGetObject(
                    path, variantSet.getId(), protocol.VariantSet)
                self.verifyVariantSetsEqual(responseObject, variantSet)
            for badId in self.getBadIds():
                variantSet = variants.AbstractVariantSet(dataset, badId)
                self.verifyGetMethodFails(path, variantSet.getId())
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testGetVariant(self):
        # get a variant from the search method
        referenceName = '1'
        start = 0
        dataset = self.backend.getDatasets()[0]
        variantSet = dataset.getVariantSets()[0]
        request = protocol.SearchVariantsRequest()
        request.variantSetId = variantSet.getId()
        request.referenceName = referenceName
        request.start = start
        request.end = 2**16
        path = '/variants/search'
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchVariantsResponse)
        variants = responseData.variants[:10]

        # get 'the same' variant using the get method
        for variant in variants:
            path = '/variants'
            responseObject = self.sendGetObject(
                path, variant.id, protocol.Variant)
            self.assertEqual(responseObject, variant)

    def testGetReferenceSet(self):
        path = "/referencesets"
        for referenceSet in self.backend.getReferenceSets():
            responseObject = self.sendGetObject(
                path, referenceSet.getId(), protocol.ReferenceSet)
            self.verifyReferenceSetsEqual(responseObject, referenceSet)
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testGetReference(self):
        path = "/references"
        for referenceSet in self.backend.getReferenceSets():
            for reference in referenceSet.getReferences():
                responseObject = self.sendGetObject(
                    path, reference.getId(), protocol.Reference)
                self.verifyReferencesEqual(responseObject, reference)
            for badId in self.getBadIds():
                referenceSet = references.AbstractReferenceSet(badId)
                self.verifyGetMethodFails(path, referenceSet.getId())
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testGetCallSet(self):
        path = "/callsets"
        for dataset in self.backend.getDatasets():
            for variantSet in dataset.getVariantSets():
                for callSet in variantSet.getCallSets():
                    responseObject = self.sendGetObject(
                        path, callSet.getId(), protocol.CallSet)
                    self.verifyCallSetsEqual(responseObject, callSet)
                for badId in self.getBadIds():
                    callSet = variants.CallSet(variantSet, badId)
                    self.verifyGetMethodFails(path, callSet.getId())
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testGetReadGroup(self):
        path = "/readgroups"
        for dataset in self.backend.getDatasets():
            for readGroupSet in dataset.getReadGroupSets():
                for readGroup in readGroupSet.getReadGroups():
                    responseObject = self.sendGetObject(
                        path, readGroup.getId(), protocol.ReadGroup)
                    self.verifyReadGroupsEqual(responseObject, readGroup)
                for badId in self.getBadIds():
                    readGroup = reads.AbstractReadGroup(readGroupSet, badId)
                    self.verifyGetMethodFails(path, readGroup.getId())
            for badId in self.getBadIds():
                readGroupSet = reads.AbstractReadGroupSet(dataset, badId)
                self.verifyGetMethodFails(path, readGroupSet.getId())
        for badId in self.getBadIds():
            self.verifyGetMethodFails(path, badId)

    def testVariantsSearch(self):
        dataset = self.backend.getDatasets()[0]
        variantSet = dataset.getVariantSets()[0]
        referenceName = '1'

        request = protocol.SearchVariantsRequest()
        request.referenceName = referenceName
        request.start = 0
        request.end = 0
        request.variantSetId = variantSet.getId()

        # Request windows is too small, no results
        path = '/variants/search'
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchVariantsResponse)
        self.assertIsNone(responseData.nextPageToken)
        self.assertEqual([], responseData.variants)

        # Larger request window, expect results
        request.end = 2 ** 16
        responseData = self.sendSearchRequest(
            path, request, protocol.SearchVariantsResponse)
        self.assertTrue(protocol.SearchVariantsResponse.validate(
            responseData.toJsonDict()))
        self.assertGreater(len(responseData.variants), 0)

        # Verify all results are in the correct range, set and reference
        for variant in responseData.variants:
            self.assertGreaterEqual(variant.start, 0)
            self.assertLessEqual(variant.end, 2 ** 16)
            self.assertEqual(variant.variantSetId, variantSet.getId())
            self.assertEqual(variant.referenceName, referenceName)

        # TODO: Add more useful test scenarios, including some covering
        # pagination behavior.

    def testListReferenceBases(self):
        for referenceSet in self.backend.getReferenceSets():
            for reference in referenceSet.getReferences():
                id_ = reference.getId()
                length = reference.getLength()
                sequence = reference.getBases(0, length)
                # fetch the bases
                args = protocol.ListReferenceBasesRequest()
                response = self.sendListReferenceBasesRequest(id_, args)
                self.assertEqual(response.sequence, sequence)
                # Try some simple slices.
                ranges = [(0, length), (0, 1), (length - 1, length), (0, 0)]
                for start, end in ranges:
                    args = protocol.ListReferenceBasesRequest()
                    args.start, args.end = start, end
                    response = self.sendListReferenceBasesRequest(id_, args)
                    self.assertEqual(response.sequence, sequence[start:end])
                    self.assertIsNone(response.nextPageToken)
                    self.assertEqual(response.offset, start)

    def testListReferenceBasesErrors(self):
        referenceSet = self.backend.getReferenceSets()[0]
        for badId in self.getBadIds():
            path = '/references/{}/bases'.format(badId)
            response = self.app.get(path)
            self.assertEqual(response.status_code, 404)
            reference = references.AbstractReference(referenceSet, badId)
            path = '/references/{}/bases'.format(reference.getId())
            response = self.app.get(path)
            self.assertEqual(response.status_code, 404)
        reference = referenceSet.getReferences()[0]
        path = '/references/{}/bases'.format(reference.getId())
        length = reference.getLength()
        badRanges = [(-1, 0), (-1, -1), (length, 0), (0, length + 1)]
        for start, end in badRanges:
            args = protocol.ListReferenceBasesRequest()
            args.start, args.end = start, end
            response = self.app.get(path, query_string=args.toJsonDict())
            self.assertEqual(response.status_code, 416)

    def testListReferenceBasesPaging(self):
        referenceSet = self.backend.getReferenceSets()[0]
        reference = referenceSet.getReferences()[0]
        id_ = reference.getId()
        length = reference.getLength()
        completeSequence = reference.getBases(0, length)
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
                while response.nextPageToken is not None:
                    args = protocol.ListReferenceBasesRequest()
                    args.pageToken = response.nextPageToken
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
        for dataset in self.backend.getDatasets():
            for readGroupSet in dataset.getReadGroupSets():
                referenceSet = readGroupSet.getReferenceSet()
                for reference in referenceSet.getReferences():
                    for readGroup in readGroupSet.getReadGroups():
                        # search reads
                        request = protocol.SearchReadsRequest()
                        request.readGroupIds = [readGroup.getId()]
                        request.referenceId = reference.getId()
                        responseData = self.sendSearchRequest(
                            path, request, protocol.SearchReadsResponse)
                        alignments = responseData.alignments
                        self.assertGreater(len(alignments), 0)
                        for alignment in alignments:
                            # TODO more tests here: this is very weak.
                            self.assertEqual(
                                alignment.readGroupId, readGroup.getId())
