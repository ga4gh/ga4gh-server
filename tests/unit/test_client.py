"""
Tests for the client
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import json

import mock

import ga4gh.protocol as protocol
import tests.utils as utils


class DummyRequest(protocol.ProtocolElement):

    __slots__ = ["stringVal", "intVal", "arrayVal", "pageToken"]

    def __init__(self):
        self.stringVal = "stringVal"
        self.intVal = 1
        self.arrayVal = [1, 2, 3]
        self.pageToken = None

    def __eq__(self, other):
        for field in self.__slots__:
            if getattr(self, field) != getattr(other, field):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)


class DummyResponse(object):

    def __init__(self, text=None):
        self.status_code = 200
        if text is None:
            self.text = self._getText()
        else:
            self.text = text

    def _getText(self):
        txt = {
            "nextPageToken": "xyz",
            "referenceSets": [
                {"id": "refA", "md5checksum": "abc"},
                {"id": "refB"}
            ]
        }
        return json.dumps(txt)

    def raise_for_status(self):
        pass


class TestSearchMethodsCallRunRequest(unittest.TestCase):
    """
    Test that search methods call lower-level functionality correctly
    """
    def setUp(self):
        self.httpClient = utils.makeHttpClient()
        self.httpClient.runSearchRequest = mock.Mock()
        self.httpClient.runGetRequest = mock.Mock()
        self.objectId = "SomeId"
        self.objectName = "objectName"
        self.datasetId = "datasetId"
        self.variantSetId = "variantSetId"
        self.referenceSetId = "referenceSetId"
        self.referenceId = "referenceId"
        self.readGroupIds = ["readGroupId"]
        self.referenceName = "referenceName"
        self.start = 100
        self.end = 101
        self.referenceName = "referenceName"
        self.callSetIds = ["id1", "id2"]
        self.pageSize = 1000
        self.httpClient.setPageSize(self.pageSize)
        self.assemblyId = "assemblyId"
        self.accession = "accession"
        self.md5checksum = "md5checksum"

    def testSetPageSize(self):
        httpClient = utils.makeHttpClient()
        # pageSize is None by default
        self.assertIsNone(httpClient.getPageSize())
        for pageSize in [1, 10, 100]:
            httpClient.setPageSize(pageSize)
            self.assertEqual(httpClient.getPageSize(), pageSize)

    def testSearchVariants(self):
        request = protocol.SearchVariantsRequest()
        request.referenceName = self.referenceName
        request.start = self.start
        request.end = self.end
        request.variantSetId = self.variantSetId
        request.callSetIds = self.callSetIds
        request.pageSize = self.pageSize
        self.httpClient.searchVariants(
            self.variantSetId, start=self.start, end=self.end,
            referenceName=self.referenceName, callSetIds=self.callSetIds)
        self.httpClient.runSearchRequest.assert_called_once_with(
            request, "variants",
            protocol.SearchVariantsResponse)

    def testSearchDatasets(self):
        request = protocol.SearchDatasetsRequest()
        request.pageSize = self.pageSize
        self.httpClient.searchDatasets()
        self.httpClient.runSearchRequest.assert_called_once_with(
            request, "datasets", protocol.SearchDatasetsResponse)

    def testSearchVariantSets(self):
        request = protocol.SearchVariantSetsRequest()
        request.datasetId = self.datasetId
        request.pageSize = self.pageSize
        self.httpClient.searchVariantSets(self.datasetId)
        self.httpClient.runSearchRequest.assert_called_once_with(
            request, "variantsets", protocol.SearchVariantSetsResponse)

    def testSearchReferenceSets(self):
        request = protocol.SearchReferenceSetsRequest()
        request.pageSize = self.pageSize
        request.accession = self.accession
        request.md5checksum = self.md5checksum
        request.assemblyId = self.assemblyId
        self.httpClient.searchReferenceSets(
            accession=self.accession, md5checksum=self.md5checksum,
            assemblyId=self.assemblyId)
        self.httpClient.runSearchRequest.assert_called_once_with(
            request, "referencesets", protocol.SearchReferenceSetsResponse)

    def testSearchReferences(self):
        request = protocol.SearchReferencesRequest()
        request.referenceSetId = self.referenceSetId
        request.pageSize = self.pageSize
        request.accession = self.accession
        request.md5checksum = self.md5checksum
        self.httpClient.searchReferences(
            self.referenceSetId, accession=self.accession,
            md5checksum=self.md5checksum)
        self.httpClient.runSearchRequest.assert_called_once_with(
            request, "references", protocol.SearchReferencesResponse)

    def testSearchReadGroupSets(self):
        request = protocol.SearchReadGroupSetsRequest()
        request.datasetId = self.datasetId
        request.name = self.objectName
        request.pageSize = self.pageSize
        self.httpClient.searchReadGroupSets(
            self.datasetId, name=self.objectName)
        self.httpClient.runSearchRequest.assert_called_once_with(
            request, "readgroupsets", protocol.SearchReadGroupSetsResponse)

    def testSearchCallSets(self):
        request = protocol.SearchCallSetsRequest()
        request.variantSetId = self.variantSetId
        request.name = self.objectName
        request.pageSize = self.pageSize
        self.httpClient.searchCallSets(
            self.variantSetId, name=self.objectName)
        self.httpClient.runSearchRequest.assert_called_once_with(
            request, "callsets", protocol.SearchCallSetsResponse)

    def testSearchReads(self):
        request = protocol.SearchReadsRequest()
        request.readGroupIds = self.readGroupIds
        request.referenceId = self.referenceId
        request.start = self.start
        request.end = self.end
        request.pageSize = self.pageSize
        self.httpClient.searchReads(
            self.readGroupIds, referenceId=self.referenceId,
            start=self.start, end=self.end)
        self.httpClient.runSearchRequest.assert_called_once_with(
            request, "reads", protocol.SearchReadsResponse)

    def testGetReferenceSet(self):
        self.httpClient.getReferenceSet(self.objectId)
        self.httpClient.runGetRequest.assert_called_once_with(
            "referencesets", protocol.ReferenceSet, self.objectId)

    def testGetVariantSet(self):
        self.httpClient.getVariantSet(self.objectId)
        self.httpClient.runGetRequest.assert_called_once_with(
            "variantsets", protocol.VariantSet, self.objectId)

    def testGetReference(self):
        self.httpClient.getReference(self.objectId)
        self.httpClient.runGetRequest.assert_called_once_with(
            "references", protocol.Reference, self.objectId)

    def testGetReadGroupSets(self):
        self.httpClient.getReadGroupSet(self.objectId)
        self.httpClient.runGetRequest.assert_called_once_with(
            "readgroupsets", protocol.ReadGroupSet, self.objectId)

    def testGetReadGroup(self):
        self.httpClient.getReadGroup(self.objectId)
        self.httpClient.runGetRequest.assert_called_once_with(
            "readgroups", protocol.ReadGroup, self.objectId)

    def testGetCallsets(self):
        self.httpClient.getCallset(self.objectId)
        self.httpClient.runGetRequest.assert_called_once_with(
            "callsets", protocol.CallSet, self.objectId)

    def testGetDatasets(self):
        self.httpClient.getDataset(self.objectId)
        self.httpClient.runGetRequest.assert_called_once_with(
            "datasets", protocol.Dataset, self.objectId)

    def testGetVariant(self):
        self.httpClient.getVariant(self.objectId)
        self.httpClient.runGetRequest.assert_called_once_with(
            "variants", protocol.Variant, self.objectId)


class TestRunRequest(unittest.TestCase):
    """
    Test the logic of the run*Request methods
    """
    def setUp(self):
        self.httpClient = utils.makeHttpClient()

    def testRunSearchRequest(self):
        # setup
        mockPost = mock.Mock()
        with mock.patch('requests.request', mockPost):
            mockPost.side_effect = [DummyResponse(), DummyResponse('{}')]
            protocolRequest = DummyRequest()
            objectName = "referencesets"
            protocolResponseClass = protocol.SearchReferenceSetsResponse

            # invoke SUT
            result = [refSet for refSet in self.httpClient.runSearchRequest(
                protocolRequest, objectName, protocolResponseClass)]

            # verify results of invocation
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0].id, "refA")
            self.assertEqual(result[0].md5checksum, "abc")
            self.assertEqual(result[1].id, "refB")

            # verify requests.post called correctly
            httpMethod = 'POST'
            url = "http://example.com/referencesets/search"
            data = protocolRequest.toJsonString()
            headers = {"Content-type": "application/json"}
            params = {u'key': u'KEY'}
            self.assertEqual(len(mockPost.call_args_list), 2)

            # assert first call correct
            firstCall = mockPost.call_args_list[0]
            self.assertRequestsCall(
                firstCall, httpMethod, url, headers, data, params, False)

            # assert second call correct
            protocolRequest.pageToken = "xyz"
            data = protocolRequest.toJsonString()
            secondCall = mockPost.call_args_list[1]
            self.assertRequestsCall(
                secondCall, httpMethod, url, headers, data, params, False)

    def testRunGetRequest(self):
        # setup
        mockGet = mock.Mock()
        with mock.patch('requests.request', mockGet):
            text = {
                "id": "gaid",
                "md5checksum": "def",
            }
            mockGet.side_effect = [DummyResponse(json.dumps(text))]
            objectName = "reference"
            protocolResponseClass = protocol.Reference
            id_ = 'anId'

            # invoke SUT
            result = self.httpClient.runGetRequest(
                objectName, protocolResponseClass, id_)

            # verify results of invocation
            self.assertEqual(result.id, "gaid")
            self.assertEqual(result.md5checksum, "def")

            # verify requests.get called correctly
            url = "http://example.com/reference/anId"
            params = {'key': 'KEY'}
            httpMethod = 'GET'
            headers = {}
            data = None
            mockGet.assert_called_once_with(
                httpMethod, url, params=params, headers=headers, data=data,
                verify=False)

    def testRunListReferenceBases(self):
        # setup
        mockGet = mock.Mock()
        with mock.patch('requests.request', mockGet):
            text = {
                "offset": 123,
                "sequence": "sequence",
                "nextPageToken": "pageTok",
            }
            mockGet.side_effect = [DummyResponse(json.dumps(text))]
            id_ = 'myId'

            # invoke SUT
            result = [chunk for chunk in self.httpClient.listReferenceBases(
                id_, start=1, end=5)]

            # verify results of invocation
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0], "sequence")

            # verify requests.get called correctly
            httpMethod = 'GET'
            url = "http://example.com/references/myId/bases"
            params = {
                'start': 1, 'end': 5, 'key': 'KEY', 'pageToken': None}
            headers = {}
            data = None
            self.assertEqual(len(mockGet.call_args_list), 2)

            # assert first call correct
            firstCall = mockGet.call_args_list[0]
            self.assertRequestsCall(
                firstCall, httpMethod, url, headers, data, params, False)

    def assertRequestsCall(
            self, call, httpMethod, url,
            headers, data, params, verify):
        self.assertEqual(call[0], (httpMethod, url))
        self.assertEqual(call[1]['headers'], headers)
        self.assertEqual(call[1]['data'], data)
        self.assertEqual(call[1]['params'], params)
        self.assertEqual(call[1]['verify'], False)
