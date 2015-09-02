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
        self.protocolRequest = DummyRequest()
        self.httpClient.runSearchRequest = mock.Mock()
        self.httpClient.runListRequest = mock.Mock()
        self.httpClient.runGetRequest = mock.Mock()
        self._id = "SomeId"

    def testSearchVariants(self):
        self.httpClient.searchVariants(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "variants",
            protocol.SearchVariantsResponse)

    def testSearchVariantSets(self):
        self.httpClient.searchVariantSets(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "variantsets",
            protocol.SearchVariantSetsResponse)

    def testSearchReferenceSets(self):
        self.httpClient.searchReferenceSets(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "referencesets",
            protocol.SearchReferenceSetsResponse)

    def testSearchReferences(self):
        self.httpClient.searchReferences(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "references",
            protocol.SearchReferencesResponse)

    def testSearchReadGroupSets(self):
        self.httpClient.searchReadGroupSets(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "readgroupsets",
            protocol.SearchReadGroupSetsResponse)

    def testSearchCallSets(self):
        self.httpClient.searchCallSets(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "callsets",
            protocol.SearchCallSetsResponse)

    def testSearchReads(self):
        self.httpClient.searchReads(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "reads",
            protocol.SearchReadsResponse)

    def testSearchDatasets(self):
        self.httpClient.searchDatasets(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "datasets",
            protocol.SearchDatasetsResponse)

    def testGetReferenceSet(self):
        self.httpClient.getReferenceSet(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "referencesets", protocol.ReferenceSet, self._id)

    def testGetVariantSet(self):
        self.httpClient.getVariantSet(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "variantsets", protocol.VariantSet, self._id)

    def testGetReference(self):
        self.httpClient.getReference(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "references", protocol.Reference, self._id)

    def testGetReadGroupSets(self):
        self.httpClient.getReadGroupSet(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "readgroupsets", protocol.ReadGroupSet, self._id)

    def testGetReadGroup(self):
        self.httpClient.getReadGroup(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "readgroups", protocol.ReadGroup, self._id)

    def testGetCallsets(self):
        self.httpClient.getCallset(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "callsets", protocol.CallSet, self._id)

    def testGetDatasets(self):
        self.httpClient.getDataset(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "datasets", protocol.Dataset, self._id)

    def testGetVariant(self):
        self.httpClient.getVariant(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "variants", protocol.Variant, self._id)

    def testListReferenceBases(self):
        self.httpClient.listReferenceBases(self.protocolRequest, self._id)
        self.httpClient.runListRequest.assert_called_once_with(
            self.protocolRequest, "references/{id}/bases",
            protocol.ListReferenceBasesResponse, self._id)


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

    def testRunListRequest(self):
        # setup
        mockGet = mock.Mock()
        with mock.patch('requests.request', mockGet):
            text = {
                "offset": 123,
                "sequence": "sequence",
                "nextPageToken": "pageTok",
            }
            mockGet.side_effect = [
                DummyResponse(json.dumps(text)), DummyResponse('{}')]
            protocolRequest = protocol.ListReferenceBasesRequest()
            protocolRequest.start = 1
            protocolRequest.end = 5
            url = "references/{id}/bases"
            protocolResponseClass = protocol.ListReferenceBasesResponse
            id_ = 'myId'

            # invoke SUT
            result = [base for base in self.httpClient.runListRequest(
                protocolRequest, url, protocolResponseClass, id_)]

            # verify results of invocation
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0].offset, 123)
            self.assertEqual(result[0].sequence, "sequence")

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

            # assert second call correct
            params['pageToken'] = 'pageTok'
            secondCall = mockGet.call_args_list[1]
            self.assertRequestsCall(
                secondCall, httpMethod, url, headers, data, params, False)

    def assertRequestsCall(
            self, call, httpMethod, url,
            headers, data, params, verify):
        self.assertEqual(call[0], (httpMethod, url))
        self.assertEqual(call[1]['headers'], headers)
        self.assertEqual(call[1]['data'], data)
        self.assertEqual(call[1]['params'], params)
        self.assertEqual(call[1]['verify'], False)
