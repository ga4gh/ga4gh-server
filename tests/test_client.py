"""
Tests for the client
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import json

import mock

import ga4gh.client as client
import ga4gh.protocol as protocol


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


def makeHttpClient():
    url = "http://example.com"
    debugLevel = 1
    workarounds = set()
    key = "KEY"
    httpClient = client.HttpClient(url, debugLevel, workarounds, key)
    return httpClient


class TestSearchMethodsCallRunRequest(unittest.TestCase):
    """
    Test that search methods call lower-level functionality correctly
    """
    def setUp(self):
        self.httpClient = makeHttpClient()
        self.protocolRequest = DummyRequest()
        self.httpClient.runSearchRequest = mock.Mock()
        self.httpClient.runListRequest = mock.Mock()
        self.httpClient.runGetRequest = mock.Mock()
        self._id = "SomeId"

    def testSearchVariants(self):
        self.httpClient.searchVariants(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "variants",
            protocol.GASearchVariantsResponse, "variants")

    def testSearchVariantSets(self):
        self.httpClient.searchVariantSets(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "variantsets",
            protocol.GASearchVariantSetsResponse, "variantSets")

    def testSearchReferenceSets(self):
        self.httpClient.searchReferenceSets(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "referencesets",
            protocol.GASearchReferenceSetsResponse, "referenceSets")

    def testSearchReferences(self):
        self.httpClient.searchReferences(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "references",
            protocol.GASearchReferencesResponse, "references")

    def testSearchReadGroupSets(self):
        self.httpClient.searchReadGroupSets(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "readgroupsets",
            protocol.GASearchReadGroupSetsResponse, "readGroupSets")

    def testSearchCallSets(self):
        self.httpClient.searchCallSets(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "callsets",
            protocol.GASearchCallSetsResponse, "callSets")

    def testSearchReads(self):
        self.httpClient.searchReads(self.protocolRequest)
        self.httpClient.runSearchRequest.assert_called_once_with(
            self.protocolRequest, "reads",
            protocol.GASearchReadsResponse, "alignments")

    def testGetReferenceSet(self):
        self.httpClient.getReferenceSet(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "referencesets", protocol.GAReferenceSet, self._id)

    def testGetReference(self):
        self.httpClient.getReference(self._id)
        self.httpClient.runGetRequest.assert_called_once_with(
            "references", protocol.GAReference, self._id)

    def testListReferenceBases(self):
        self.httpClient.listReferenceBases(self.protocolRequest, self._id)
        self.httpClient.runListRequest.assert_called_once_with(
            self.protocolRequest, "references/{id}/bases",
            protocol.GAListReferenceBasesResponse, self._id)


class TestRunRequest(unittest.TestCase):
    """
    Test the logic of the run*Request methods
    """
    def setUp(self):
        self.httpClient = makeHttpClient()

    def testRunSearchRequest(self):
        # setup
        mockPost = mock.Mock()
        with mock.patch('requests.request', mockPost):
            mockPost.side_effect = [DummyResponse(), DummyResponse('{}')]
            protocolRequest = DummyRequest()
            objectName = "referencesets"
            protocolResponseClass = protocol.GASearchReferenceSetsResponse
            listAttr = "referenceSets"

            # invoke SUT
            result = [refSet for refSet in self.httpClient.runSearchRequest(
                protocolRequest, objectName,
                protocolResponseClass, listAttr)]

            # verify results of invocation
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0].id, "refA")
            self.assertEqual(result[0].md5checksum, "abc")
            self.assertEqual(result[1].id, "refB")

            # verify requests.post called correctly
            url = "http://example.com/referencesets/search"
            jsonString = protocolRequest.toJsonString()
            headers = {"Content-type": "application/json"}
            httpMethod = 'POST'
            mockPost.assert_called_twice_with(
                httpMethod, url, jsonString, headers=headers, verify=False)

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
            protocolResponseClass = protocol.GAReference
            id_ = 'anId'

            # invoke SUT
            result = self.httpClient.runGetRequest(
                objectName, protocolResponseClass, id_)

            # verify results of invocation
            self.assertEqual(result.id, "gaid")
            self.assertEqual(result.md5checksum, "def")

            # verify requests.get called correctly
            url = "http://example.com/reference/anId"
            params = {}
            httpMethod = 'GET'
            headers = {}
            data = None
            mockGet.assert_called_once_with(
                httpMethod, url, params=params, data=data, headers=headers)

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
            protocolRequest = protocol.GAListReferenceBasesRequest()
            protocolRequest.start = 1
            protocolRequest.end = 5
            url = "references/{id}/bases"
            protocolResponseClass = protocol.GAListReferenceBasesResponse
            id_ = 'myId'

            # invoke SUT
            result = [base for base in self.httpClient.runListRequest(
                protocolRequest, url, protocolResponseClass, id_)]

            # verify results of invocation
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0].offset, 123)
            self.assertEqual(result[0].sequence, "sequence")

            # verify requests.get called correctly
            url = "http://example.com/references/myId/bases"
            params = {"start": 1, "end": 5}
            httpMethod = 'GET'
            mockGet.assert_called_twice_with(httpMethod, url, params=params)
