"""
Tests for the client
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import json

import mock
import requests

import tests
import ga4gh.client as client
import ga4gh.protocol as protocol


class DummyRequest(protocol.ProtocolElement):

    def __init__(self):
        self.stringVal = "stringVal"
        self.intVal = 1
        self.arrayVal = [1, 2, 3]
        self.pageToken = None


class DummyResponse(object):

    def __init__(self, text=None):
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
    httpClient = client.HTTPClient(url, debugLevel, workarounds, key)
    return httpClient


class TestSearchMethodsCallRunRequest(unittest.TestCase):
    """
    Test that search methods call lower-level functionality correctly
    """
    def setUp(self):
        self.httpClient = makeHttpClient()
        self.request = DummyRequest()
        self.httpClient.runRequest = mock.Mock()

    def testSearchVariants(self):
        self.httpClient.searchVariants(self.request)
        self.httpClient.runRequest.assert_called_once_with(
            self.request, "variants/search",
            protocol.GASearchVariantsResponse, "variants")

    def testSearchVariantSets(self):
        self.httpClient.searchVariantSets(self.request)
        self.httpClient.runRequest.assert_called_once_with(
            self.request, "variantsets/search",
            protocol.GASearchVariantSetsResponse, "variantSets")

    def testSearchReferenceSets(self):
        self.httpClient.searchReferenceSets(self.request)
        self.httpClient.runRequest.assert_called_once_with(
            self.request, "referencesets/search",
            protocol.GASearchReferenceSetsResponse, "referenceSets")

    def testSearchReferences(self):
        self.httpClient.searchReferences(self.request)
        self.httpClient.runRequest.assert_called_once_with(
            self.request, "references/search",
            protocol.GASearchReferencesResponse, "references")

    def testSearchReadGroupSets(self):
        self.httpClient.searchReadGroupSets(self.request)
        self.httpClient.runRequest.assert_called_once_with(
            self.request, "readgroupsets/search",
            protocol.GASearchReadGroupSetsResponse, "readGroupSets")

    def testSearchCallSets(self):
        self.httpClient.searchCallSets(self.request)
        self.httpClient.runRequest.assert_called_once_with(
            self.request, "callsets/search",
            protocol.GASearchCallSetsResponse, "callSets")

    def testSearchReads(self):
        self.httpClient.searchReads(self.request)
        self.httpClient.runRequest.assert_called_once_with(
            self.request, "reads/search",
            protocol.GASearchReadsResponse, "alignments")


class TestRunRequest(unittest.TestCase):
    """
    Test the logic of runRequest
    """
    def testRunRequest(self):
        # setup
        mockPost = mock.Mock()
        with mock.patch('requests.post', mockPost):
            mockPost.side_effect = [DummyResponse(), DummyResponse('{}')]
            request = DummyRequest()

            # invoke SUT
            url = "referencesets/search"
            protocolClass = protocol.GASearchReferenceSetsResponse
            listAttr = "referenceSets"
            httpClient = makeHttpClient()
            result = [refSet for refSet in httpClient.runRequest(
                request, url, protocolClass, listAttr)]

            # verify results of invocation
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0].id, "refA")
            self.assertEqual(result[0].md5checksum, "abc")
            self.assertEqual(result[1].id, "refB")

            # verify requests.post called correctly
            url = "http://example.com/referencesets/search"
            jsonString = request.toJSONString()
            headers = {"Content-type": "application/json"}
            mockPost.assert_called_twice_with(
                url, jsonString, headers=headers, verify=False)
