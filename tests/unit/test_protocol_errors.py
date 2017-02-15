"""
Unit tests for frontend error conditions.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.server.frontend as frontend
import ga4gh.server.exceptions as exceptions

import ga4gh.schemas.protocol as protocol


class TestFrontendErrors(unittest.TestCase):
    """
    Tests the frontend for various errors that can occur and verify
    that the correct exception was raised by the error code sent
    back.
    """
    @classmethod
    def setUpClass(cls):
        frontend.reset()
        frontend.configure(baseConfig="TestConfig")
        cls.app = frontend.app.test_client()

    @classmethod
    def tearDownClass(cls):
        cls.app = None

    def setUp(self):
        # TODO replace this with ALL post methods once the rest of the
        # end points have been implemented. This should also add an API
        # to protocol.py to simplify and document the process of getting
        # the correct API endpoints and classes. That is, we shouldn't
        # use protocol.postMethods directly, but instead call a function.
        supportedMethods = set([
            protocol.SearchCallSetsRequest,
            protocol.SearchVariantSetsRequest,
            protocol.SearchVariantsRequest,
        ])
        self.endPointMap = {}
        for endPoint, requestClass, responseClass in protocol.postMethods:
            if requestClass in supportedMethods:
                self.endPointMap[endPoint] = requestClass

    def assertRawRequestRaises(self, exceptionClass, url, requestString):
        """
        Verifies that the specified request string returns a protocol
        exception corresponding to the specified class when applied to
        all POST endpoints.
        """
        response = self.app.post(
            url, headers={'Content-type': 'application/json'},
            data=requestString)
        self.assertEqual(response.status_code, exceptionClass.httpStatus)
        error = protocol.fromJson(response.data, protocol.GAException)
        self.assertEqual(
            error.error_code, exceptionClass.getErrorCode())
        self.assertGreater(len(error.message), 0)

    def assertRequestRaises(self, exceptionClass, url, request):
        """
        Verifies that the specified request returns a protocol exception
        corresponding to the specified exception class.
        """
        self.assertRawRequestRaises(
            exceptionClass, url, protocol.toJson(request))

    def testPageSize(self):
        for url, requestClass in self.endPointMap.items():
            for badSize in [-100, -1]:
                request = requestClass()
                request.page_size = badSize
                self.assertRequestRaises(
                    exceptions.BadPageSizeException, url, request)

    @unittest.skip("Gets caught by the protocol buffer checkers")
    def testPageToken(self):
        for url, requestClass in self.endPointMap.items():
            for badType in [0, 0.0, 1e-3, {}, [], [None]]:
                request = requestClass()
                request.page_token = badType
                self.assertRequestRaises(
                    exceptions.RequestValidationFailureException, url, request)

    @unittest.skip("TODO: create invalid JSON to test validation")
    def testInvalidFields(self):
        for url, requestClass in self.endPointMap.items():
            request = self._createInvalidInstance(requestClass)
            self.assertRequestRaises(
                exceptions.RequestValidationFailureException, url, request)
