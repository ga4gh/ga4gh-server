"""
Unit tests for frontend error conditions.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.frontend as frontend
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol
import ga4gh.avrotools as avrotools


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

    def _createInstance(self, requestClass):
        """
        Returns a valid instance of the specified class.
        """
        creator = avrotools.Creator(requestClass)
        instance = creator.getTypicalInstance()
        return instance

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
        error = protocol.GAException.fromJsonString(response.data)
        self.assertEqual(
            error.errorCode, exceptionClass.getErrorCode())
        self.assertGreater(len(error.message), 0)

    def assertRequestRaises(self, exceptionClass, url, request):
        """
        Verifies that the specified request returns a protocol exception
        corresponding to the specified exception class.
        """
        self.assertRawRequestRaises(
            exceptionClass, url, request.toJsonString())

    def testPageSize(self):
        for url, requestClass in self.endPointMap.items():
            for badType in ["", "1", "None", 0.0, 1e3]:
                request = self._createInstance(requestClass)
                request.pageSize = badType
                self.assertRequestRaises(
                    exceptions.RequestValidationFailureException, url, request)
            for badSize in [-100, -1, 0]:
                request = self._createInstance(requestClass)
                request.pageSize = badSize
                self.assertRequestRaises(
                    exceptions.BadPageSizeException, url, request)

    def testPageToken(self):
        for url, requestClass in self.endPointMap.items():
            for badType in [0, 0.0, 1e-3, {}, [], [None]]:
                request = self._createInstance(requestClass)
                request.pageToken = badType
                self.assertRequestRaises(
                    exceptions.RequestValidationFailureException, url, request)
