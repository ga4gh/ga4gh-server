"""
Unit tests for the frontend code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import shutil
import unittest
import logging

import ga4gh.server.frontend as frontend
import ga4gh.server.exceptions as exceptions
import ga4gh.server.auth as auth

import ga4gh.schemas.protocol as protocol


class TestAuth0(unittest.TestCase):
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
            "CACHE_DIRECTORY": "/tmp/ga4gh-test"
        }
        frontend.reset()
        frontend.configure(
            baseConfig="TestAuth0Config", extraConfig=config)
        cls.app = frontend.app
        # silence usually unhelpful CORS log
        logging.getLogger('ga4gh.frontend.cors').setLevel(logging.CRITICAL)
        cls.backend = frontend.app.backend
        cls.client = frontend.app.test_client()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree('/tmp/ga4gh-test', True)
        frontend.reset()

    def sendPostRequest(self, path, request, extraHeaders=None):
        """
        Sends the specified GA request object and returns the response.
        """
        headers = {
            'Content-type': 'application/json',
            'Origin': self.exampleUrl,
        }
        if extraHeaders:
            headers.update(extraHeaders)
        return self.client.post(
            path, headers=headers, data=protocol.toJson(request))

    def sendGetRequest(self, path):
        """
        Sends a get request to the specified URL and returns the response.
        """
        headers = {
            'Origin': self.exampleUrl,
        }
        return self.client.get(path, headers=headers)

    def testCallback(self):
        response = self.sendGetRequest("callback")
        self.assertEqual(
            response.status_code,
            401, "Ensure that when the callback is hit without a code"
                 "it will return 401 but got {}".format(response.status_code))
        response = self.sendGetRequest("callback?code=abc")
        self.assertEqual(
            response.status_code,
            401, "Ensure that when the callback is hit without a code"
                 "it will return 401 but got {}".format(response.status_code))

    def testLogin(self):
        response = self.sendGetRequest("login")
        self.assertEqual(
            response.status_code,
            200, "Ensure that when Auth0 is turned on the login page"
                 "returns 200 {}".format(response.status_code))

    def testBadBearer(self):
        """
        Tests to see if a malformed bearer token fails in expected ways
        """
        response = self.client.get('/')
        self.assertEqual(response.status_code, 401)
        protectedPath = "datasets/search"
        request = protocol.SearchDatasetsRequest()
        headers = {"Authorization": ""}
        response = self.sendPostRequest(protectedPath, request, headers)
        self.assertEquals(response.status_code, 401, "No bearer should fail"
                                                     "with 401")
        headers = {"Authorization": "Bearer"}
        response = self.sendPostRequest(protectedPath, request, headers)
        self.assertEquals(response.status_code, 401, "")

    def testProtected(self):
        protectedPath = "datasets/search"
        request = protocol.SearchDatasetsRequest()
        response = self.sendPostRequest(protectedPath, request)
        self.assertEquals(
            response.status_code,
            401, "If Auth0 is enabled this endpoint "
                 "is not accessible without auth.")

    def testDecodeExceptions(self):
        """
        Unit tests the header parsing functions.
        """
        with self.assertRaises(exceptions.NotAuthorizedException):
            auth._has_header('')
        with self.assertRaises(exceptions.NotAuthorizedException):
            auth._has_bearer('empty')
        with self.assertRaises(exceptions.NotAuthorizedException):
            auth._has_token('Bearer')
        with self.assertRaises(exceptions.NotAuthorizedException):
            auth._has_token('Bearer123')
        with self.assertRaises(exceptions.NotAuthorizedException):
            auth._well_formed('Bearer 123 456')
        with self.assertRaises(exceptions.NotAuthorizedException):
            client_id = self.app.config.get("AUTH0_CLIENT_ID")
            client_secret = self.app.config.get("AUTH0_CLIENT_SECRET")
            token, payload = auth._decode_header(
                'Bearer 123', client_id, client_secret)

    def testRenderLogin(self):
        """
        Tests that the login template renders without failing
        """
        return auth.render_login(
            app=self.app,
            scopes=self.app.config.get('AUTH0_SCOPES'),
            redirect_uri=self.app.config.get('AUTH0_CALLBACK_URL'),
            domain=self.app.config.get('AUTH0_HOST'),
            client_id=self.app.config.get('AUTH0_CLIENT_ID'))

    def testAuthorizeEmail(self):
        """
        Tests that the email is set to authorized in the cache.
        """
        email = "test@test.it"
        auth.authorize_email(email, self.app.cache)
        entry = self.app.cache.get(email)
        self.assertTrue(entry['authorized'])
        self.assertTrue(auth.is_authorized(self.app.cache, email))

    def testIsActive(self):
        """
        Tests that is active throws an exception when the token is not found.
        """
        token = "123"
        email = "test@test.com"
        with self.assertRaises(exceptions.NotAuthenticatedException):
            auth.is_active(self.app.cache, token)
        self.app.cache.set("123", {"email": email})
        # Even though the token is set, in the cache, it
        self.assertEquals(
            email, auth.is_active(self.app.cache, token)['email'])

    def testCallbackMaker(self):
        """
        Tests that the callback maker returns a function and runs without
         failure.
        """
        callback = auth.callback_maker(
            cache=self.app.cache,
            domain=self.app.config.get('AUTH0_HOST'),
            client_id=self.app.config.get('AUTH0_CLIENT_ID'),
            client_secret=self.app.config.get('AUTH0_CLIENT_SECRET'),
            redirect_uri=self.app.config.get('AUTH0_CALLBACK_URL'))
        self.assertTrue(callable(callback))
