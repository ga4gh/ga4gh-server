"""
Unit tests for the frontend code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import logging

import tests.paths as paths

import ga4gh.server.frontend as frontend


class TestAuth0(unittest.TestCase):
    """
    Tests the basic routing and HTTP handling for the Flask app.
    """
    exampleUrl = 'www.example.com'

    def sendGetRequest(self, path):
        """
        Sends a get request to the specified URL and returns the response.
        """
        headers = {
            'Origin': self.exampleUrl,
        }
        return self.app.get(path, headers=headers)

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
            baseConfig="Auth0Config", extraConfig=config)
        cls.app = frontend.app.test_client()
        # silence usually unhelpful CORS log
        logging.getLogger('ga4gh.frontend.cors').setLevel(logging.CRITICAL)

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
