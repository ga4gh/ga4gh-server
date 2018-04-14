"""
Tests related to the server start script
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.server.cli.server as server
import ga4gh.server.frontend as frontend


class TestExceptionHandler(unittest.TestCase):
    """
    Test that the server script functions behave in expected ways.
    """
    def testGetServerParser(self):
        self.assertIsNotNone(server.getServerParser(),
                             "The server parser should be returned so "
                             "that we can create docs for sphinx")

    def test_number_of_workers(self):
        self.assertTrue(type(server.number_of_workers()) == int,
                        "The number of workers function should return an "
                        "integer.")

    def testStandaloneApplicationInstance(self):
        app = server.StandaloneApplication(frontend.app)
        self.assertIsNotNone(app.run, "Ensures the class instantiates from "
                                      "our WSGI app properly. The run "
                                      "function will spawn many processes.")
