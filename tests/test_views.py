"""
Tests for the frontend.
"""

import unittest

import ga4gh.server as server
import ga4gh.cli

class TestFrontend(unittest.TestCase):

    def setUp(self):
        self.app = server.app.test_client()

    def testServer(self):
        self.assertTrue('404' in self.app.get('/').data)

    # TODO: Fill in actual test cases here...

if __name__ == '__main__':
    unittest.main()
