"""
Tests that the server starts up correctly with a file system backend.

Necessary because all other tests that start up the server do so
with non-file backends (e.g. simulated backends).
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import server as server


class TestServerStart(unittest.TestCase):

    @unittest.skip("Skip test until test data is converted to DB repo")
    def testServerStart(self):
        dataDir = "tests/data"
        app = server.Ga4ghServerForTestingDataSource(dataDir)
        try:
            app.start()
        finally:
            app.shutdown()
