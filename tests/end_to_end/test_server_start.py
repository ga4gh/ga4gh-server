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

import tests.paths as paths


class TestServerStart(unittest.TestCase):

    def testServerStart(self):
        app = server.Ga4ghServerForTestingDataSource(paths.testDataRepo)
        try:
            app.start()
        finally:
            app.shutdown()
