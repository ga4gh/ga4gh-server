"""
Tests the ability of the server to fetch remote data files
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import client as client
import server as server
import server_test as server_test


baseDir = "tests/end_to_end/"
dataDir = os.path.join(baseDir, "remoteTestData")
remoteDataDir = "tests/data/reads/wgBam/"


class Ga4ghServerForTestingDataSource(server.Ga4ghServerForTesting):
    """
    A test server that reads data from a data source
    """
    def getConfig(self):
        config = """
DATA_SOURCE = "{}"
DEBUG = True""".format(dataDir)
        return config


class RemoteServerTest(server_test.ServerTest):
    """
    A test that uses the data source server
    """
    def otherSetup(self):
        self.remoteServer = server.RemoteServerForTesting(remoteDataDir)
        self.remoteServer.start()

    def otherTeardown(self):
        self.remoteServer.shutdown()

    def otherExceptionHandling(self):
        self.remoteServer.printDebugInfo()

    def getServer(self):
        return Ga4ghServerForTestingDataSource()


class TestRemoteFetch(RemoteServerTest):
    """
    Tests fetching genomics data from files on a remote server
    """
    def testBamFetch(self):
        self.client = client.ClientForTesting(self.server.getUrl())
        self.readGroupIds = \
            'remoteTest:wgEncodeUwRepliSeqBg02esG1bAlnRep1_sample'
        self.runClientCmd(
            self.client,
            "reads-search --readGroupIds '{}'".format(
                self.readGroupIds))
        self._assertLogsWritten()
        self.client.cleanup()

    def _assertLogsWritten(self):
        # client should have gotten some reads back
        self.assertEqual(len(self.client.getOutLines()), 10)
        # hard to know what else to test here since the
        # remote server logs are out of our control
