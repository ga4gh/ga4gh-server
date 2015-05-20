"""
Tests the ability of the server to fetch remote data files
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import client as client
import server as server
import server_test as server_test


class RemoteServerTestReads(server_test.RemoteServerTest):
    """
    Configures a RemoteServerTest subclass with read data
    """
    def getRemoteDataDir(self):
        return "tests/data/dataset1/reads/wgBam/"

    def getServer(self):
        dataDir = "tests/end_to_end/remoteTestDataReads"
        return server.Ga4ghServerForTestingDataSource(dataDir)


class TestRemoteFetchReads(RemoteServerTestReads):
    """
    Tests fetching read data from files on a remote server
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


class RemoteServerTestVariants(server_test.RemoteServerTest):
    """
    Configures a RemoteServerTest subclass with variant data
    """
    def getRemoteDataDir(self):
        return "tests/data/dataset1/variants/1kgPhase1"

    def getServer(self):
        dataDir = "tests/end_to_end/remoteTestDataVariants"
        return server.Ga4ghServerForTestingDataSource(dataDir)


class TestRemoteFetchVariants(RemoteServerTestVariants):
    """
    Tests fetching variant data from files on a remote server
    """
    def testVcfFetch(self):
        self.client = client.ClientForTesting(self.server.getUrl())
        self.runClientCmd(self.client, "variants-search")
        self._assertLogsWritten()
        self.client.cleanup()

    def _assertLogsWritten(self):
        # client should have gotten some variants back
        self.assertEqual(len(self.client.getOutLines()), 100)
