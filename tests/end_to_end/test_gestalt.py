"""
An end to end test which tests:
- client cmd line parsing
- client operation
- client logging
- server cmd line parsing
- server operation
- simulated variantSet backend
- server logging
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import server_test as server_test
import client as client


class TestGestalt(server_test.ServerTest):
    """
    An end-to-end test of the client and server
    """
    def testEndToEnd(self):
        self.client = client.ClientForTesting(self.server.getUrl())
        self.runVariantsRequest()
        self.assertLogsWritten()
        self.runReadsRequest()
        self.runReferencesRequest()
        self.runVariantSetsRequestDatasetTwo()
        self.client.cleanup()

    def assertLogsWritten(self):
        serverOutLines = self.server.getOutLines()
        serverErrLines = self.server.getErrLines()
        clientOutLines = self.client.getOutLines()
        clientErrLines = self.client.getErrLines()

        # nothing should be written to server stdout
        self.assertEqual(
            [], serverOutLines,
            "Server stdout log not empty")

        # server stderr should log at least one response success
        responseFound = False
        for line in serverErrLines:
            if ' 200 ' in line:
                responseFound = True
                break
        self.assertTrue(
            responseFound,
            "No successful server response logged to stderr")

        # client stdout should not be empty
        self.assertNotEqual(
            [], clientOutLines,
            "Client stdout log is empty")

        # num of client stdout should be twice the value of
        # SIMULATED_BACKEND_NUM_VARIANT_SETS
        # times the number of datasets (currently 2)
        expectedNumClientOutLines = 40
        self.assertEqual(len(clientOutLines), expectedNumClientOutLines)

        # client stderr should log at least one post
        requestFound = False
        for line in clientErrLines:
            if 'POST' in line:
                requestFound = True
                break
        self.assertTrue(
            requestFound,
            "No request logged from the client to stderr")

    def runVariantsRequest(self):
        self.runClientCmd(self.client, "variants-search -s0 -e2")

    def runReadsRequest(self):
        cmd = "reads-search --read_group_ids 'aReadGroupSet:one'"
        self.runClientCmd(self.client, cmd)

    def runReferencesRequest(self):
        reference_set_id = 'referenceSet0'
        reference_id = '{}:srs0'.format(reference_set_id)
        cmd = "referencesets-search"
        self.runClientCmd(self.client, cmd)
        cmd = "references-search"
        self.runClientCmd(self.client, cmd)
        cmd = "referencesets-get {}".format(reference_set_id)
        self.runClientCmd(self.client, cmd)
        cmd = "references-get {}".format(reference_id)
        self.runClientCmd(self.client, cmd)
        cmd = "references-list-bases {}".format(reference_id)
        self.runClientCmd(self.client, cmd)

    def runVariantSetsRequestDatasetTwo(self):
        dataset_id = "simulatedDataset2"
        cmd = "variantsets-search --dataset_ids {}".format(dataset_id)
        self.runClientCmd(self.client, cmd)
