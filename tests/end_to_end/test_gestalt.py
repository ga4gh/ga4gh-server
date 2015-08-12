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

import server_test
import client


class TestGestalt(server_test.ServerTest):
    """
    An end-to-end test of the client and server
    """
    def testEndToEnd(self):
        self.simulatedDatasetId = "c2ltdWxhdGVkRGF0YXNldDA="
        self.simulatedVariantSetId = "c2ltdWxhdGVkRGF0YXNldDA6c2ltVnMw"
        self.simulatedReadGroupId = "c2ltdWxhdGVkRGF0YXNldDA6c2ltUmdzMDpyZzA="
        self.simulatedReferenceSetId = "cmVmZXJlbmNlU2V0MA=="
        self.simulatedReferenceId = "cmVmZXJlbmNlU2V0MDpzcnMw"
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

        # number of variants to expect
        expectedNumClientOutLines = 2
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
        self.runClientCmd(
            self.client,
            "variants-search -s 0 -e 2 -V {}".format(
                self.simulatedVariantSetId))

    def runReadsRequest(self):
        cmd = (
            "reads-search --readGroupIds {} "
            "--referenceId chr1".format(
                self.simulatedReadGroupId))
        self.runClientCmd(self.client, cmd)

    def runReferencesRequest(self):
        referenceSetId = self.simulatedReferenceSetId
        referenceId = self.simulatedReferenceId
        cmd = "referencesets-search"
        self.runClientCmd(self.client, cmd)
        cmd = "references-search --referenceSetId={}".format(referenceSetId)
        self.runClientCmd(self.client, cmd)
        cmd = "referencesets-get {}".format(referenceSetId)
        self.runClientCmd(self.client, cmd)
        cmd = "references-get {}".format(referenceId)
        self.runClientCmd(self.client, cmd)
        cmd = "references-list-bases {}".format(referenceId)
        self.runClientCmd(self.client, cmd)

    def runVariantSetsRequestDatasetTwo(self):
        datasetId = self.simulatedDatasetId
        cmd = "variantsets-search --datasetId {}".format(datasetId)
        self.runClientCmd(self.client, cmd)
