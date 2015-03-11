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

import tempfile
import shlex
import subprocess

import ga4gh.protocol as protocol
import server as server


class TestGestalt(server.ServerTestConfigFile):
    """
    An end-to-end test of the client and server
    """
    def setConfig(self):
        self.config = """
SIMULATED_BACKEND_NUM_VARIANT_SETS = 10
SIMULATED_BACKEND_VARIANT_DENSITY = 1
DATA_SOURCE = "__SIMULATED__"
"""

    def setUp(self):
        self.client = None
        self.clientOutFile = None
        self.clientErrFile = None

    def testEndToEnd(self):
        self.createLogFiles()
        self.runRequest()
        self.assertLogsWritten()
        self.removeLogFiles()

    def createLogFiles(self):
        self.clientOutFile = tempfile.TemporaryFile()
        self.clientErrFile = tempfile.TemporaryFile()

    def removeLogFiles(self):
        self.clientOutFile.close()
        self.clientErrFile.close()

    def assertLogsWritten(self):
        self.clientOutFile.seek(0)
        self.clientErrFile.seek(0)
        serverOutLines = self.getServerOutLines()
        serverErrLines = self.getServerErrLines()
        clientOutLines = self.clientOutFile.readlines()
        clientErrLines = self.clientErrFile.readlines()

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
        self.assertEqual(len(clientOutLines), 20)

        # client stderr should log at least one post
        requestFound = False
        for line in clientErrLines:
            if 'POST' in line:
                requestFound = True
                break
        self.assertTrue(
            requestFound,
            "No request logged from the client to stderr")

    def runRequest(self):
        clientCmdLine = """python client_dev.py -v -O variants-search
            -s0 -e2 {}/v{}""".format(
            self.serverUrl, protocol.version)
        splits = shlex.split(clientCmdLine)
        self.client = subprocess.check_call(
            splits, stdout=self.clientOutFile, stderr=self.clientErrFile)
