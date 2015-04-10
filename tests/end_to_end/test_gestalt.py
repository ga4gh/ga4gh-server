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
import unittest

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
DEBUG = True
"""

    def setUp(self):
        self.clientOutFile = None
        self.clientErrFile = None

    @unittest.skipIf(protocol.version.startswith("0.6"), "")
    def testEndToEnd(self):
        self.createLogFiles()
        self.runVariantSetRequest()
        self.assertLogsWritten()
        self.runReadsRequest()
        self.removeLogFiles()

    def getClientOutLines(self):
        self.clientOutFile.flush()
        self.clientOutFile.seek(0)
        clientOutLines = self.clientOutFile.readlines()
        return clientOutLines

    def getClientErrLines(self):
        self.clientErrFile.flush()
        self.clientErrFile.seek(0)
        clientErrLines = self.clientErrFile.readlines()
        return clientErrLines

    def createLogFiles(self):
        self.clientOutFile = tempfile.TemporaryFile()
        self.clientErrFile = tempfile.TemporaryFile()

    def removeLogFiles(self):
        self.clientOutFile.close()
        self.clientErrFile.close()

    def assertLogsWritten(self):
        serverOutLines = self.getServerOutLines()
        serverErrLines = self.getServerErrLines()
        clientOutLines = self.getClientOutLines()
        clientErrLines = self.getClientErrLines()

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
        expectedNumClientOutLines = 20
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

    def runVariantSetRequest(self):
        self._runClientCmdLine("variants-search -s0 -e2")

    def runReadsRequest(self):
        cmd = "reads-search --readGroupIds 'aReadGroupSet:one'"
        self._runClientCmdLine(cmd)

    def _runClientCmdLine(self, cmd):
        clientCmdLine = "python client_dev.py -v -O {} {}/v{}".format(
            cmd, self.serverUrl, protocol.version)
        splits = shlex.split(clientCmdLine)
        try:
            subprocess.check_call(
                splits,
                stdout=self.clientOutFile,
                stderr=self.clientErrFile)
        except subprocess.CalledProcessError as error:
            serverOutLines = self.getServerOutLines()
            serverErrLines = self.getServerErrLines()
            clientOutLines = self.getClientOutLines()
            clientErrLines = self.getClientErrLines()
            print('\n')
            print('*** Server CONFIG ***')
            print(self.config)
            print('*** Server CMD ***')
            print(self.serverCmdLine)
            print('*** Client CMD ***')
            print(clientCmdLine)
            print('*** Server STDOUT ***')
            print(''.join(serverOutLines))
            print('*** Server STDERR ***')
            print(''.join(serverErrLines))
            print('*** Client STDOUT ***')
            print(''.join(clientOutLines))
            print('*** Client STDERR ***')
            print(''.join(clientErrLines))
            raise error
