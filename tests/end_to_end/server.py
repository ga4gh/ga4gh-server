"""
Assists in executing tests with a running server on localhost
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import logging
import requests
import shlex
import subprocess
import tempfile
import unittest

import tests.utils as utils


class ServerTest(unittest.TestCase):
    """
    Manages server bootup and shutdown for a test
    """
    def run(self, *args, **kwargs):
        self.setupServer()
        try:
            self.startServer()
            super(ServerTest, self).run(*args, **kwargs)
        finally:
            self.shutDownServer()

    def setupServer(self):
        # suppress requests package log messages
        logging.getLogger("requests").setLevel(logging.CRITICAL)

        # ensure another server isn't running, bail if it is
        self.port = 8001
        self.serverUrl = "http://localhost:{}".format(self.port)
        self.assertFalse(
            self._isServerRunning(),
            "Another server is running")

        # init test vars
        self.server = None
        self.serverCmdLine = None
        self.serverOutFile = None
        self.serverErrFile = None

    def startServer(self):
        self.serverOutFile = tempfile.TemporaryFile()
        self.serverErrFile = tempfile.TemporaryFile()
        splits = shlex.split(self.getServerCmdLine())
        self.server = subprocess.Popen(
            splits, stdout=self.serverOutFile,
            stderr=self.serverErrFile)
        self._waitForServerStartup()

    @utils.Timeout()
    @utils.Repeat()
    def _waitForServerStartup(self):
        self.server.poll()
        if self.server.returncode is not None:
            self._waitForServerErrLines()
            message = "Server process unexpectedly died; stderr: {0}"
            failMessage = message.format(self.getServerErrLines())
            self.fail(failMessage)
        return not self._isServerRunning()

    @utils.Timeout()
    @utils.Repeat()
    def _waitForServerErrLines(self):
        # not sure why there's some delay in getting the server
        # process' stderr...
        return self.getServerErrLines() == []

    def shutDownServer(self):
        if self._isServerRunning():
            self.server.kill()
        if self.server is not None:
            self.server.wait()
            self._assertServerShutdown()
        self.serverOutFile.close()
        self.serverErrFile.close()

    def _assertServerShutdown(self):
        shutdownString = "Server did not shut down correctly"
        self.assertIsNotNone(self.server.returncode, shutdownString)
        self.assertFalse(self._isServerRunning(), shutdownString)

    def _isServerRunning(self):
        try:
            response = self._pingServer()
            self.assertEqual(response.status_code, 200)
            return True
        except requests.ConnectionError:
            return False

    def _pingServer(self):
        response = requests.get(self.serverUrl)
        return response

    def getServerOutLines(self):
        self.serverOutFile.flush()
        self.serverOutFile.seek(0)
        serverOutLines = self.serverOutFile.readlines()
        return serverOutLines

    def getServerErrLines(self):
        self.serverErrFile.flush()
        self.serverErrFile.seek(0)
        serverErrLines = self.serverErrFile.readlines()
        return serverErrLines

    def getServerCmdLine(self):
        self.serverCmdLine = """
python server_dev.py
--dont-use-reloader
--config TestConfig
--port {} """.format(self.port)
        return self.serverCmdLine


class ServerTestConfigFile(ServerTest):
    """
    Launch a server test, but with a custom configuration file
    """
    def run(self, *args, **kwargs):
        self.config = None
        self.configFile = None
        self.configFilePath = None
        super(ServerTestConfigFile, self).run(*args, **kwargs)

    def setConfig(self):
        raise NotImplementedError("ServerTestConfigFile must be subclassed")

    def getServerCmdLine(self):
        self.setConfig()
        self.configFile = tempfile.NamedTemporaryFile()
        self.configFile.write(self.config)
        self.configFile.flush()
        self.configFilePath = self.configFile.name
        self.serverCmdLine = """
python server_dev.py
--dont-use-reloader
--config TestConfig
--config-file {}
--port {} """.format(self.configFilePath, self.port)
        return self.serverCmdLine
