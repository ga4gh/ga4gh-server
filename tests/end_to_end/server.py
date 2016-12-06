"""
Servers to assist in testing
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import logging
import tempfile
import shlex
import subprocess
import socket

import requests

import ga4gh.common.utils as utils


ga4ghPort = 8001
remotePort = 8002
oidcOpPort = 8443


class ServerForTesting(object):
    """
    The base class of a test server
    """
    def __init__(self, port, protocol='http',
                 subdirectory=None, pingStatusCode=200):
        # suppress requests package log messages
        logging.getLogger("requests").setLevel(logging.CRITICAL)
        self.port = port
        self.subdirectory = subdirectory
        self.pingStatusCode = pingStatusCode
        self.outFile = None
        self.errFile = None
        self.server = None
        self.serverUrl = "{}://{}:{}".format(protocol,
                                             socket.gethostname(),
                                             self.port)

    def getUrl(self):
        """
        Return the url at which the server is configured to run
        """
        return self.serverUrl

    def getCmdLine(self):
        """
        Return the command line string used to launch the server.
        Subclasses must override this method.
        """
        raise NotImplementedError()

    def start(self):
        """
        Start the server
        """
        assert not self.isRunning(), "Another server is running at {}".format(
            self.serverUrl)
        self.outFile = tempfile.TemporaryFile()
        self.errFile = tempfile.TemporaryFile()
        splits = shlex.split(self.getCmdLine())
        self.server = subprocess.Popen(
            splits, stdout=self.outFile,
            stderr=self.errFile,
            cwd=self.subdirectory)
        self._waitForServerStartup()

    def shutdown(self):
        """
        Shut down the server
        """
        if self.isRunning():
            self.server.kill()
        if self.server is not None:
            self.server.wait()
            self._assertServerShutdown()
        if self.outFile is not None:
            self.outFile.close()
        if self.errFile is not None:
            self.errFile.close()

    def restart(self):
        """
        Restart the server
        """
        self.shutdown()
        self.start()

    def isRunning(self):
        """
        Returns true if the server is running, false otherwise
        """
        try:
            response = self.ping()
            if response.status_code != self.pingStatusCode:
                msg = ("Ping of server {} returned unexpected status code "
                       "({})").format(self.serverUrl, response.status_code)
                assert False, msg
            return True
        except requests.ConnectionError:
            return False

    def ping(self):
        """
        Pings the server by doing a GET request to /
        """
        response = requests.get(self.serverUrl, verify=False)
        return response

    def getOutLines(self):
        """
        Return the lines of the server stdout file
        """
        return utils.getLinesFromLogFile(self.outFile)

    def getErrLines(self):
        """
        Return the lines of the server stderr file
        """
        return utils.getLinesFromLogFile(self.errFile)

    def printDebugInfo(self):
        """
        Print debugging information about the server
        """
        className = self.__class__.__name__
        print('\n')
        print('*** {} CMD ***'.format(className))
        print(self.getCmdLine())
        print('*** {} STDOUT ***'.format(className))
        print(''.join(self.getOutLines()))
        print('*** {} STDERR ***'.format(className))
        print(''.join(self.getErrLines()))

    @utils.Timeout()
    @utils.Repeat()
    def _waitForServerStartup(self):
        self.server.poll()
        if self.server.returncode is not None:
            self._waitForErrLines()
            message = "Server process unexpectedly died; stderr: {0}"
            failMessage = message.format(''.join(self.getErrLines()))
            assert False, failMessage
        return not self.isRunning()

    @utils.Timeout()
    @utils.Repeat()
    def _waitForErrLines(self):
        # not sure why there's some delay in getting the server
        # process' stderr (at least for the ga4gh server)...
        return self.getErrLines() == []

    def _assertServerShutdown(self):
        shutdownString = "Server did not shut down correctly"
        assert self.server.returncode is not None, shutdownString
        assert not self.isRunning(), shutdownString


class Ga4ghServerForTesting(ServerForTesting):
    """
    A ga4gh test server
    """
    def __init__(self, useOidc=False):
        protocol = 'https' if useOidc else 'http'
        super(Ga4ghServerForTesting, self).__init__(ga4ghPort, protocol)
        self.configFile = None
        self.useOidc = useOidc

    def getConfig(self):
        config = """
SIMULATED_BACKEND_NUM_VARIANT_SETS = 10
SIMULATED_BACKEND_VARIANT_DENSITY = 1
DATA_SOURCE = "simulated://"
DEBUG = True
"""
        if self.useOidc:
            config += """
TESTING = True
OIDC_PROVIDER = "https://localhost:{0}"
SECRET_KEY = "secret"
""".format(oidcOpPort)
        return config

    def getCmdLine(self):
        if self.configFile is None:
            self.configFile = tempfile.NamedTemporaryFile()
        config = self.getConfig()
        self.configFile.write(config)
        self.configFile.flush()
        configFilePath = self.configFile.name
        cmdLine = """
python server_dev.py
--dont-use-reloader
--disable-urllib-warnings
--host 0.0.0.0
--config TestConfig
--config-file {}
--port {} """.format(configFilePath, self.port)
        return cmdLine

    def shutdown(self):
        super(Ga4ghServerForTesting, self).shutdown()
        if self.configFile is not None:
            self.configFile.close()

    def printDebugInfo(self):
        super(Ga4ghServerForTesting, self).printDebugInfo()
        className = self.__class__.__name__
        print('*** {} CONFIG ***'.format(className))
        print(self.getConfig())


class Ga4ghServerForTestingDataSource(Ga4ghServerForTesting):
    """
    A test server that reads data from a data source
    """
    def __init__(self, dataDir):
        super(Ga4ghServerForTestingDataSource, self).__init__()
        self.dataDir = dataDir

    def getConfig(self):
        config = """
DATA_SOURCE = "{}"
DEBUG = True""".format(self.dataDir)
        return config


class OidcOpServerForTesting(ServerForTesting):
    """
    Runs a test OP server on localhost
    """
    def __init__(self):
        super(OidcOpServerForTesting, self).__init__(
            oidcOpPort, protocol="https",
            subdirectory="oidc-provider/simple_op",
            pingStatusCode=404)

    def getCmdLine(self):
        return ("python src/run.py --base https://localhost:{}" +
                " -p {} -d settings.yaml").format(oidcOpPort, oidcOpPort)
