"""
Assists in executing tests with a running server on localhost
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import subprocess
import unittest

import server as server


class ServerTest(unittest.TestCase):
    """
    Manages server bootup and shutdown for a test
    """
    def otherSetup(self):
        """
        We need this method if we want to, say, start up another server
        in the same test before the Ga4gh server because setUp is called
        immediately prior to the running of the test method, which
        is too late.
        """
        pass

    def otherExceptionHandling(self):
        pass

    def otherTeardown(self):
        pass

    def getServer(self):
        return server.Ga4ghServerForTesting()

    def run(self, *args, **kwargs):
        self.otherSetup()
        self.server = self.getServer()
        try:
            self.server.start()
            super(ServerTest, self).run(*args, **kwargs)
        except Exception as e:
            self.server.printDebugInfo()
            self.otherExceptionHandling()
            raise e
        finally:
            self.server.shutdown()
            self.otherTeardown()

    def runClientCmd(self, client, command):
        try:
            client.runCommand(command)
        except subprocess.CalledProcessError as error:
            self.server.printDebugInfo()
            raise error
