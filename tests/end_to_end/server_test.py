"""
Assists in executing tests with a running server on localhost
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import subprocess

import server as server


class ClientHelperMixin(object):
    """
    Helper methods involving the client for server tests
    """
    def runClientCmd(self, client, command, arguments=""):
        try:
            client.runCommand(command, arguments)
        except subprocess.CalledProcessError as error:
            self.server.printDebugInfo()
            raise error


class ServerTest(ClientHelperMixin, unittest.TestCase):
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
        except Exception as exception:
            self.server.printDebugInfo()
            self.otherExceptionHandling()
            raise exception
        finally:
            self.server.shutdown()
            self.otherTeardown()


class ServerTestClass(ClientHelperMixin, unittest.TestCase):
    """
    Like ServerTest, except starts and stops the server at the class
    level instead of the method level.
    """
    @classmethod
    def setUpClass(cls):
        cls.otherSetup()
        cls.server = cls.getServer()
        cls.server.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.otherTeardown()

    @classmethod
    def otherSetup(cls):
        pass

    @classmethod
    def otherTeardown(cls):
        pass

    @classmethod
    def otherExceptionHandling(cls):
        pass

    @classmethod
    def getServer(cls):
        return server.Ga4ghServerForTesting()

    def run(self, *args, **kwargs):
        try:
            super(ServerTestClass, self).run(*args, **kwargs)
        except Exception as exception:
            self.server.printDebugInfo()
            self.otherExceptionHandling()
            raise exception
