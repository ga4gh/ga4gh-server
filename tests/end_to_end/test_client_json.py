"""
Performs an end-to-end test where we verify that the data
output by the client command line interface is equal to
the values we expect using the test dataset.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import unittest

import ga4gh.client as client
import ga4gh.backend as backend
import ga4gh.cli as cli
import tests.utils as utils


class TestClientJson(unittest.TestCase):
    """
    Tests that the JSON output by the client on the command line for
    various options is equal to the values we find using the Python
    client API.
    """
    def setUp(self):
        self._dataDir = "tests/data"
        self._dataUrl = "file://{}".format(self._dataDir)
        self._backend = backend.FileSystemBackend(self._dataDir)
        self._client = client.LocalClient(self._backend)

    def captureJsonOutput(self, command, arguments=""):
        """
        Runs the specified command add the JSON output option and
        returns the result as a list of JSON parsed dictionaries.
        """
        clientCommand = "{} {} {} -O json".format(
            command, self._dataUrl, arguments)
        stdout, stderr = utils.captureOutput(
            cli.client_main, clientCommand.split())
        cliOutput = []
        self.assertEqual(len(stderr), 0)
        for line in stdout.splitlines():
            cliOutput.append(json.loads(line))
        return cliOutput

    def verifyParsedOutputsEqual(
            self, clientIterator, cliCommand, cliArguments=""):
        """
        Verify that the parsed JSON of all the objects in the specified
        client iterator are equal to the parsed JSON from the specified
        CLI command.
        """
        cliOutput = self.captureJsonOutput(cliCommand, cliArguments)
        clientOutput = [gaObject.toJsonDict() for gaObject in clientIterator]
        self.assertEqual(clientOutput, cliOutput)

    def testSearchAllDatasets(self):
        iterator = self._client.searchDatasets()
        self.verifyParsedOutputsEqual(iterator, "datasets-search")

    def testSearchAllReferenceSets(self):
        iterator = self._client.searchReferenceSets()
        self.verifyParsedOutputsEqual(iterator, "referencesets-search")

    def testSearchAllReferences(self):
        for referenceSet in self._client.searchReferenceSets():
            iterator = self._client.searchReferences(
                referenceSetId=referenceSet.id)
            args = "--referenceSetId={}".format(referenceSet.id)
            self.verifyParsedOutputsEqual(iterator, "references-search", args)

    def testGetDataset(self):
        for dataset in self._client.searchDatasets():
            self.verifyParsedOutputsEqual(
                [dataset], "datasets-get", dataset.id)
