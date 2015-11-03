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


class TestClientOutput(unittest.TestCase):
    """
    Base class for client output tests
    """
    def setUp(self):
        self._dataDir = "tests/data"
        self._dataUrl = "file://{}".format(self._dataDir)
        self._backend = backend.FileSystemBackend(self._dataDir)
        self._client = client.LocalClient(self._backend)

    def captureCliOutput(self, command, arguments, outputFormat):
        clientCommand = "{} {} {} -O {}".format(
            command, self._dataUrl, arguments, outputFormat)
        stdout, stderr = utils.captureOutput(
            cli.client_main, clientCommand.split())
        self.assertEqual(len(stderr), 0)
        return stdout


class TestClientFasta(TestClientOutput):
    """
    Tests client FASTA output
    """
    def captureFastaOutput(self, command, arguments=""):
        stdout = self.captureCliOutput(command, arguments, "fasta")
        lines = stdout.split()
        return lines

    def testListReferenceBases(self):
        referenceSetIterator = self._client.searchReferenceSets()
        referenceSet = next(referenceSetIterator)
        referencesIterator = self._client.searchReferences(referenceSet.id)
        reference = next(referencesIterator)
        start = 1
        end = 5
        lines = self.captureFastaOutput(
            "references-list-bases --start {} --end {}".format(start, end),
            reference.id)
        self.assertEqual(
            lines[0], ">{}:{}-{}".format(reference.id, start, end))
        cliBases = ''.join(lines[1:])
        bases = self._client.listReferenceBases(reference.id, start, end)
        self.assertEqual(cliBases, bases)


class TestClientJson(TestClientOutput):
    """
    Tests that the JSON output by the client on the command line for
    various options is equal to the values we find using the Python
    client API.
    """
    def captureJsonOutput(self, command, arguments=""):
        """
        Runs the specified command add the JSON output option and
        returns the result as a list of JSON parsed dictionaries.
        """
        stdout = self.captureCliOutput(command, arguments, "json")
        cliOutput = []
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

    def testSearchAllReadGroups(self):
        # TODO: add more rigorous testing here
        cliOutput = self.captureJsonOutput("reads-search")
        self.assertGreater(len(cliOutput), 0)
