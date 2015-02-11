"""
Tests the cli
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import argparse

import ga4gh.cli as cli
import ga4gh.protocol as protocol
import ga4gh.client as client


class TestNoInput(unittest.TestCase):
    """
    Test the cli as if there was no input; print_help should be called
    """
    def testServerNoInput(self):
        parser = StubArgumentParser(self)
        cli.server_main(parser)
        parser.assertParseArgsCountEquals(1)
        parser.assertHelpCountEquals(1)

    def testCliNoInput(self):
        parser = StubArgumentParser(self)
        cli.client_main(parser)
        parser.assertParseArgsCountEquals(1)
        parser.assertHelpCountEquals(1)


class TestClientArguments(unittest.TestCase):
    """
    Tests the client cli can parse all arguments it is supposed to
    and can initialize the runner in preparation for a request
    """
    def setUp(self):
        # initialize the client parser
        self.parser = StubArgumentParser(self)
        cli.client_main(self.parser)

    def verifyInput(self):
        # include arguments common to all commands
        inputStr = "--verbose --workarounds WORK,AROUND --key KEY {0} URL"
        cliInput = inputStr.format(self.cliInput)
        splits = cliInput.split()

        # parse the arguments
        args = self.parser.parser.parse_args(splits)

        # invoke the initializer of the runner (which also parses args)
        runner = args.runner(args)

        # ensure the correct attributes on the runner are set
        if hasattr(runner, '_request'):
            self.assertIsInstance(runner._request, protocol.ProtocolElement)
        self.assertIsInstance(runner._httpClient, client.HttpClient)

    def run(self, *args, **kwargs):
        super(TestClientArguments, self).run(*args, **kwargs)
        self.verifyInput()

    def testVariantsSearchArguments(self):
        self.cliInput = """variants-search --referenceName REFERENCENAME
        --variantName VARIANTNAME --callSetIds CALL,SET,IDS --start 0
        --end 1 --pageSize 2"""

    def testVariantSetsSearchArguments(self):
        self.cliInput = """variantsets-search --pageSize 1 --datasetIds
        DATA,SET,IDS"""

    def testReferenceSetsSearchArguments(self):
        self.cliInput = """referencesets-search --pageSize 1 --accessions
        ACC,ESS,IONS --md5checksums MD5,CHECKSUMS --assemblyId ASSEMBLYID"""

    def testReferencesSearchArguments(self):
        self.cliInput = """references-search --pageSize 1 --accessions
        ACC,ESS,IONS --md5checksums MD5,CHECKSUMS"""

    def testReadGroupSetsSearchArguments(self):
        self.cliInput = """readgroupsets-search --pageSize 1 --datasetIds
        DATA,SET,IDS --name NAME"""

    def testCallSetsSearchArguments(self):
        self.cliInput = """callsets-search --pageSize 1 --name NAME
        --variantSetIds VARIANT,SET,IDS"""

    def testReadsSearchArguments(self):
        self.cliInput = """reads-search --pageSize 1 --start 2 --end 3
        --readGroupIds READ,GROUP,IDS --referenceId REFERENCEID
        --referenceName REFERENCENAME"""

    def testReferenceSetGetArguments(self):
        self.cliInput = """referencesets-get --id ID"""

    def testReferenceGetArguments(self):
        self.cliInput = """references-get --id ID"""

    def testReferenceBasesListArguments(self):
        self.cliInput = """references-list-bases --id ID
        --start 1 --end 2"""


class StubArgumentParser(object):
    """
    A stand-in object for an ArgumentParser that intercepts calls
    to parse_args and print_help, but otherwise provides normal
    behavior via passing through calls to an attribute ArgumentParser
    """
    def __init__(self, currentTest):
        self.parser = argparse.ArgumentParser(description="Stub")
        self.helpCount = 0
        self.parseArgsCount = 0
        self.parse_args = self._parseArgs
        self.print_help = self._help
        self.currentTest = currentTest

    def _help(self):
        self.helpCount += 1

    def _parseArgs(self):
        self.parseArgsCount += 1
        return ()

    def assertParseArgsCountEquals(self, parseArgsCount):
        self.currentTest.assertEquals(self.parseArgsCount, parseArgsCount)

    def assertHelpCountEquals(self, helpCount):
        self.currentTest.assertEquals(self.helpCount, helpCount)

    def __getattr__(self, name):
        return getattr(self.parser, name)
