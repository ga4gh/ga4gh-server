"""
Tests the cli
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import unittest
import mock

import ga4gh.cli as cli
import ga4gh.protocol as protocol
import ga4gh.client as client


class TestNoInput(unittest.TestCase):
    """
    Test the cli as if there was no input; print_help should be called.
    The development server is a special case here, as it works without
    arguments.
    """
    def setUp(self):
        self.parser = StubArgumentParser(self)

    def run(self, *args, **kwargs):
        super(TestNoInput, self).run(*args, **kwargs)
        self.verifyInput()

    def verifyInput(self):
        self.parser.assertParseArgsCountEquals(1)
        self.parser.assertHelpCountEquals(1)

    def testGa2VcfNoInput(self):
        cli.ga2vcf_main(self.parser)

    def testGa2SamNoInput(self):
        cli.ga2sam_main(self.parser)

    def testCliNoInput(self):
        cli.client_main(self.parser)


class TestGa2VcfArguments(unittest.TestCase):
    """
    Tests the ga2vcf cli can parse all arguments it is supposed to
    """
    def testVariantsSearchArguments(self):
        cliInput = """--workarounds WORK,AROUND
        --key KEY -O --outputFile /dev/null
        --referenceName REFERENCENAME
        --variantName VARIANTNAME --callSetIds CALL,SET,IDS --start 0
        --end 1 --pageSize 2 BASEURL VARIANTSETID"""
        stubConverterModule = StubConverterModuleVcf(self)
        with mock.patch(
                'ga4gh.converters.VcfConverter',
                stubConverterModule):
            parser = StubArgumentParserCli(self, cliInput)
            cli.ga2vcf_main(parser)
            parser.assertParseArgsCountEquals(1)
            stubConverterModule.assertVcfConvertCountEquals(1)


class TestGa2SamArguments(unittest.TestCase):
    """
    Tests the ga2sam cli can parse all arguments it is supposed to
    """
    def testReadsSearchArguments(self):
        cliInput = """--workarounds WORK,AROUND --key KEY -O
        --pageSize 1 --start 2 --end 3 --outputFile OUT.SAM
        --readGroupIds READ,GROUP,IDS --referenceId REFERENCEID
        --binaryOutput BASEURL"""
        stubConverterModule = StubConverterModuleSam(self)
        with mock.patch(
                'ga4gh.converters.SamConverter',
                stubConverterModule):
            parser = StubArgumentParserCli(self, cliInput)
            cli.ga2sam_main(parser)
            parser.assertParseArgsCountEquals(1)
            stubConverterModule.assertSamConvertCountEquals(1)


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
        --end 1 --pageSize 2 --variantSetId VARIANTSETIDS"""

    def testVariantSetsSearchArguments(self):
        self.cliInput = """variantsets-search --pageSize 1 --datasetId
        DATASETID"""

    def testReferenceSetsSearchArguments(self):
        self.cliInput = """referencesets-search --pageSize 1 --accessions
        ACC,ESS,IONS --md5checksums MD5,CHECKSUMS --assemblyId ASSEMBLYID"""

    def testReferencesSearchArguments(self):
        self.cliInput = """references-search --pageSize 1 --accessions
        ACC,ESS,IONS --md5checksums MD5,CHECKSUMS"""

    def testReadGroupSetsSearchArguments(self):
        self.cliInput = """readgroupsets-search --pageSize 1 --datasetId
        DATASETID --name NAME"""

    def testCallSetsSearchArguments(self):
        self.cliInput = """callsets-search --pageSize 1 --name NAME
        --variantSetId VARIANTSETID"""

    def testReadsSearchArguments(self):
        self.cliInput = """reads-search --pageSize 1 --start 2 --end 3
        --readGroupIds READ,GROUP,IDS --referenceId REFERENCEID"""

    def testDatasetsSearchArguments(self):
        self.cliInput = """datasets-search"""

    def testReferenceSetGetArguments(self):
        self.cliInput = """referencesets-get ID"""

    def testReferenceGetArguments(self):
        self.cliInput = """references-get ID"""

    def testReadGroupSetGetArguments(self):
        self.cliInput = """readgroupsets-get ID"""

    def testReadGroupGetArguments(self):
        self.cliInput = """readgroups-get ID"""

    def testCallSetGetArguments(self):
        self.cliInput = """callsets-get ID"""

    def testDatasetsGetArguments(self):
        self.cliInput = """datasets-get ID"""

    def testVariantGetArguments(self):
        self.cliInput = """variants-get ID"""

    def testReferenceBasesListArguments(self):
        self.cliInput = """references-list-bases ID
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


class StubArgumentParserCli(StubArgumentParser):
    """
    Like StubArgumentParser, but returns real arguments from the
    parse_args call (that the user provides)
    """
    def __init__(self, currentTest, cliInput):
        super(StubArgumentParserCli, self).__init__(currentTest)
        self.cliInput = cliInput

    def _parseArgs(self):
        self.parseArgsCount += 1
        splits = self.cliInput.split()
        return self.parser.parse_args(splits)


class StubConverterModule(object):
    """
    A stand-in object for the converter module.
    Just provides access to dummy objects.
    """
    def __init__(self, currentTest):
        self.currentTest = currentTest
        self.VcfConverter = StubConverter(self.currentTest)
        self.SamConverter = StubConverter(self.currentTest)

    def assertVcfConvertCountEquals(self, convertCount):
        self.VcfConverter.assertConvertCountEquals(convertCount)

    def assertSamConvertCountEquals(self, convertCount):
        self.SamConverter.assertConvertCountEquals(convertCount)


class StubConverterModuleSam(StubConverterModule):
    """
    The StubConverterModule for Sam tests
    """
    def __call__(self, *args):
        return self.SamConverter


class StubConverterModuleVcf(StubConverterModule):
    """
    The StubConverterModule for Vcf tests
    """
    def __call__(self, *args):
        return self.VcfConverter


class StubConverter(object):
    """
    A stand-in object for a converter that does nothing
    """
    def __init__(self, currentTest):
        self.currentTest = currentTest
        self.convertCount = 0

    def convert(self):
        self.convertCount += 1

    def assertConvertCountEquals(self, convertCount):
        self.currentTest.assertEquals(self.convertCount, convertCount)
