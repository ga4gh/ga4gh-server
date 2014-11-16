"""
Command line interface programs for the GA4GH reference implementation.

TODO: document how to use these for development and simple deployment.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import time
import argparse
import werkzeug.serving

import ga4gh
import ga4gh.server
import ga4gh.client
import ga4gh.protocol

##############################################################################
# Server
##############################################################################

class ServerRunner(object):
    """
    Superclass of server runner; takes care of functionality common to
    all backends.
    """
    def __init__(self, args):
        backend = self.getBackend(args)
        self._port = args.port
        self._httpHandler = ga4gh.server.HTTPHandler(backend)

    def run(self):
        werkzeug.serving.run_simple(
            '', self._port, self._httpHandler.wsgiApplication,
            use_reloader=True)


class WormtableRunner(ServerRunner):
    """
    Runner class to run the server using the wormtable based backend.
    """
    def getBackend(self, args):
        backend = ga4gh.server.WormtableBackend(args.dataDir)
        return backend


class TabixRunner(ServerRunner):
    """
    Runner class to start the server using a tabix backend.
    """
    def getBackend(self, args):
        backend = ga4gh.server.TabixBackend(args.dataDir)
        return backend


def server_main():
    parser = argparse.ArgumentParser(description="GA4GH reference server")
    # Add global options
    parser.add_argument(
        "--port", "-P", default=8000, type=int,
        help="The port to listen on")
    parser.add_argument('--verbose', '-v', action='count', default=0)
    subparsers = parser.add_subparsers(title='subcommands',)

    # help
    helpParser = subparsers.add_parser(
        "help",
        description="ga4gh_server help",
        help="show this help message and exit")
    # Wormtable backend
    wtbParser = subparsers.add_parser(
        "wormtable",
        description="Serve the API using a wormtable based backend.",
        help="Serve data from tables.")
    wtbParser.add_argument(
        "dataDir",
        help="The directory containing the wormtables to be served.")
    wtbParser.set_defaults(runner=WormtableRunner)
    # Tabix
    tabixParser = subparsers.add_parser(
        "tabix",
        description="Serve the API using a tabix based backend.",
        help="Serve data from Tabix indexed VCFs")
    tabixParser.add_argument(
        "dataDir",
        help="The directory containing VCFs")
    tabixParser.set_defaults(runner=TabixRunner)

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

##############################################################################
# Client
##############################################################################

class VariantSetSearchRunner(object):
    """
    Runner class for the variantsets/search method.
    """
    def __init__(self, args):
        svsr = ga4gh.protocol.GASearchVariantSetsRequest()
        svsr.pageSize = args.pageSize
        self._request = svsr
        self._verbosity = args.verbose
        self._httpClient = ga4gh.client.HTTPClient(args.baseUrl, args.verbose)

    def run(self):
        for v in self._httpClient.searchVariantSets(self._request):
            print(v.datasetId, v.id)


class VariantSearchRunner(object):
    """
    Runner class for the variants/search method.
    """
    def __init__(self, args):
        svr = ga4gh.protocol.GASearchVariantsRequest()
        svr.referenceName = args.referenceName
        svr.variantName = args.variantName
        svr.start = args.start
        svr.end = args.end
        svr.pageSize = args.pageSize
        if args.callSetIds is not None:
            svr.callSetIds = args.callSetIds.split(",")
        svr.variantSetIds = args.variantSetIds.split(",")
        self._request = svr
        self._verbosity = args.verbose
        self._httpClient = ga4gh.client.HTTPClient(args.baseUrl, args.verbose)

    def run(self):
        for v in self._httpClient.searchVariants(self._request):
            self.printVariant(v)

    def printVariant(self, v):
        """
        Prints out the specified GAVariant object in a VCF-like form.
        """
        print(
            v.id, v.variantSetId, v.names,
            v.referenceName, v.start, v.end, v.referenceBases,
            v.alternateBases, sep="\t", end="\t")
        for key, value in v.info.items():
            print(key, value, sep="=", end=";")
        print("\t", end="")
        for c in v.calls:
            print(
                c.callSetId, c.genotype, c.genotypeLikelihood, c.info,
                c.phaseset, sep=":", end="\t")
        print()


class BenchmarkRunner(VariantSearchRunner):
    """
    Runner class for the client side benchmarking. This is intended to give
    rough figures on protocol throughput on the server side over various
    requests.
    """
    def run(self):
        numVariants = 0
        beforeCpu = time.clock()
        beforeWall = time.time()
        try:
            for v in self._httpClient.searchVariants(self._request):
                numVariants += 1
        except KeyboardInterrupt:
            pass
        cpuTime = time.clock() - beforeCpu
        wallTime = time.time() - beforeWall
        totalBytes = self._httpClient.getBytesRead()
        totalBytes /= 1024 * 1024
        s = "read {0} variants in {1:.2f} seconds; CPU time {2:.2f}".format(
            numVariants, wallTime, cpuTime)
        s += "; {0:.2f} MB @ {1:.2f} MB/s; {2:.2f} vars/s".format(
            totalBytes, totalBytes / wallTime, numVariants / wallTime)
        print(s)


def addOptions(parser):
    """
    Adds common options to a command line parser.
    """
    parser.add_argument(
        "variantSetIds",
        help="The variant set id(s) to search over")
    parser.add_argument(
        "--referenceName", "-r", default="chrSim",
        help="Only return variants on this reference.")
    parser.add_argument(
        "--variantName", "-n", default=None,
        help="Only return variants which have exactly this name.")
    parser.add_argument(
        "--callSetIds", "-c", default=None,
        help="""Only return variant calls which belong to call sets
            with these IDs (comma seperated list).""")
    parser.add_argument(
        "--start", "-s", default=0, type=int,
        help="The start of the search range (inclusive).")
    parser.add_argument(
        "--end", "-e", default=1, type=int,
        help="The end of the search range (exclusive).")
    parser.add_argument(
        "--pageSize", "-m", default=100, type=int,
        help="The maximum number of variants returned in one response.")


def addUrlArgument(parser):
    """
    Adds the URL endpoint argument to the specified parser.
    """
    parser.add_argument("baseUrl", help="The URL of the API endpoint")


def client_main():
    parser = argparse.ArgumentParser(description="GA4GH reference client")
    # Add global options
    parser.add_argument('--verbose', '-v', action='count', default=0)
    subparsers = parser.add_subparsers(title='subcommands',)

    # help
    helpParser = subparsers.add_parser(
        "help", description="ga4gh_client help",
        help="show this help message and exit")
    # variants/search
    vsParser = subparsers.add_parser(
        "variants-search",
        description="Search for variants",
        help="Search for variants.")
    vsParser.set_defaults(runner=VariantSearchRunner)
    addUrlArgument(vsParser)
    addOptions(vsParser)
    # benchmarking
    bmParser = subparsers.add_parser(
        "benchmark",
        description="Run simple benchmarks on the various methods",
        help="Benchmark server performance")
    addUrlArgument(bmParser)
    bmParser.set_defaults(runner=BenchmarkRunner)
    addOptions(bmParser)
    # variantsets/search
    vssParser = subparsers.add_parser(
        "variantsets-search",
        description="Search for variantSets",
        help="Search for variantSets.")
    vssParser.set_defaults(runner=VariantSetSearchRunner)
    addUrlArgument(vssParser)
    vssParser.add_argument(
        "--pageSize", "-m", default=100, type=int,
        help="The maximum number of variants returned in one response.")

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()
