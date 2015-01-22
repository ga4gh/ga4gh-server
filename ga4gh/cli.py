"""
Command line interface programs for the GA4GH reference implementation.

TODO: document how to use these for development and simple deployment.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import time
import argparse

import ga4gh.frontend as frontend
import ga4gh.client as client
import ga4gh.backend as backend
import ga4gh.protocol as protocol
import ga4gh.datamodel.variants as variants

##############################################################################
# Server
##############################################################################


def server_main():
    parser = argparse.ArgumentParser(description="GA4GH reference server")
    # Add global options
    parser.add_argument(
        "--port", "-P", default=8000, type=int,
        help="The port to listen on")
    parser.add_argument(
        "--config", "-C", default='DefaultConfig', type=str,
        help="The configuration to use")

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
    wtbParser.set_defaults(variantSetClass=variants.WormtableVariantSet)
    # Tabix
    tabixParser = subparsers.add_parser(
        "tabix",
        description="Serve the API using a tabix based backend.",
        help="Serve data from Tabix indexed VCFs")
    tabixParser.add_argument(
        "dataDir",
        help="The directory containing VCFs")
    tabixParser.set_defaults(variantSetClass=variants.TabixVariantSet)

    args = parser.parse_args()
    if "variantSetClass" not in args:
        parser.print_help()
    else:
        frontend.configure(args.config)
        frontend.app.backend = backend.Backend(
            args.dataDir, args.variantSetClass)
        frontend.app.run(host="0.0.0.0", port=args.port, debug=True)

##############################################################################
# Client
##############################################################################


class AbstractSearchRunner(object):
    """
    Abstract base class for search runner classes
    """
    def __init__(self, args):
        self._workarounds = set(args.workarounds.split(','))
        self._key = args.key

    def usingWorkaroundsFor(self, workaround):
        return workaround in self._workarounds


class VariantSetSearchRunner(AbstractSearchRunner):
    """
    Runner class for the variantsets/search method.
    """
    def __init__(self, args):
        super(VariantSetSearchRunner, self).__init__(args)
        svsr = protocol.GASearchVariantSetsRequest()
        svsr.pageSize = args.pageSize
        self._request = svsr
        self._verbosity = args.verbose
        self._httpClient = client.HTTPClient(
            args.baseUrl, args.verbose, self._workarounds, self._key)

    def run(self):
        for variantSet in self._httpClient.searchVariantSets(self._request):
            print(variantSet.datasetId, variantSet.id)


class VariantSearchRunner(AbstractSearchRunner):
    """
    Runner class for the variants/search method.
    """
    def __init__(self, args):
        super(VariantSearchRunner, self).__init__(args)
        svr = protocol.GASearchVariantsRequest()
        svr.referenceName = args.referenceName
        svr.variantName = args.variantName
        svr.start = args.start
        svr.end = args.end
        svr.pageSize = args.pageSize
        if args.callSetIds == []:
            svr.callSetIds = []
        elif args.callSetIds == '*':
            svr.callSetIds = None
        else:
            svr.callSetIds = args.callSetIds.split(",")
        svr.variantSetIds = args.variantSetIds.split(",")
        self._request = svr
        self._verbosity = args.verbose
        self._httpClient = client.HTTPClient(
            args.baseUrl, args.verbose, self._workarounds, self._key)

    def run(self):
        for variant in self._httpClient.searchVariants(self._request):
            self.printVariant(variant)

    def printVariant(self, variant):
        """
        Prints out the specified GAVariant object in a VCF-like form.
        """
        print(
            variant.id, variant.variantSetId, variant.names,
            variant.referenceName, variant.start, variant.end,
            variant.referenceBases, variant.alternateBases,
            sep="\t", end="\t")
        for key, value in variant.info.items():
            print(key, value, sep="=", end=";")
        print("\t", end="")
        for c in variant.calls:
            print(
                c.callSetId, c.genotype, c.genotypeLikelihood, c.info,
                c.phaseset, sep=":", end="\t")
        print()


class ReferenceSetSearchRunner(AbstractSearchRunner):
    """
    Runner class for the referencesets/search method.
    """
    def __init__(self, args):
        super(ReferenceSetSearchRunner, self).__init__(args)
        request = protocol.GASearchReferenceSetsRequest()
        if args.accessions is not None:
            request.accessions = args.accessions.split(",")
        if self.usingWorkaroundsFor(client.HTTPClient.workaroundGoogle):
            # google says assemblyId not a valid field
            del request.__dict__['assemblyId']
        if args.md5checksums is not None:
            request.md5checksums = args.md5checksums.split(",")
        request.pageSize = args.pageSize
        self._request = request
        self._verbosity = args.verbose
        self._httpClient = client.HTTPClient(
            args.baseUrl, args.verbose, self._workarounds, self._key)

    def run(self):
        for referenceSetId in self._httpClient.searchReferenceSets(
                self._request):
            print(referenceSetId.id)


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
            for variant in self._httpClient.searchVariants(self._request):
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
        "--callSetIds", "-c", default=[],
        help="""Return variant calls which belong to call sets
            with these IDs. Pass in IDs as a comma separated list (no spaces),
            or '*' (with the single quotes!) to indicate 'all call sets'.
            Omit this option to indicate 'no call sets'.
            """)
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
    parser.add_argument(
        "--workarounds", "-w", default='', help="The workarounds to use")
    parser.add_argument(
        "--key", "-k", help="The auth key to use")
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
    bmParser.set_defaults(runner=BenchmarkRunner)
    addUrlArgument(bmParser)
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
    # referencesets/search
    rssParser = subparsers.add_parser(
        "referencesets-search",
        description="Search for referenceSets",
        help="Search for referenceSets")
    rssParser.set_defaults(runner=ReferenceSetSearchRunner)
    addUrlArgument(rssParser)
    rssParser.add_argument(
        "--accessions", default=None,
        help="The accessions to search over")
    rssParser.add_argument(
        "--md5checksums", default=None,
        help="The md5checksums to search over")
    rssParser.add_argument(
        "--assembly-id",
        help="The assembly id to search over")
    rssParser.add_argument(
        "--pageSize", "-m", default=100, type=int,
        help="The maximum number of variants returned in one response.")

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()
