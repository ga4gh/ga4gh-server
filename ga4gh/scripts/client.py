"""
Command line interface for the ga4gh reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import time
import argparse

import ga4gh
import ga4gh.client
import ga4gh.protocol


class VariantSetSearchRunner(object):
    """
    Runner class for the variantsets/search method.
    """
    def __init__(self, args):
        svsr = ga4gh.protocol.GASearchVariantSetsRequest()
        svsr.pageSize = args.pageSize
        self._request = svsr
        self._verbosity = args.verbose
        self._httpClient = ga4gh.client.HTTPClient(
            args.hostname, args.port, args.verbose)

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
        self._httpClient = ga4gh.client.HTTPClient(
            args.hostname, args.port, args.verbose)

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
        for kv in v.info:
            print(kv.key, kv.value, sep="=", end=";")
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


def main():
    parser = argparse.ArgumentParser(description="GA4GH reference client")
    # Add global options
    parser.add_argument(
        "--hostname", "-H", default="localhost", help="The server host")
    parser.add_argument(
        "--port", "-P", default=8000, type=int, help="The server port")
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
    addOptions(vsParser)
    # benchmarking
    bmParser = subparsers.add_parser(
        "benchmark",
        description="Run simple benchmarks on the various methods",
        help="Benchmark server performance")
    bmParser.set_defaults(runner=BenchmarkRunner)
    addOptions(bmParser)
    # variantsets/search
    vssParser = subparsers.add_parser(
        "variantsets-search",
        description="Search for variantSets",
        help="Search for variantSets.")
    vssParser.set_defaults(runner=VariantSetSearchRunner)
    vssParser.add_argument(
        "--pageSize", "-m", default=100, type=int,
        help="The maximum number of variants returned in one response.")

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

if __name__ == "__main__":
    main()
