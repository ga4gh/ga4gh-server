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


class BenchmarkRunner(object):
    """
    Runner class for the client side benchmarking. This is intended to give
    rough figures on protocol throughput on the server side over various
    requests.
    """
    def __init__(self, args):
        self._maxResults = args.maxResults
        self._verbosity = args.verbose
        self._httpClient = ga4gh.client.HTTPClient(
            args.hostname, args.port, args.verbose)

    def run(self):
        svr = ga4gh.protocol.GASearchVariantsRequest()
        svr.start = 0
        svr.end = 2**32
        numVariants = 0
        beforeCpu = time.clock()
        beforeWall = time.time()
        try:
            for v in self._httpClient.searchVariants(svr):
                numVariants += 1
        except KeyboardInterrupt:
            pass
        cpuTime = time.clock() - beforeCpu
        wallTime = time.time() - beforeWall
        totalBytes = self._httpClient.getBytesRead()
        totalBytes /= 1024 * 1024
        s = "read {0} variants in {1:.2f} seconds; CPU time {2:.2f}".format(
            numVariants, wallTime, cpuTime)
        s += "; {0:.2f} MB @ {1:.2f} MB/s".format(
            totalBytes, totalBytes / wallTime)
        print(s)


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
        svr.maxResults = args.maxResults
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
        # TODO insert info fields
        for kv in v.info:
            print(kv.key, kv.value, sep="=", end=";")
        print("\t", end="")
        for c in v.calls:
            print(c.genotype, c.genotypeLikelihood, sep=":", end="\t")
        print()


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
    vsParser.add_argument(
        "--referenceName", "-r", default="chrSim",
        help="Only return variants on this reference.")
    vsParser.add_argument(
        "--variantName", "-n", default=None,
        help="Only return variants which have exactly this name. TODO")
    vsParser.add_argument(
        "--start", "-s", default=0, type=int,
        help="The start of the search range (inclusive).")
    vsParser.add_argument(
        "--end", "-e", default=1, type=int,
        help="The end of the search range (exclusive).")
    vsParser.add_argument(
        "--maxResults", "-m", default=10, type=int,
        help="The maximum number of variants returned in one response.")
    # benchmarking
    bmParser = subparsers.add_parser(
        "benchmark",
        description="Run simple benchmarks on the various methods",
        help="Benchmark server performance")
    bmParser.set_defaults(runner=BenchmarkRunner)
    bmParser.add_argument(
        "--maxResults", "-m", default=10, type=int,
        help="The maximum number of variants returned in one response.")

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

if __name__ == "__main__":
    main()
