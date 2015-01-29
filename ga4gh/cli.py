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


def server_main(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="GA4GH reference server")
    # Add global options
    parser.add_argument(
        "--port", "-P", default=8000, type=int,
        help="The port to listen on")
    parser.add_argument(
        "--config", "-C", default='DefaultConfig', type=str,
        help="The configuration to use")
    parser.add_argument(
        "--config-file", "-F", type=str,
        help="The configuration file to use")

    subparsers = parser.add_subparsers(title='subcommands',)

    # help
    subparsers.add_parser(
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
        frontend.configure(args.config, args.config_file)
        frontend.app.backend = backend.Backend(
            args.dataDir, args.variantSetClass)
        frontend.app.run(host="0.0.0.0", port=args.port, debug=True)

##############################################################################
# Client
##############################################################################


class AbstractQueryRunner(object):
    """
    Abstract base class for runner classes
    """
    def __init__(self, args):
        self._workarounds = set(args.workarounds.split(','))
        self._key = args.key
        self._verbosity = args.verbose
        self._httpClient = client.HttpClient(
            args.baseUrl, args.verbose, self._workarounds, self._key)

    def usingWorkaroundsFor(self, workaround):
        """
        Returns true if we are using the passed-in workaround
        """
        return workaround in self._workarounds


class AbstractGetRunner(AbstractQueryRunner):
    """
    Abstract base class for get runner classes
    """
    def __init__(self, args):
        super(AbstractGetRunner, self).__init__(args)
        self._id = args.id
        self._httpClient = client.HttpClient(
            args.baseUrl, args.verbose, self._workarounds, self._key)

    def _run(self, method):
        response = method(self._id)
        print(response.id)


class AbstractSearchRunner(AbstractQueryRunner):
    """
    Abstract base class for search runner classes
    """
    def __init__(self, args):
        super(AbstractSearchRunner, self).__init__(args)

    def _setRequest(self, request, args):
        """
        Sets the _httpClient and other common attributes
        """
        self._minimalOutput = args.minimalOutput
        if 'pageSize' in args:
            # GAListReferenceBasesRequest does not have a pageSize attr
            request.pageSize = args.pageSize
        self._request = request

    def _run(self, method, attrName=None):
        """
        Runs the request given methodname and prints out
        the each result's attrName attribute if it is provided.
        If not, prints each entire result object.
        """
        results = method(self._request)
        for result in results:
            if attrName is None:
                print(result)
            else:
                attr = getattr(result, attrName)
                print(attr)


class SearchVariantSetsRunner(AbstractSearchRunner):
    """
    Runner class for the variantsets/search method.
    """
    def __init__(self, args):
        super(SearchVariantSetsRunner, self).__init__(args)
        request = protocol.GASearchVariantSetsRequest()
        setCommaSeparatedAttribute(request, args, 'datasetIds')
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchVariantSets, 'datasetId')


class SearchVariantsRunner(AbstractSearchRunner):
    """
    Runner class for the variants/search method.
    """
    def __init__(self, args):
        super(SearchVariantsRunner, self).__init__(args)
        request = protocol.GASearchVariantsRequest()
        request.referenceName = args.referenceName
        request.variantName = args.variantName
        request.start = args.start
        request.end = args.end
        if self.usingWorkaroundsFor(client.HttpClient.workaroundGoogle):
            request.maxCalls = args.maxCalls
        if args.callSetIds == []:
            request.callSetIds = []
        elif args.callSetIds == '*':
            request.callSetIds = None
        else:
            request.callSetIds = args.callSetIds.split(",")
        setCommaSeparatedAttribute(request, args, 'variantSetIds')
        self._setRequest(request, args)

    def run(self):
        if self._minimalOutput:
            self._run(self._httpClient.searchVariants, 'id')
        else:
            results = self._httpClient.searchVariants(self._request)
            for result in results:
                self.printVariant(result)

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


class SearchReferenceSetsRunner(AbstractSearchRunner):
    """
    Runner class for the referencesets/search method.
    """
    def __init__(self, args):
        super(SearchReferenceSetsRunner, self).__init__(args)
        request = protocol.GASearchReferenceSetsRequest()
        setCommaSeparatedAttribute(request, args, 'accessions')
        setCommaSeparatedAttribute(request, args, 'md5checksums')
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchReferenceSets, 'id')


class SearchReferencesRunner(AbstractSearchRunner):
    """
    Runner class for the references/search method
    """
    def __init__(self, args):
        super(SearchReferencesRunner, self).__init__(args)
        request = protocol.GASearchReferencesRequest()
        setCommaSeparatedAttribute(request, args, 'accessions')
        setCommaSeparatedAttribute(request, args, 'md5checksums')
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchReferences, 'id')


class SearchReadGroupSetsRunner(AbstractSearchRunner):
    """
    Runner class for the readgroupsets/search method
    """
    def __init__(self, args):
        super(SearchReadGroupSetsRunner, self).__init__(args)
        request = protocol.GASearchReadGroupSetsRequest()
        setCommaSeparatedAttribute(request, args, 'datasetIds')
        request.name = args.name
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchReadGroupSets, 'id')


class SearchCallSetsRunner(AbstractSearchRunner):
    """
    Runner class for the callsets/search method
    """
    def __init__(self, args):
        super(SearchCallSetsRunner, self).__init__(args)
        request = protocol.GASearchCallSetsRequest()
        setCommaSeparatedAttribute(request, args, 'variantSetIds')
        request.name = args.name
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchCallSets, 'id')


class SearchReadsRunner(AbstractSearchRunner):
    """
    Runner class for the reads/search method
    """
    class GASearchReadsRequestGoogle(protocol.ProtocolElement):

        __slots__ = ['end', 'pageSize', 'pageToken', 'readGroupIds',
                     'referenceName', 'start']

        def __init__(self):
            self.end = None
            self.pageSize = None
            self.pageToken = None
            self.readGroupIds = []
            self.referenceName = None
            self.start = 0

    def __init__(self, args):
        super(SearchReadsRunner, self).__init__(args)
        request = protocol.GASearchReadsRequest()
        if self.usingWorkaroundsFor(client.HttpClient.workaroundGoogle):
            # google says referenceId not a valid field
            request = self.GASearchReadsRequestGoogle()
        setCommaSeparatedAttribute(request, args, 'readGroupIds')
        request.start = args.start
        request.end = args.end
        request.referenceId = args.referenceId
        request.referenceName = args.referenceName
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchReads, 'id')


class ListReferenceBasesRunner(AbstractSearchRunner):
    """
    Runner class for the references/{id}/bases method
    """
    def __init__(self, args):
        super(ListReferenceBasesRunner, self).__init__(args)
        request = protocol.GAListReferenceBasesRequest()
        request.start = args.start
        request.end = args.end
        self._id = args.id
        self._setRequest(request, args)

    def run(self):
        method = self._httpClient.listReferenceBases
        for base in method(self._request, self._id):
            print(base.sequence)


class GetReferenceSetRunner(AbstractGetRunner):
    """
    Runner class for the referencesets/{id} method
    """
    def __init__(self, args):
        super(GetReferenceSetRunner, self).__init__(args)

    def run(self):
        self._run(self._httpClient.getReferenceSet)


class GetReferenceRunner(AbstractGetRunner):
    """
    Runner class for the references/{id} method
    """
    def __init__(self, args):
        super(GetReferenceRunner, self).__init__(args)

    def run(self):
        self._run(self._httpClient.getReference)


class BenchmarkRunner(SearchVariantsRunner):
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


def addVariantSearchOptions(parser):
    """
    Adds common options to a variant searches command line parser.
    """
    addVariantSetIdsArgument(parser)
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
    addStartArgument(parser)
    addEndArgument(parser)
    addPageSizeArgument(parser)
    # maxCalls not in protocol; supported by google
    parser.add_argument(
        "--maxCalls", default=1,
        help="The maxiumum number of calls to return")


def addVariantSetIdsArgument(parser):
    parser.add_argument(
        "--variantSetIds", "-V",
        help="The variant set id(s) to search over")


def addStartArgument(parser):
    parser.add_argument(
        "--start", "-s", default=0, type=int,
        help="The start of the search range (inclusive).")


def addEndArgument(parser):
    parser.add_argument(
        "--end", "-e", default=1, type=int,
        help="The end of the search range (exclusive).")


def addIdArgument(parser):
    parser.add_argument("--id", default=None, help="The id of the object")


def addGetArguments(parser):
    addIdArgument(parser)
    addUrlArgument(parser)


def addUrlArgument(parser):
    """
    Adds the URL endpoint argument to the specified parser.
    """
    parser.add_argument("baseUrl", help="The URL of the API endpoint")


def addAccessionsArgument(parser):
    parser.add_argument(
        "--accessions", default=None,
        help="The accessions to search over")


def addMd5ChecksumsArgument(parser):
    parser.add_argument(
        "--md5checksums", default=None,
        help="The md5checksums to search over")


def addPageSizeArgument(parser):
    parser.add_argument(
        "--pageSize", "-m", default=100, type=int,
        help="The maximum number of results returned in one response.")


def addDatasetIdsArgument(parser):
    parser.add_argument(
        "--datasetIds", default=None,
        help="The datasetIds to search over")


def addNameArgument(parser):
    parser.add_argument(
        "--name", default=None,
        help="The name to search over")


def setCommaSeparatedAttribute(request, args, attr):
    attribute = getattr(args, attr)
    if attribute is not None:
        setattr(request, attr, attribute.split(","))


def client_main(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="GA4GH reference client")
    # Add global options
    parser.add_argument('--verbose', '-v', action='count', default=0)
    parser.add_argument(
        "--workarounds", "-w", default='', help="The workarounds to use")
    parser.add_argument(
        "--key", "-k", help="The auth key to use")
    parser.add_argument(
        "--minimalOutput", "-O", default=False,
        help="Use minimal output; default False",
        action='store_true')
    subparsers = parser.add_subparsers(title='subcommands')

    # help
    subparsers.add_parser(
        "help", description="ga4gh_client help",
        help="show this help message and exit")

    # benchmarking
    bmParser = subparsers.add_parser(
        "benchmark",
        description="Run simple benchmarks on the various methods",
        help="Benchmark server performance")
    bmParser.set_defaults(runner=BenchmarkRunner)
    addUrlArgument(bmParser)
    addVariantSearchOptions(bmParser)

    # variants/search
    vsParser = subparsers.add_parser(
        "variants-search",
        description="Search for variants",
        help="Search for variants.")
    vsParser.set_defaults(runner=SearchVariantsRunner)
    addUrlArgument(vsParser)
    addVariantSearchOptions(vsParser)

    # variantsets/search
    vssParser = subparsers.add_parser(
        "variantsets-search",
        description="Search for variantSets",
        help="Search for variantSets.")
    vssParser.set_defaults(runner=SearchVariantSetsRunner)
    addUrlArgument(vssParser)
    addPageSizeArgument(vssParser)
    addDatasetIdsArgument(vssParser)

    # referencesets/search
    rssParser = subparsers.add_parser(
        "referencesets-search",
        description="Search for referenceSets",
        help="Search for referenceSets")
    rssParser.set_defaults(runner=SearchReferenceSetsRunner)
    addUrlArgument(rssParser)
    addPageSizeArgument(rssParser)
    addAccessionsArgument(rssParser)
    addMd5ChecksumsArgument(rssParser)
    rssParser.add_argument(
        "--assemblyId",
        help="The assembly id to search over")

    # references/search
    rsParser = subparsers.add_parser(
        "references-search",
        description="Search for references",
        help="Search for references")
    rsParser.set_defaults(runner=SearchReferencesRunner)
    addUrlArgument(rsParser)
    addPageSizeArgument(rsParser)
    addAccessionsArgument(rsParser)
    addMd5ChecksumsArgument(rsParser)

    # readgroupsets/search
    rgsParser = subparsers.add_parser(
        "readgroupsets-search",
        description="Search for readGroupSets",
        help="Search for readGroupSets")
    rgsParser.set_defaults(runner=SearchReadGroupSetsRunner)
    addUrlArgument(rgsParser)
    addPageSizeArgument(rgsParser)
    addDatasetIdsArgument(rgsParser)
    addNameArgument(rgsParser)

    # callsets/search
    csParser = subparsers.add_parser(
        "callsets-search",
        description="Search for callSets",
        help="Search for callSets")
    csParser.set_defaults(runner=SearchCallSetsRunner)
    addUrlArgument(csParser)
    addPageSizeArgument(csParser)
    addNameArgument(csParser)
    addVariantSetIdsArgument(csParser)

    # reads/search
    rParser = subparsers.add_parser(
        "reads-search",
        description="Search for reads",
        help="Search for reads")
    rParser.set_defaults(runner=SearchReadsRunner)
    addUrlArgument(rParser)
    addPageSizeArgument(rParser)
    addStartArgument(rParser)
    addEndArgument(rParser)
    rParser.add_argument(
        "--readGroupIds", default=None,
        help="The readGroupIds to search over")
    rParser.add_argument(
        "--referenceId", default=None,
        help="The referenceId to search over")
    rParser.add_argument(
        "--referenceName", default=None,
        help="The referenceName to search over")

    # referencesets/{id}
    rsidParser = subparsers.add_parser(
        "referencesets-get",
        description="Get a referenceset",
        help="Get a referenceset")
    rsidParser.set_defaults(runner=GetReferenceSetRunner)
    addGetArguments(rsidParser)

    # references/{id}
    ridParser = subparsers.add_parser(
        "references-get",
        description="Get a reference",
        help="Get a reference")
    ridParser.set_defaults(runner=GetReferenceRunner)
    addGetArguments(ridParser)

    # referencesets/{id}/bases
    basesParser = subparsers.add_parser(
        "references-list-bases",
        description="List bases of a reference",
        help="List bases of a reference")
    basesParser.set_defaults(runner=ListReferenceBasesRunner)
    addGetArguments(basesParser)
    addStartArgument(basesParser)
    addEndArgument(basesParser)

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()
