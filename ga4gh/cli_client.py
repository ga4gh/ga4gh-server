"""
The client command line interface
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import logging

import requests

import ga4gh.backend as backend
import ga4gh.client as client
import ga4gh.exceptions as exceptions
import ga4gh.datarepo as datarepo
import ga4gh.cli_common as cli_common


# the maximum value of a long type in avro = 2**63 - 1
# (64 bit signed integer)
# http://avro.apache.org/docs/1.7.7/spec.html#schema_primitive
# AVRO_LONG_MAX = (1 << 63) - 1
# TODO in the meantime, this is the max value pysam can handle
# This should be removed once pysam input sanitisation has been
# implemented.
AVRO_LONG_MAX = 2**31 - 1


def verbosityToLogLevel(verbosity):
    """
    Returns the specfied verbosity level interpreted as a logging level.
    """
    ret = 0
    if verbosity == 1:
        ret = logging.INFO
    elif verbosity >= 2:
        ret = logging.DEBUG
    return ret


class AbstractQueryRunner(object):
    """
    Abstract base class for runner classes
    """
    def __init__(self, args):
        self._key = args.key
        # TODO this is an experimental addition which is useful for
        # testing. We should think about this and document it if we
        # this it's a useful feature. There is an argument for pushing
        # the backend instantiation into the client, and letting the
        # client be a factory, instantiating the correct Client class
        # depending on the prefix.
        filePrefix = "file://"
        if args.baseUrl.startswith(filePrefix):
            dataDir = args.baseUrl[len(filePrefix):]
            theBackend = backend.Backend(
                datarepo.FileSystemDataRepository(dataDir))
            self._client = client.LocalClient(theBackend)
        else:
            self._client = client.HttpClient(
                args.baseUrl, verbosityToLogLevel(args.verbose), self._key)


class FormattedOutputRunner(AbstractQueryRunner):
    """
    Superclass of runners that support output in common formats.
    """
    def __init__(self, args):
        super(FormattedOutputRunner, self).__init__(args)
        self._output = self._textOutput
        if args.outputFormat == "json":
            self._output = self._jsonOutput

    def _jsonOutput(self, gaObjects):
        """
        Outputs the specified protocol objects as one JSON string per
        line.
        """
        for gaObject in gaObjects:
            print(gaObject.toJsonString())

    def _textOutput(self, gaObjects):
        """
        Outputs a text summary of the specified protocol objects, one
        per line.
        """
        for gaObject in gaObjects:
            print(gaObject.id, gaObject.name, sep="\t")


class AbstractGetRunner(FormattedOutputRunner):
    """
    Abstract base class for get runner classes
    """
    def __init__(self, args):
        super(AbstractGetRunner, self).__init__(args)
        self._id = args.id

    def run(self):
        response = self._method(self._id)
        self._output([response])


class AbstractSearchRunner(FormattedOutputRunner):
    """
    Abstract base class for search runner classes
    """

    def __init__(self, args):
        super(AbstractSearchRunner, self).__init__(args)
        self._pageSize = args.pageSize
        self._client.setPageSize(self._pageSize)

    def getAllDatasets(self):
        """
        Returns all datasets on the server.
        """
        return self._client.searchDatasets()

    def getAllVariantSets(self):
        """
        Returns all variant sets on the server.
        """
        for dataset in self.getAllDatasets():
            iterator = self._client.searchVariantSets(datasetId=dataset.id)
            for variantSet in iterator:
                yield variantSet

    def getAllReadGroupSets(self):
        """
        Returns all readgroup sets on the server.
        """
        for dataset in self.getAllDatasets():
            iterator = self._client.searchReadGroupSets(
                datasetId=dataset.id)
            for readGroupSet in iterator:
                yield readGroupSet

    def getAllReadGroups(self):
        """
        Get all read groups in a read group set
        """
        for dataset in self.getAllDatasets():
            iterator = self._client.searchReadGroupSets(
                datasetId=dataset.id)
            for readGroupSet in iterator:
                readGroupSet = self._client.getReadGroupSet(readGroupSet.id)
                for readGroup in readGroupSet.readGroups:
                    yield readGroup.id

    def getAllReferenceSets(self):
        """
        Returns all reference sets on the server.
        """
        return self._client.searchReferenceSets()


# Runners for the various search methods

class SearchDatasetsRunner(AbstractSearchRunner):
    """
    Runner class for the datasets/search method
    """
    def __init__(self, args):
        super(SearchDatasetsRunner, self).__init__(args)

    def run(self):
        iterator = self._client.searchDatasets()
        self._output(iterator)


class SearchReferenceSetsRunner(AbstractSearchRunner):
    """
    Runner class for the referencesets/search method.
    """
    def __init__(self, args):
        super(SearchReferenceSetsRunner, self).__init__(args)
        self._accession = args.accession
        self._md5checksum = args.md5checksum

    def run(self):
        iterator = self._client.searchReferenceSets(
            accession=self._accession, md5checksum=self._md5checksum)
        self._output(iterator)


class SearchReferencesRunner(AbstractSearchRunner):
    """
    Runner class for the references/search method
    """
    def __init__(self, args):
        super(SearchReferencesRunner, self).__init__(args)
        self._referenceSetId = args.referenceSetId
        self._accession = args.accession
        self._md5checksum = args.md5checksum

    def _run(self, referenceSetId):
        iterator = self._client.searchReferences(
            accession=self._accession, md5checksum=self._md5checksum,
            referenceSetId=referenceSetId)
        self._output(iterator)

    def run(self):
        if self._referenceSetId is None:
            for referenceSet in self.getAllReferenceSets():
                self._run(referenceSet.id)
        else:
            self._run(self._referenceSetId)


class SearchVariantSetsRunner(AbstractSearchRunner):
    """
    Runner class for the variantsets/search method.
    """
    def __init__(self, args):
        super(SearchVariantSetsRunner, self).__init__(args)
        self._datasetId = args.datasetId

    def _run(self, datasetId):
        iterator = self._client.searchVariantSets(datasetId=datasetId)
        self._output(iterator)

    def run(self):
        if self._datasetId is None:
            for dataset in self.getAllDatasets():
                self._run(dataset.id)
        else:
            self._run(self._datasetId)


class SearchReadGroupSetsRunner(AbstractSearchRunner):
    """
    Runner class for the readgroupsets/search method
    """
    def __init__(self, args):
        super(SearchReadGroupSetsRunner, self).__init__(args)
        self._datasetId = args.datasetId
        self._name = args.name

    def _run(self, datasetId):
        iterator = self._client.searchReadGroupSets(
            datasetId=datasetId, name=self._name)
        self._output(iterator)

    def run(self):
        if self._datasetId is None:
            for dataset in self.getAllDatasets():
                self._run(dataset.id)
        else:
            self._run(self._datasetId)


class SearchCallSetsRunner(AbstractSearchRunner):
    """
    Runner class for the callsets/search method
    """
    def __init__(self, args):
        super(SearchCallSetsRunner, self).__init__(args)
        self._variantSetId = args.variantSetId
        self._name = args.name

    def _run(self, variantSetId):
        iterator = self._client.searchCallSets(
            variantSetId=variantSetId, name=self._name)
        self._output(iterator)

    def run(self):
        if self._variantSetId is None:
            for variantSet in self.getAllVariantSets():
                self._run(variantSet.id)
        else:
            self._run(self._variantSetId)


class VariantFormatterMixin(object):
    """
    Simple mixin to format variant objects.
    """
    def _textOutput(self, gaObjects):
        """
        Prints out the specified Variant objects in a VCF-like form.
        """
        for variant in gaObjects:
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


class SearchVariantsRunner(VariantFormatterMixin, AbstractSearchRunner):
    """
    Runner class for the variants/search method.
    """
    def __init__(self, args):
        super(SearchVariantsRunner, self).__init__(args)
        self._referenceName = args.referenceName
        self._variantSetId = args.variantSetId
        self._start = args.start
        self._end = args.end
        if args.callSetIds == []:
            self._callSetIds = []
        elif args.callSetIds == '*':
            self._callSetIds = None
        else:
            self._callSetIds = args.callSetIds.split(",")

    def _run(self, variantSetId):
        iterator = self._client.searchVariants(
            start=self._start, end=self._end,
            referenceName=self._referenceName,
            variantSetId=variantSetId, callSetIds=self._callSetIds)
        self._output(iterator)

    def run(self):
        if self._variantSetId is None:
            for variantSet in self.getAllVariantSets():
                self._run(variantSet.id)
        else:
            self._run(self._variantSetId)


class SearchReadsRunner(AbstractSearchRunner):
    """
    Runner class for the reads/search method
    """
    def __init__(self, args):
        super(SearchReadsRunner, self).__init__(args)
        self._start = args.start
        self._end = args.end
        self._referenceId = args.referenceId
        self._readGroupIds = None
        if args.readGroupIds is not None:
            self._readGroupIds = args.readGroupIds.split(",")

    def _run(self, referenceGroupId, referenceId=None):
        """
        automatically guess reference id if not passed
        """
        # check if we can get reference id from rg
        if referenceId is None:
            referenceId = self._referenceId
        if referenceId is None:
            rg = self._client.getReadGroup(readGroupId=referenceGroupId)
            iterator = self._client.searchReferences(rg.referenceSetId)
            for reference in iterator:
                self._run(referenceGroupId, reference.id)
        else:
            iterator = self._client.searchReads(
                readGroupIds=[referenceGroupId], referenceId=referenceId,
                start=self._start, end=self._end)
            self._output(iterator)

    def run(self):
        """
        Iterate passed read group ids, or go through all available read groups
        """
        if not self._readGroupIds:
            for referenceGroupId in self.getAllReadGroups():
                self._run(referenceGroupId)
        else:
            for referenceGroupId in self._readGroupIds:
                self._run(referenceGroupId)

    def _textOutput(self, gaObjects):
        """
        Prints out the specified Variant objects in a VCF-like form.
        """
        for read in gaObjects:
            # TODO add in some more useful output here.
            print(read.id)


# ListReferenceBases is an oddball, and doesn't fit either get or
# search patterns.

class ListReferenceBasesRunner(AbstractQueryRunner):
    """
    Runner class for the references/{id}/bases method
    """
    def __init__(self, args):
        super(ListReferenceBasesRunner, self).__init__(args)
        self._referenceId = args.id
        self._start = args.start
        self._end = args.end
        self._outputFormat = args.outputFormat

    def run(self):
        sequence = self._client.listReferenceBases(
            self._referenceId, self._start, self._end)
        if self._outputFormat == "text":
            print(sequence)
        else:
            start = self._start if self._start else ""
            end = self._end if self._end else ""
            print(">{}:{}-{}".format(self._referenceId, start, end))

            textWidth = 70
            for index in xrange(0, len(sequence), textWidth):
                print(sequence[index: index+textWidth])


# Runners for the various GET methods.

class GetReferenceSetRunner(AbstractGetRunner):
    """
    Runner class for the referencesets/{id} method
    """
    def __init__(self, args):
        super(GetReferenceSetRunner, self).__init__(args)
        self._method = self._client.getReferenceSet


class GetReferenceRunner(AbstractGetRunner):
    """
    Runner class for the references/{id} method
    """
    def __init__(self, args):
        super(GetReferenceRunner, self).__init__(args)
        self._method = self._client.getReference


class GetReadGroupSetRunner(AbstractGetRunner):
    """
    Runner class for the readgroupsets/{id} method
    """
    def __init__(self, args):
        super(GetReadGroupSetRunner, self).__init__(args)
        self._method = self._client.getReadGroupSet


class GetReadGroupRunner(AbstractGetRunner):
    """
    Runner class for the references/{id} method
    """
    def __init__(self, args):
        super(GetReadGroupRunner, self).__init__(args)
        self._method = self._client.getReadGroup


class GetCallSetRunner(AbstractGetRunner):
    """
    Runner class for the callsets/{id} method
    """
    def __init__(self, args):
        super(GetCallSetRunner, self).__init__(args)
        self._method = self._client.getCallSet


class GetDatasetRunner(AbstractGetRunner):
    """
    Runner class for the datasets/{id} method
    """
    def __init__(self, args):
        super(GetDatasetRunner, self).__init__(args)
        self._method = self._client.getDataset


class GetVariantRunner(VariantFormatterMixin, AbstractGetRunner):
    """
    Runner class for the variants/{id} method
    """
    def __init__(self, args):
        super(GetVariantRunner, self).__init__(args)
        self._method = self._client.getVariant


class GetVariantSetRunner(AbstractGetRunner):
    """
    Runner class for the variantsets/{id} method
    """
    def __init__(self, args):
        super(GetVariantSetRunner, self).__init__(args)
        self._method = self._client.getVariantSet


def addVariantSearchOptions(parser):
    """
    Adds common options to a variant searches command line parser.
    """
    addVariantSetIdArgument(parser)
    addReferenceNameArgument(parser)
    addCallSetIdsArgument(parser)
    addStartArgument(parser)
    addEndArgument(parser)
    addPageSizeArgument(parser)


def addVariantSetIdArgument(parser):
    parser.add_argument(
        "--variantSetId", "-V", default=None,
        help="The variant set id to search over")


def addReferenceNameArgument(parser):
    parser.add_argument(
        "--referenceName", "-r", default="1",
        help="Only return variants on this reference.")


def addCallSetIdsArgument(parser):
    parser.add_argument(
        "--callSetIds", "-c", default=[],
        help="""Return variant calls which belong to call sets
            with these IDs. Pass in IDs as a comma separated list (no spaces).
            Use '*' to request all call sets (the quotes are important!).
            """)


def addStartArgument(parser):
    parser.add_argument(
        "--start", "-s", default=0, type=int,
        help="The start of the search range (inclusive).")


def addEndArgument(parser, defaultValue=AVRO_LONG_MAX):
    parser.add_argument(
        "--end", "-e", default=defaultValue, type=int,
        help="The end of the search range (exclusive).")


def addIdArgument(parser):
    parser.add_argument("id", default=None, help="The id of the object")


def addGetArguments(parser):
    addUrlArgument(parser)
    addIdArgument(parser)
    addOutputFormatArgument(parser)


def addUrlArgument(parser):
    """
    Adds the URL endpoint argument to the specified parser.
    """
    parser.add_argument("baseUrl", help="The URL of the API endpoint")


def addOutputFormatArgument(parser):
    parser.add_argument(
        "--outputFormat", "-O", choices=['text', 'json'], default="text",
        help=(
            "The format for object output. Currently supported are "
            "'text' (default), which gives a short summary of the object and "
            "'json', which outputs each object in line-delimited JSON"))


def addAccessionArgument(parser):
    parser.add_argument(
        "--accession", default=None,
        help="The accession to search for")


def addMd5ChecksumArgument(parser):
    parser.add_argument(
        "--md5checksum", default=None,
        help="The md5checksum to search for")


def addPageSizeArgument(parser):
    parser.add_argument(
        "--pageSize", "-m", default=None, type=int,
        help=(
            "The maximum number of results returned in one page. "
            "The default is to let the server decide how many "
            "results to return in a single page."))


def addDatasetIdArgument(parser):
    parser.add_argument(
        "--datasetId", default=None,
        help="The datasetId to search over")


def addReferenceSetIdArgument(parser):
    parser.add_argument(
        "--referenceSetId", default=None,
        help="The referenceSet to search over")


def addNameArgument(parser):
    parser.add_argument(
        "--name", default=None,
        help="The name to search over")


def addClientGlobalOptions(parser):
    parser.add_argument(
        '--verbose', '-v', action='count', default=0,
        help="Increase verbosity; can be supplied multiple times")
    parser.add_argument(
        "--key", "-k", default='invalid',
        help="Auth Key. Found on server index page.")
    cli_common.addDisableUrllibWarningsArgument(parser)


def addHelpParser(subparsers):
    parser = subparsers.add_parser(
        "help", description="ga4gh_client help",
        help="show this help message and exit")
    return parser


def addVariantsSearchParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "variants-search", "Search for variants")
    parser.set_defaults(runner=SearchVariantsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addVariantSearchOptions(parser)
    return parser


def addVariantSetsSearchParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "variantsets-search", "Search for variantSets")
    parser.set_defaults(runner=SearchVariantSetsRunner)
    addOutputFormatArgument(parser)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addDatasetIdArgument(parser)
    return parser


def addVariantSetsGetParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "variantsets-get", "Get a variantSet")
    parser.set_defaults(runner=GetVariantSetRunner)
    addGetArguments(parser)


def addReferenceSetsSearchParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "referencesets-search", "Search for referenceSets")
    parser.set_defaults(runner=SearchReferenceSetsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addPageSizeArgument(parser)
    addAccessionArgument(parser)
    addMd5ChecksumArgument(parser)
    parser.add_argument(
        "--assemblyId",
        help="The assembly id to search for")
    return parser


def addReferencesSearchParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "references-search", "Search for references")
    parser.set_defaults(runner=SearchReferencesRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addPageSizeArgument(parser)
    addAccessionArgument(parser)
    addMd5ChecksumArgument(parser)
    addReferenceSetIdArgument(parser)
    return parser


def addReadGroupSetsSearchParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "readgroupsets-search", "Search for readGroupSets")
    parser.set_defaults(runner=SearchReadGroupSetsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addPageSizeArgument(parser)
    addDatasetIdArgument(parser)
    addNameArgument(parser)
    return parser


def addCallSetsSearchParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "callsets-search", "Search for callSets")
    parser.set_defaults(runner=SearchCallSetsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addPageSizeArgument(parser)
    addNameArgument(parser)
    addVariantSetIdArgument(parser)
    return parser


def addReadsSearchParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "reads-search", "Search for reads")
    parser.set_defaults(runner=SearchReadsRunner)
    addOutputFormatArgument(parser)
    addReadsSearchParserArguments(parser)
    return parser


def addDatasetsGetParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "datasets-get", "Get a dataset")
    parser.set_defaults(runner=GetDatasetRunner)
    addGetArguments(parser)


def addDatasetsSearchParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "datasets-search", "Search for datasets")
    parser.set_defaults(runner=SearchDatasetsRunner)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addOutputFormatArgument(parser)
    return parser


def addReadsSearchParserArguments(parser):
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addStartArgument(parser)
    addEndArgument(parser)
    parser.add_argument(
        "--readGroupIds", default=None,
        help="The readGroupIds to search over")
    parser.add_argument(
        "--referenceId", default=None,
        help="The referenceId to search over")


def addReferenceSetsGetParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "referencesets-get", "Get a referenceset")
    parser.set_defaults(runner=GetReferenceSetRunner)
    addGetArguments(parser)


def addReferencesGetParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "references-get", "Get a reference")
    parser.set_defaults(runner=GetReferenceRunner)
    addGetArguments(parser)


def addReadGroupSetsGetParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "readgroupsets-get", "Get a read group set")
    parser.set_defaults(runner=GetReadGroupSetRunner)
    addGetArguments(parser)


def addReadGroupsGetParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "readgroups-get", "Get a read group")
    parser.set_defaults(runner=GetReadGroupRunner)
    addGetArguments(parser)


def addCallSetsGetParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "callsets-get", "Get a callSet")
    parser.set_defaults(runner=GetCallSetRunner)
    addGetArguments(parser)


def addVariantsGetParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "variants-get", "Get a variant")
    parser.set_defaults(runner=GetVariantRunner)
    addGetArguments(parser)


def addReferencesBasesListParser(subparsers):
    parser = cli_common.addSubparser(
        subparsers, "references-list-bases", "List bases of a reference")
    parser.add_argument(
        "--outputFormat", "-O", choices=['text', 'fasta'], default="text",
        help=(
            "The format for sequence output. Currently supported are "
            "'text' (default), which prints the sequence out directly and "
            "'fasta', which formats the sequence into fixed width FASTA"))
    parser.set_defaults(runner=ListReferenceBasesRunner)
    addUrlArgument(parser)
    addIdArgument(parser)
    addStartArgument(parser)
    addEndArgument(parser, defaultValue=None)


def getClientParser():
    parser = argparse.ArgumentParser(
        description="GA4GH reference client")
    addClientGlobalOptions(parser)
    subparsers = parser.add_subparsers(title='subcommands',)
    addHelpParser(subparsers)
    addVariantsSearchParser(subparsers)
    addVariantSetsSearchParser(subparsers)
    addVariantSetsGetParser(subparsers)
    addReferenceSetsSearchParser(subparsers)
    addReferencesSearchParser(subparsers)
    addReadGroupSetsSearchParser(subparsers)
    addCallSetsSearchParser(subparsers)
    addReadsSearchParser(subparsers)
    addDatasetsSearchParser(subparsers)
    addReferenceSetsGetParser(subparsers)
    addReferencesGetParser(subparsers)
    addReadGroupSetsGetParser(subparsers)
    addReadGroupsGetParser(subparsers)
    addCallSetsGetParser(subparsers)
    addVariantsGetParser(subparsers)
    addDatasetsGetParser(subparsers)
    addReferencesBasesListParser(subparsers)
    return parser


def client_main(args=None):
    parser = getClientParser()
    parsedArgs = parser.parse_args(args)
    if "runner" not in parsedArgs:
        parser.print_help()
    else:
        if parsedArgs.disable_urllib_warnings:
            requests.packages.urllib3.disable_warnings()
        try:
            runner = parsedArgs.runner(parsedArgs)
            runner.run()
        except (exceptions.BaseClientException,
                requests.exceptions.RequestException) as exception:
            # TODO suppress exception unless debug settings are enabled
            raise exception
