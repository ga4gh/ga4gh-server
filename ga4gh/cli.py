"""
Command line interface programs for the GA4GH reference implementation.

TODO: document how to use these for development and simple deployment.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import logging
import unittest
import unittest.loader
import unittest.suite

import requests

import ga4gh.backend as backend
import ga4gh.client as client
import ga4gh.converters as converters
import ga4gh.frontend as frontend
import ga4gh.configtest as configtest
import ga4gh.exceptions as exceptions
import ga4gh.datarepo as datarepo
import ga4gh.repo_manager as repo_manager
import ga4gh.protocol as protocol


# the maximum value of a long type in avro = 2**63 - 1
# (64 bit signed integer)
# http://avro.apache.org/docs/1.7.7/spec.html#schema_primitive
# AVRO_LONG_MAX = (1 << 63) - 1
# TODO in the meantime, this is the max value pysam can handle
# This should be removed once pysam input sanitisation has been
# implemented.
AVRO_LONG_MAX = 2**31 - 1


##############################################################################
# common
##############################################################################


def addSubparser(subparsers, subcommand, description):
    parser = subparsers.add_parser(
        subcommand, description=description, help=description)
    return parser


##############################################################################
# Server
##############################################################################


def addServerOptions(parser):
    parser.add_argument(
        "--port", "-P", default=8000, type=int,
        help="The port to listen on")
    parser.add_argument(
        "--host", "-H", default="127.0.0.1",
        help="The server host string; use 0.0.0.0 to allow all connections.")
    parser.add_argument(
        "--config", "-c", default='DevelopmentConfig', type=str,
        help="The configuration to use")
    parser.add_argument(
        "--config-file", "-f", type=str, default=None,
        help="The configuration file to use")
    parser.add_argument(
        "--tls", "-t", action="store_true", default=False,
        help="Start in TLS (https) mode.")
    parser.add_argument(
        "--dont-use-reloader", default=False, action="store_true",
        help="Don't use the flask reloader")
    addDisableUrllibWarningsArgument(parser)


def getServerParser():
    parser = argparse.ArgumentParser(
        description="GA4GH reference server")
    addServerOptions(parser)
    return parser


def server_main(args=None):
    parser = getServerParser()
    parsedArgs = parser.parse_args(args)
    if parsedArgs.disable_urllib_warnings:
        requests.packages.urllib3.disable_warnings()
    frontend.configure(
        parsedArgs.config_file, parsedArgs.config, parsedArgs.port)
    sslContext = None
    if parsedArgs.tls or ("OIDC_PROVIDER" in frontend.app.config):
        sslContext = "adhoc"
    frontend.app.run(
        host=parsedArgs.host, port=parsedArgs.port,
        use_reloader=not parsedArgs.dont_use_reloader,
        ssl_context=sslContext)


##############################################################################
# Client
##############################################################################


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
            if hasattr(gaObject, 'name'):
                print(gaObject.id, gaObject.name, sep="\t")
            else:
                print(gaObject.id, sep="\t")


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


class SearchVariantAnnotationSetsRunner(AbstractSearchRunner):
    """
    Runner class for the variantannotationsets/search method.
    """
    def __init__(self, args):
        super(SearchVariantAnnotationSetsRunner, self).__init__(args)
        self._variantSetId = args.variantSetId

    def _run(self, variantSetId):
        iterator = self._client.searchVariantAnnotationSets(
            variantSetId=variantSetId)
        self._output(iterator)

    def run(self):
        self._run(self._variantSetId)


class SearchReadGroupSetsRunner(AbstractSearchRunner):
    """
    Runner class for the readgroupsets/search method
    """
    def __init__(self, args):
        super(SearchReadGroupSetsRunner, self).__init__(args)
        self._datasetId = args.datasetId
        self._name = args.name
        self._bioSampleId = args.bioSampleId

    def _run(self, datasetId):
        iterator = self._client.searchReadGroupSets(
            datasetId=datasetId,
            name=self._name,
            bioSampleId=self._bioSampleId)
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
        self._bioSampleId = args.bioSampleId

    def _run(self, variantSetId):
        iterator = self._client.searchCallSets(
            variantSetId=variantSetId,
            name=self._name,
            bioSampleId=self._bioSampleId)
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


class AnnotationFormatterMixin(object):
    """
    Simple mixin to format variant objects.
    """
    def _textOutput(self, gaObjects):
        """
        Prints out the specified Variant objects in a VCF-like form.
        """
        for variantAnnotation in gaObjects:
            print(
                variantAnnotation.id, variantAnnotation.variantId,
                variantAnnotation.variantAnnotationSetId,
                variantAnnotation.createDateTime, sep="\t", end="\t")
            for effect in variantAnnotation.transcriptEffects:
                print(effect.alternateBases, sep="|", end="|")
                for so in effect.effects:
                    print(so.term, sep="&", end="|")
                    print(so.id, sep="&", end="|")
                print(effect.hgvsAnnotation.transcript,
                      effect.hgvsAnnotation.protein, sep="|", end="\t")
            print()


class SearchBioSamplesRunner(AbstractSearchRunner):
    def __init__(self, args):
        super(SearchBioSamplesRunner, self).__init__(args)
        self._datasetId = args.datasetId
        self._name = args.name

    def _run(self, datasetId):
        iterator = self._client.searchBioSamples(
            datasetId=datasetId,
            name=self._name)
        self._output(iterator)

    def run(self):
        if self._datasetId is None:
            for dataset in self.getAllDatasets():
                self._run(dataset.id)
        else:
            self._run(self._datasetId)


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


class SearchVariantAnnotationsRunner(
        AnnotationFormatterMixin, AbstractSearchRunner):
    """
    Runner class for the variantannotations/search method.
    """
    def __init__(self, args):
        super(SearchVariantAnnotationsRunner, self).__init__(args)
        self._referenceName = args.referenceName
        self._referenceId = args.referenceId
        self._variantAnnotationSetId = args.variantAnnotationSetId
        self._start = args.start
        self._end = args.end

        if args.effects == "":
            self._effects = []
        else:
            self._effects = []
            for eff in args.effects.split(","):
                term = protocol.OntologyTerm()
                term.id = eff
                self._effects.append(term)

    def _run(self, variantAnnotationSetId):
        iterator = self._client.searchVariantAnnotations(
            variantAnnotationSetId=variantAnnotationSetId,
            referenceName=self._referenceName, referenceId=self._referenceId,
            start=self._start, end=self._end,
            effects=self._effects)
        self._output(iterator)

    def getAllAnnotationSets(self):
        """
        Returns all variant annotation sets on the server.
        """
        for dataset in self.getAllDatasets():
            iterator = self._client.searchVariantAnnotationSets(
                datasetId=dataset.id)
            for variantAnnotationSet in iterator:
                yield variantAnnotationSet

    def run(self):
        if self._variantAnnotationSetId is None:
            for annotationSet in self.getAllAnnotationSets():
                self._run(annotationSet.id)
        else:
            self._run(self._variantAnnotationSetId)


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


class GetBioSamplesRunner(AbstractGetRunner):
    """
    Runner class for the biosamples/{id} method
    """
    def __init__(self, args):
        super(GetBioSamplesRunner, self).__init__(args)
        self._method = self._client.getBioSample


def addDisableUrllibWarningsArgument(parser):
    parser.add_argument(
        "--disable-urllib-warnings", default=False, action="store_true",
        help="Disable urllib3 warnings")


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


def addAnnotationsSearchOptions(parser):
    """
    Adds common options to a annotation searches command line parser.
    """
    addAnnotationSetIdArgument(parser)
    addReferenceNameArgument(parser)
    addReferenceIdArgument(parser)
    addStartArgument(parser)
    addEndArgument(parser)
    addFeatureIdsArgument(parser)
    addEffectsArgument(parser)
    addPageSizeArgument(parser)


def addVariantSetIdArgument(parser):
    parser.add_argument(
        "--variantSetId", "-V", default=None,
        help="The variant set id to search over")


def addVariantSetIdMandatoryArgument(parser):
    parser.add_argument(
        "variantSetId", help="The variant set id to search over")


def addAnnotationSetIdArgument(parser):
    parser.add_argument(
        "--variantAnnotationSetId", "-V", default=None,
        help="The annotation set id to search over")


def addReferenceNameArgument(parser):
    parser.add_argument(
        "--referenceName", "-r", default="1",
        help="Only return variants on this reference.")


def addReferenceIdArgument(parser):
    parser.add_argument(
        "--referenceId", "-c", default="",
        help="Only return variants on this reference ID.")


def addCallSetIdsArgument(parser):
    parser.add_argument(
        "--callSetIds", "-c", default=[],
        help="""Return variant calls which belong to call sets
            with these IDs. Pass in IDs as a comma separated list (no spaces).
            Use '*' to request all call sets (the quotes are important!).
            """)


def addFeatureIdsArgument(parser):
    parser.add_argument(
        "--featureIds", "-f", default=[],
        help="""Return annotations on any of the feature IDs.
            Pass in IDs as a comma separated list (no spaces).
            """)


def addEffectsArgument(parser):
    parser.add_argument(
        "--effects", "-effs", default="",
        help="""Return annotations having any of these effects.
            Pass in IDs as a comma separated list (no spaces).
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


def addBioSampleIdArgument(subparser):
    subparser.add_argument(
        "--bioSampleId", default=None,
        help="the id of a BioSample")


def addClientGlobalOptions(parser):
    parser.add_argument(
        '--verbose', '-v', action='count', default=0,
        help="Increase verbosity; can be supplied multiple times")
    parser.add_argument(
        "--key", "-k", default='invalid',
        help="Auth Key. Found on server index page.")
    addDisableUrllibWarningsArgument(parser)


def addHelpParser(subparsers):
    parser = subparsers.add_parser(
        "help", description="ga4gh_client help",
        help="show this help message and exit")
    return parser


def addVariantsSearchParser(subparsers):
    parser = addSubparser(
        subparsers, "variants-search", "Search for variants")
    parser.set_defaults(runner=SearchVariantsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addVariantSearchOptions(parser)
    return parser


def addVariantSetsSearchParser(subparsers):
    parser = addSubparser(
        subparsers, "variantsets-search", "Search for variantSets")
    parser.set_defaults(runner=SearchVariantSetsRunner)
    addOutputFormatArgument(parser)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addDatasetIdArgument(parser)
    return parser


def addVariantAnnotationSearchParser(subparsers):
    parser = subparsers.add_parser(
        "variantannotations-search",
        description="Search for variant annotaions",
        help="Search for variantAnnotaions.")
    parser.set_defaults(runner=SearchVariantAnnotationsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addAnnotationsSearchOptions(parser)
    return parser


def addVariantAnnotationSetsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "variantannotationsets-search",
        description="Search for variant annotation sets",
        help="Search for variantAnnotationSets.")
    parser.set_defaults(runner=SearchVariantAnnotationSetsRunner)
    addOutputFormatArgument(parser)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addVariantSetIdMandatoryArgument(parser)
    return parser


def addVariantSetsGetParser(subparsers):
    parser = addSubparser(
        subparsers, "variantsets-get", "Get a variantSet")
    parser.set_defaults(runner=GetVariantSetRunner)
    addGetArguments(parser)


def addReferenceSetsSearchParser(subparsers):
    parser = addSubparser(
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
    parser = addSubparser(
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
    parser = addSubparser(
        subparsers, "readgroupsets-search", "Search for readGroupSets")
    parser.set_defaults(runner=SearchReadGroupSetsRunner)
    addUrlArgument(parser)
    addBioSampleIdArgument(parser)
    addOutputFormatArgument(parser)
    addPageSizeArgument(parser)
    addDatasetIdArgument(parser)
    addNameArgument(parser)
    return parser


def addCallSetsSearchParser(subparsers):
    parser = addSubparser(
        subparsers, "callsets-search", "Search for callSets")
    parser.set_defaults(runner=SearchCallSetsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addPageSizeArgument(parser)
    addNameArgument(parser)
    addBioSampleIdArgument(parser)
    addVariantSetIdArgument(parser)
    return parser


def addReadsSearchParser(subparsers):
    parser = addSubparser(
        subparsers, "reads-search", "Search for reads")
    parser.set_defaults(runner=SearchReadsRunner)
    addOutputFormatArgument(parser)
    addReadsSearchParserArguments(parser)
    return parser


def addBioSampleSearchParser(subparsers):
    parser = addSubparser(
        subparsers, "biosamples-search", "Search for BioSamples")
    parser.set_defaults(runner=SearchBioSamplesRunner)
    addUrlArgument(parser)
    addDatasetIdArgument(parser)
    addOutputFormatArgument(parser)
    addNameArgument(parser)
    addPageSizeArgument(parser)
    return parser


def addBioSampleGetParser(subparsers):
    parser = addSubparser(
        subparsers, "biosamples-get", "Get a BioSamples")
    parser.set_defaults(runner=GetBioSamplesRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addIdArgument(parser)
    addPageSizeArgument(parser)
    return parser


def addDatasetsGetParser(subparsers):
    parser = addSubparser(
        subparsers, "datasets-get", "Get a dataset")
    parser.set_defaults(runner=GetDatasetRunner)
    addGetArguments(parser)


def addDatasetsSearchParser(subparsers):
    parser = addSubparser(
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
    parser = addSubparser(
        subparsers, "referencesets-get", "Get a referenceset")
    parser.set_defaults(runner=GetReferenceSetRunner)
    addGetArguments(parser)


def addReferencesGetParser(subparsers):
    parser = addSubparser(
        subparsers, "references-get", "Get a reference")
    parser.set_defaults(runner=GetReferenceRunner)
    addGetArguments(parser)


def addReadGroupSetsGetParser(subparsers):
    parser = addSubparser(
        subparsers, "readgroupsets-get", "Get a read group set")
    parser.set_defaults(runner=GetReadGroupSetRunner)
    addGetArguments(parser)


def addReadGroupsGetParser(subparsers):
    parser = addSubparser(
        subparsers, "readgroups-get", "Get a read group")
    parser.set_defaults(runner=GetReadGroupRunner)
    addGetArguments(parser)


def addCallSetsGetParser(subparsers):
    parser = addSubparser(
        subparsers, "callsets-get", "Get a callSet")
    parser.set_defaults(runner=GetCallSetRunner)
    addGetArguments(parser)


def addVariantsGetParser(subparsers):
    parser = addSubparser(
        subparsers, "variants-get", "Get a variant")
    parser.set_defaults(runner=GetVariantRunner)
    addGetArguments(parser)


def addReferencesBasesListParser(subparsers):
    parser = addSubparser(
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
    addVariantAnnotationSearchParser(subparsers)
    addVariantAnnotationSetsSearchParser(subparsers)
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
    addBioSampleSearchParser(subparsers)
    addBioSampleGetParser(subparsers)
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


##############################################################################
# ga2vcf
##############################################################################

class Ga2VcfRunner(SearchVariantsRunner):
    """
    Runner class for the ga2vcf
    """
    def __init__(self, args):
        super(Ga2VcfRunner, self).__init__(args)
        self._outputFile = args.outputFile
        self._binaryOutput = False
        if args.outputFormat == "bcf":
            self._binaryOutput = True

    def run(self):
        variantSet = self._client.getVariantSet(self._variantSetId)
        iterator = self._client.searchVariants(
            start=self._start, end=self._end,
            referenceName=self._referenceName,
            variantSetId=self._variantSetId,
            callSetIds=self._callSetIds)
        # do conversion
        vcfConverter = converters.VcfConverter(
            variantSet, iterator, self._outputFile, self._binaryOutput)
        vcfConverter.convert()


def addOutputFileArgument(parser):
    parser.add_argument(
        "--outputFile", "-o", default=None,
        help="the file to write the output to")


def getGa2VcfParser():
    parser = argparse.ArgumentParser(
        description=(
            "GA4GH VCF conversion tool. Converts variant information "
            "stored in a GA4GH repository into VCF format."))
    addClientGlobalOptions(parser)
    addOutputFileArgument(parser)
    addUrlArgument(parser)
    parser.add_argument("variantSetId", help="The variant set to convert")
    parser.add_argument(
        "--outputFormat", "-O", choices=['vcf', 'bcf'], default="vcf",
        help=(
            "The format for object output. Currently supported are "
            "'vcf' (default), which is a text-based format and "
            "'bcf', which is the binary equivalent"))
    addReferenceNameArgument(parser)
    addCallSetIdsArgument(parser)
    addStartArgument(parser)
    addEndArgument(parser)
    addPageSizeArgument(parser)
    return parser


def ga2vcf_main():
    parser = getGa2VcfParser()
    args = parser.parse_args()
    if "baseUrl" not in args:
        parser.print_help()
    else:
        runner = Ga2VcfRunner(args)
        runner.run()


##############################################################################
# ga2sam
##############################################################################


class Ga2SamRunner(SearchReadsRunner):
    """
    Runner class for the ga2vcf
    """
    def __init__(self, args):
        args.readGroupIds = args.readGroupId
        super(Ga2SamRunner, self).__init__(args)
        self._outputFile = args.outputFile
        self._binaryOutput = False
        if args.outputFormat == "bam":
            self._binaryOutput = True

    def run(self):
        samConverter = converters.SamConverter(
            self._client, readGroupId=self._readGroupIds[0],
            referenceId=self._referenceId, start=self._start, end=self._end,
            outputFileName=self._outputFile, binaryOutput=self._binaryOutput)
        samConverter.convert()


def getGa2SamParser():
    parser = argparse.ArgumentParser(
        description="GA4GH SAM conversion tool")
    addClientGlobalOptions(parser)
    addUrlArgument(parser)
    parser.add_argument(
        "readGroupId",
        help="The ReadGroup to convert to SAM/BAM format.")
    addPageSizeArgument(parser)
    addStartArgument(parser)
    addEndArgument(parser)
    parser.add_argument(
        "--referenceId", default=None,
        help="The referenceId to search over")
    parser.add_argument(
        "--outputFormat", "-O", default="sam", choices=["sam", "bam"],
        help=(
            "The format for object output. Currently supported are "
            "'sam' (default), which is a text-based format and "
            "'bam', which is the binary equivalent"))
    addOutputFileArgument(parser)
    return parser


def ga2sam_main():
    parser = getGa2SamParser()
    args = parser.parse_args()
    if "baseUrl" not in args:
        parser.print_help()
    else:
        runner = Ga2SamRunner(args)
        runner.run()


##############################################################################
# Configuration testing
##############################################################################


class SimplerResult(unittest.TestResult):
    """
    The TestResult class gives formatted tracebacks as error messages, which
    is not what we want. Instead we just want the error message from the
    err praram. Hence this subclass.
    """
    def addError(self, test, err):
        self.errors.append((test,
                            "{0}: {1}".format(err[0].__name__, err[1])))

    def addFailure(self, test, err):
        self.failures.append((test,
                              "{0}: {1}".format(err[0].__name__, err[1])))


def configtest_main(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="GA4GH server configuration validator")
    parser.add_argument(
        "--config", "-c", default='DevelopmentConfig', type=str,
        help="The configuration to use")
    parser.add_argument(
        "--config-file", "-f", type=str, default=None,
        help="The configuration file to use")
    args = parser.parse_args()
    configStr = 'ga4gh.serverconfig:{0}'.format(args.config)

    configtest.TestConfig.configStr = configStr
    configtest.TestConfig.configFile = args.config_file
    configtest.TestConfig.configEnv = "GA4GH_CONFIGURATION"

    loader = unittest.TestLoader()
    tests = loader.loadTestsFromModule(configtest)
    results = SimplerResult()
    tests.run(results)

    logging.basicConfig(level=logging.INFO)
    log = logging.getLogger(__name__)
    log.info('{0} Tests run. {1} errors, {2} failures, {3} skipped'.
             format(results.testsRun,
                    len(results.errors),
                    len(results.failures),
                    len(results.skipped)))
    for result in results.errors:
        if result is not None:
            log.critical('Error: {0}: {1}'.format(result[0].id(), result[1]))
    for result in results.failures:
        if result is not None:
            log.critical('Failure: {0}: {1}'.format(result[0].id(), result[1]))
    for result in results.skipped:
        if result is not None:
            log.info('Skipped: {0}: {1}'.format(result[0].id(), result[1]))

##############################################################################
# data repository management tool
##############################################################################


def getRawInput(displayString):
    userResponse = raw_input(displayString)
    return userResponse


class AbstractRepoCommandRunner(object):

    def __init__(self, args):
        self.args = args
        self.repoPath = args.repoPath
        self.repoManager = repo_manager.RepoManager(self.repoPath)
        if 'force' in args:
            self.force = args.force

    def confirmRun(self, func, deleteString):
        if self.force:
            func()
        else:
            displayString = (
                "Are you sure you want to delete data in {}? "
                "[y|N] ".format(deleteString))
            userResponse = getRawInput(displayString)
            if userResponse.strip() == 'y':
                func()
            else:
                print("Aborted")


class AbstractRepoAddCommandRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(AbstractRepoAddCommandRunner, self).__init__(args)
        self.filePath = args.filePath


class AbstractRepoAddMoveCommandRunner(AbstractRepoAddCommandRunner):

    def __init__(self, args):
        super(AbstractRepoAddMoveCommandRunner, self).__init__(args)
        self.moveMode = args.moveMode


class AbstractRepoDatasetCommandRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(AbstractRepoDatasetCommandRunner, self).__init__(args)
        self.datasetName = args.datasetName


class AbstractRepoDatasetFilepathCommandRunner(
        AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(AbstractRepoDatasetFilepathCommandRunner, self).__init__(args)
        self.filePath = args.filePath
        self.moveMode = args.moveMode


class CheckRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(CheckRunner, self).__init__(args)
        self.doConsistencyCheck = not args.skipConsistencyCheck

    def run(self):
        self.repoManager.check(self.doConsistencyCheck)


class ListRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(ListRunner, self).__init__(args)

    def run(self):
        self.repoManager.list()


class DestroyRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(DestroyRunner, self).__init__(args)

    def run(self):
        def func(): self.repoManager.destroy()
        self.confirmRun(func, 'the Data Repository')


class InitRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(InitRunner, self).__init__(args)

    def run(self):
        self.repoManager.init()


class AddDatasetRunner(AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(AddDatasetRunner, self).__init__(args)

    def run(self):
        self.repoManager.addDataset(self.datasetName)


class RemoveDatasetRunner(AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(RemoveDatasetRunner, self).__init__(args)

    def run(self):
        def func(): self.repoManager.removeDataset(self.datasetName)
        self.confirmRun(func, 'dataset {}'.format(self.datasetName))


class AddReferenceSetRunner(AbstractRepoAddMoveCommandRunner):

    def __init__(self, args):
        super(AddReferenceSetRunner, self).__init__(args)
        self.metadata = {
            'description': args.description,
        }

    def run(self):
        self.repoManager.addReferenceSet(
            self.filePath, self.moveMode, self.metadata)


class RemoveReferenceSetRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(RemoveReferenceSetRunner, self).__init__(args)
        self.referenceSetName = args.referenceSetName

    def run(self):
        def func():
            self.repoManager.removeReferenceSet(self.referenceSetName)
        self.confirmRun(
            func, 'reference set {}'.format(self.referenceSetName))


class AddReadGroupSetRunner(AbstractRepoDatasetFilepathCommandRunner):

    def __init__(self, args):
        super(AddReadGroupSetRunner, self).__init__(args)

    def run(self):
        self.repoManager.addReadGroupSet(
            self.datasetName, self.filePath, self.moveMode)


class RemoveReadGroupSetRunner(AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(RemoveReadGroupSetRunner, self).__init__(args)
        self.readGroupSetName = args.readGroupSetName

    def run(self):
        def func():
            self.repoManager.removeReadGroupSet(
                self.datasetName, self.readGroupSetName)
        self.confirmRun(
            func, 'read group set {}'.format(self.readGroupSetName))


class AddVariantSetRunner(AbstractRepoDatasetFilepathCommandRunner):

    def __init__(self, args):
        super(AddVariantSetRunner, self).__init__(args)

    def run(self):
        self.repoManager.addVariantSet(
            self.datasetName, self.filePath, self.moveMode)


class RemoveVariantSetRunner(AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(RemoveVariantSetRunner, self).__init__(args)
        self.variantSetName = args.variantSetName

    def run(self):
        def func():
            self.repoManager.removeVariantSet(
                self.datasetName, self.variantSetName)
        self.confirmRun(func, 'variant set {}'.format(self.variantSetName))


class AddBioSampleRunner(AbstractRepoDatasetFilepathCommandRunner):

    def __init__(self, args):
        super(AddBioSampleRunner, self).__init__(args)

    def run(self):
        self.repoManager.addBioSample(
            self.datasetName, self.filePath, self.moveMode)


class RemoveBioSampleRunner(AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(RemoveBioSampleRunner, self).__init__(args)
        self.bioSampleName = args.bioSampleName

    def run(self):
        def func():
            self.repoManager.removeBioSample(
                self.datasetName, self.bioSampleName)
        self.confirmRun(
            func, 'BioSample {}'.format(self.bioSampleName))


def addRepoArgument(subparser):
    subparser.add_argument(
        "repoPath", help="the file path of the data repository")


def addSkipConsistencyCheckArgument(subparser):
    subparser.add_argument(
        "-s", "--skipConsistencyCheck", action='store_true', default=False,
        help="skip the data repo consistency check")


def addForceArgument(subparser):
    subparser.add_argument(
        "-f", "--force", action='store_true',
        default=False, help="do not prompt for confirmation")


def addDatasetNameArgument(subparser):
    subparser.add_argument(
        "datasetName", help="the name of the dataset to create/modify")


def addBioSampleNameArgument(subparser):
    subparser.add_argument(
        "bioSampleName",
        help="the name of the BioSample")


def addReadGroupSetNameArgument(subparser):
    subparser.add_argument(
        "readGroupSetName",
        help="the name of the read group set")


def addVariantSetNameArgument(subparser):
    subparser.add_argument(
        "variantSetName",
        help="the name of the variant set")


def addFilePathArgument(subparser):
    subparser.add_argument(
        "filePath", help="the path of the file to be moved into the repo")


def addMoveModeArgument(subparser):
    subparser.add_argument(
        "--moveMode",
        help="move, copy or link the target file (default link)",
        choices=['move', 'copy', 'link'],
        default='link')


def addReferenceSetMetadataArguments(subparser):
    subparser.add_argument(
        "--description",
        help="description of the reference set",
        default="TODO")


def getRepoParser():
    parser = argparse.ArgumentParser(
        description="GA4GH data repository management tool")
    subparsers = parser.add_subparsers(title='subcommands',)
    parser.add_argument(
        "--loud", default=False, action="store_true",
        help="propagate exceptions from RepoManager")

    initParser = addSubparser(
        subparsers, "init", "Initialize a data repository")
    initParser.set_defaults(runner=InitRunner)
    addRepoArgument(initParser)

    checkParser = addSubparser(
        subparsers, "check", "Check to see if repo is well-formed")
    checkParser.set_defaults(runner=CheckRunner)
    addRepoArgument(checkParser)
    addSkipConsistencyCheckArgument(checkParser)

    listParser = addSubparser(
        subparsers, "list", "List the contents of the repo")
    listParser.set_defaults(runner=ListRunner)
    addRepoArgument(listParser)

    destroyParser = addSubparser(
        subparsers, "destroy", "Destroy the repo")
    destroyParser.set_defaults(runner=DestroyRunner)
    addRepoArgument(destroyParser)
    addForceArgument(destroyParser)

    addDatasetParser = addSubparser(
        subparsers, "add-dataset", "Add a dataset to the data repo")
    addDatasetParser.set_defaults(runner=AddDatasetRunner)
    addRepoArgument(addDatasetParser)
    addDatasetNameArgument(addDatasetParser)

    removeDatasetParser = addSubparser(
        subparsers, "remove-dataset",
        "Remove a dataset from the data repo")
    removeDatasetParser.set_defaults(runner=RemoveDatasetRunner)
    addRepoArgument(removeDatasetParser)
    addDatasetNameArgument(removeDatasetParser)
    addForceArgument(removeDatasetParser)

    addReferenceSetParser = addSubparser(
        subparsers, "add-referenceset",
        "Add a reference set to the data repo")
    addReferenceSetParser.set_defaults(runner=AddReferenceSetRunner)
    addRepoArgument(addReferenceSetParser)
    addFilePathArgument(addReferenceSetParser)
    addMoveModeArgument(addReferenceSetParser)
    addReferenceSetMetadataArguments(addReferenceSetParser)

    removeReferenceSetParser = addSubparser(
        subparsers, "remove-referenceset",
        "Remove a reference set from the repo")
    removeReferenceSetParser.set_defaults(runner=RemoveReferenceSetRunner)
    addRepoArgument(removeReferenceSetParser)
    removeReferenceSetParser.add_argument(
        "referenceSetName",
        help="the name of the reference set")
    addForceArgument(removeReferenceSetParser)

    addReadGroupSetParser = addSubparser(
        subparsers, "add-readgroupset",
        "Add a read group set to the data repo")
    addReadGroupSetParser.set_defaults(runner=AddReadGroupSetRunner)
    addRepoArgument(addReadGroupSetParser)
    addDatasetNameArgument(addReadGroupSetParser)
    addFilePathArgument(addReadGroupSetParser)
    addMoveModeArgument(addReadGroupSetParser)

    removeReadGroupSetParser = addSubparser(
        subparsers, "remove-readgroupset",
        "Remove a read group set from the repo")
    removeReadGroupSetParser.set_defaults(runner=RemoveReadGroupSetRunner)
    addRepoArgument(removeReadGroupSetParser)
    addDatasetNameArgument(removeReadGroupSetParser)
    addReadGroupSetNameArgument(removeReadGroupSetParser)
    addForceArgument(removeReadGroupSetParser)

    addVariantSetParser = addSubparser(
        subparsers, "add-variantset", "Add a variant set to the data repo")
    addVariantSetParser.set_defaults(runner=AddVariantSetRunner)
    addRepoArgument(addVariantSetParser)
    addDatasetNameArgument(addVariantSetParser)
    addFilePathArgument(addVariantSetParser)
    addMoveModeArgument(addVariantSetParser)

    removeVariantSetParser = addSubparser(
        subparsers, "remove-variantset",
        "Remove a variant set from the repo")
    removeVariantSetParser.set_defaults(runner=RemoveVariantSetRunner)
    addRepoArgument(removeVariantSetParser)
    addDatasetNameArgument(removeVariantSetParser)
    addVariantSetNameArgument(removeVariantSetParser)
    addForceArgument(removeVariantSetParser)

    removeBioSampleParser = addSubparser(
        subparsers, "remove-biosample",
        "Remove a biosample from the repo")
    removeBioSampleParser.set_defaults(runner=RemoveBioSampleRunner)
    addRepoArgument(removeBioSampleParser)
    addDatasetNameArgument(removeBioSampleParser)
    addBioSampleNameArgument(removeBioSampleParser)
    addForceArgument(removeBioSampleParser)

    addBioSampleParser = addSubparser(
        subparsers, "add-biosample",
        "Add a biosample to the repo")
    addBioSampleParser.set_defaults(runner=AddBioSampleRunner)
    addRepoArgument(addBioSampleParser)
    addDatasetNameArgument(addBioSampleParser)
    addFilePathArgument(addBioSampleParser)
    addMoveModeArgument(addBioSampleParser)
    return parser

    return parser


def repo_main(args=None):
    parser = getRepoParser()
    parsedArgs = parser.parse_args(args)
    if "runner" not in parsedArgs:
        parser.print_help()
    else:
        runner = parsedArgs.runner(parsedArgs)
        try:
            runner.run()
        except exceptions.RepoManagerException as exception:
            print(exception.message)
            if parsedArgs.loud:
                raise
