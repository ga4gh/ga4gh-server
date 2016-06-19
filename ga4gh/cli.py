"""
Command line interface programs for the GA4GH reference implementation.

TODO: document how to use these for development and simple deployment.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import glob
import logging
import operator
import os
import sys
import textwrap
import traceback
import unittest
import unittest.loader
import unittest.suite
import urlparse

import requests

import ga4gh
import ga4gh.backend as backend
import ga4gh.client as client
import ga4gh.converters as converters
import ga4gh.frontend as frontend
import ga4gh.configtest as configtest
import ga4gh.exceptions as exceptions
import ga4gh.datarepo as datarepo
import ga4gh.protocol as protocol
import ga4gh.datamodel.reads as reads
import ga4gh.datamodel.variants as variants
import ga4gh.datamodel.references as references
import ga4gh.datamodel.sequenceAnnotations as sequenceAnnotations
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.ontologies as ontologies


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


class SortedHelpFormatter(argparse.HelpFormatter):
    """
    An argparse HelpFormatter that sorts the flags and subcommands
    in alphabetical order
    """
    def add_arguments(self, actions):
        """
        Sort the flags alphabetically
        """
        actions = sorted(
            actions, key=operator.attrgetter('option_strings'))
        super(SortedHelpFormatter, self).add_arguments(actions)

    def _iter_indented_subactions(self, action):
        """
        Sort the subcommands alphabetically
        """
        try:
            get_subactions = action._get_subactions
        except AttributeError:
            pass
        else:
            self._indent()
            if isinstance(action, argparse._SubParsersAction):
                for subaction in sorted(
                        get_subactions(), key=lambda x: x.dest):
                    yield subaction
            else:
                for subaction in get_subactions():
                    yield subaction
            self._dedent()


def addSubparser(subparsers, subcommand, description):
    parser = subparsers.add_parser(
        subcommand, description=description, help=description)
    return parser


def createArgumentParser(description):
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=SortedHelpFormatter)
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
    addVersionArgument(parser)
    addDisableUrllibWarningsArgument(parser)


def getServerParser():
    parser = createArgumentParser("GA4GH reference server")
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
            registryPath = args.baseUrl[len(filePrefix):]
            repo = datarepo.SqlDataRepository(registryPath)
            repo.open(datarepo.MODE_READ)
            theBackend = backend.Backend(repo)
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
            print(protocol.toJson(gaObject))

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

    def getAllFeatureSets(self):
        """
        Returns all feature sets on the server.
        """
        for dataset in self.getAllDatasets():
            iterator = self._client.searchFeatureSets(datasetId=dataset.id)
            for featureSet in iterator:
                yield featureSet

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
                for readGroup in readGroupSet.read_groups:
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


class SearchFeatureSetsRunner(AbstractSearchRunner):
    """
    Runner class for the featuresets/search method.
    """
    def __init__(self, args):
        super(SearchFeatureSetsRunner, self).__init__(args)
        self._datasetId = args.datasetId

    def _run(self, datasetId):
        iterator = self._client.searchFeatureSets(datasetId=datasetId)
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
                variant.id, variant.variant_set_id, variant.names,
                variant.reference_name, variant.start, variant.end,
                variant.reference_bases, variant.alternate_bases,
                sep="\t", end="\t")
            for key, value in variant.info.items():
                print(key, value, sep="=", end=";")
            print("\t", end="")
            for c in variant.calls:
                print(
                    c.call_set_id, c.genotype, c.genotype_likelihood, c.info,
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
                variantAnnotation.id, variantAnnotation.variant_id,
                variantAnnotation.variant_annotation_set_id,
                variantAnnotation.create_date_time, sep="\t", end="\t")
            for effect in variantAnnotation.transcript_effects:
                print(effect.alternate_bases, sep="|", end="|")
                for so in effect.effects:
                    print(so.term, sep="&", end="|")
                    print(so.id, sep="&", end="|")
                print(effect.hgvs_annotation.transcript,
                      effect.hgvs_annotation.protein, sep="|", end="\t")
            print()


class FeatureFormatterMixin(object):
    """
    Mix-in class to format Feature (Sequence Annotation) objects
    """
    def _textOutput(self, gaObjects):
        for feature in gaObjects:
            print(
                feature.id, feature.parentId, feature.featureSetId,
                feature.referenceName, feature.start, feature.end,
                feature.strand, sep="\t", end="\t")
            print(
                "FeatureType:", feature.featureType.id,
                feature.featureType.term, end="\t")
            for attrkey in feature.attributes.vals.keys():
                print(
                    attrkey, feature.attributes.vals[attrkey],
                    sep=":", end="; ")
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
        for variantSet in self.getAllVariantSets():
            iterator = self._client.searchVariantAnnotationSets(
                variantSetId=variantSet.id)
            for variantAnnotationSet in iterator:
                yield variantAnnotationSet

    def run(self):
        if self._variantAnnotationSetId is None:
            for annotationSet in self.getAllAnnotationSets():
                self._run(annotationSet.id)
        else:
            self._run(self._variantAnnotationSetId)


class SearchFeaturesRunner(FeatureFormatterMixin, AbstractSearchRunner):
    """
    Runner class for the features/search method.
    """
    def __init__(self, args):
        super(SearchFeaturesRunner, self).__init__(args)
        self._referenceName = args.referenceName
        self._featureSetId = args.featureSetId
        self._parentId = args.parentId
        self._start = args.start
        self._end = args.end
        if args.featureTypes == "":
            self._featureTypes = []
        else:
            self._featureTypes = args.featureTypes.split(",")

    def _run(self, featureSetId):
        iterator = self._client.searchFeatures(
            start=self._start, end=self._end,
            referenceName=self._referenceName,
            featureSetId=featureSetId, parentId=self._parentId,
            featureTypes=self._featureTypes)
        self._output(iterator)

    def run(self):
        if self._featureSetId is None and self._parentId is None:
            for featureSet in self.getAllFeatureSets():
                self._run(featureSet)
        else:
            self._run(self._featureSetId)


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
            iterator = self._client.searchReferences(rg.reference_set_id)
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


class GetVariantAnnotationSetRunner(AbstractGetRunner):
    """
    Runner class for the variantannotationsets/{id} method
    """
    def __init__(self, args):
        super(GetVariantAnnotationSetRunner, self).__init__(args)
        self._method = self._client.getVariantAnnotationSet


class GetFeatureRunner(FeatureFormatterMixin, AbstractGetRunner):
    """
    Runner class for the features/{id} method
    """
    def __init__(self, args):
        super(GetFeatureRunner, self).__init__(args)
        self._method = self._client.getFeature


class GetFeatureSetRunner(AbstractGetRunner):
    """
    Runner class for the featuresets/{id} method
    """
    def __init__(self, args):
        super(GetFeatureSetRunner, self).__init__(args)
        self._method = self._client.getFeatureSet


def addDisableUrllibWarningsArgument(parser):
    parser.add_argument(
        "--disable-urllib-warnings", default=False, action="store_true",
        help="Disable urllib3 warnings")


def addVersionArgument(parser):
    # TODO argparse strips newlines from version output
    versionString = (
        "GA4GH Server Version {}\n"
        "(Protocol Version {})".format(
            ga4gh.__version__, protocol.version))
    parser.add_argument(
        "--version", version=versionString, action="version")


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
    addEffectsArgument(parser)
    addPageSizeArgument(parser)


def addFeaturesSearchOptions(parser):
    """
    Adds common options to a features search command line parser.
    """
    addFeatureSetIdArgument(parser)
    addReferenceNameArgument(parser)
    addStartArgument(parser)
    addEndArgument(parser)
    addParentFeatureIdArgument(parser)
    addFeatureTypesArgument(parser)


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
        help="The variant annotation set id to search over")


def addFeatureSetIdArgument(parser):
    parser.add_argument(
        "--featureSetId", "-F", default=None,
        help="The feature set id to search over")


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


def addFeatureTypesArgument(parser):
    parser.add_argument(
        "--featureTypes", "-t", default="",
        help="""Return features matching any of the supplied
            feature types (ontology terms).
            Pass in terms as a comma separated list (no spaces).
            """)


def addParentFeatureIdArgument(parser):
    parser.add_argument(
        "--parentId", "-p", default=None,
        help="Filter features by supplied parent ID")


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
    addDisableUrllibWarningsArgument(parser)
    addVersionArgument(parser)


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
        description="Search for variant annotations",
        help="Search for variant annotations.")
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


def addVariantAnnotationSetsGetParser(subparsers):
    parser = addSubparser(
        subparsers, "variantannotationsets-get", "Get a variantAnnotationSet")
    parser.set_defaults(runner=GetVariantAnnotationSetRunner)
    addGetArguments(parser)


def addVariantSetsGetParser(subparsers):
    parser = addSubparser(
        subparsers, "variantsets-get", "Get a variantSet")
    parser.set_defaults(runner=GetVariantSetRunner)
    addGetArguments(parser)


def addFeaturesGetParser(subparsers):
    parser = addSubparser(
        subparsers, "features-get", "Get a feature by ID")
    parser.set_defaults(runner=GetFeatureRunner)
    addGetArguments(parser)


def addFeatureSetsGetParser(subparsers):
    parser = addSubparser(
        subparsers, "featuresets-get", "Get a featureSet by ID")
    parser.set_defaults(runner=GetFeatureSetRunner)
    addGetArguments(parser)


def addFeaturesSearchParser(subparsers):
    parser = subparsers.add_parser(
        "features-search",
        description="Search for sequence annotation features",
        help="Search for sequence annotation features.")
    parser.set_defaults(runner=SearchFeaturesRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addPageSizeArgument(parser)
    addFeaturesSearchOptions(parser)
    return parser


def addFeatureSetsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "featuresets-search",
        description="Search for sequence annotation feature sets",
        help="Search for featureSets.")
    parser.set_defaults(runner=SearchFeatureSetsRunner)
    addOutputFormatArgument(parser)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addDatasetIdArgument(parser)
    return parser


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
    addVariantSetIdArgument(parser)
    return parser


def addReadsSearchParser(subparsers):
    parser = addSubparser(
        subparsers, "reads-search", "Search for reads")
    parser.set_defaults(runner=SearchReadsRunner)
    addOutputFormatArgument(parser)
    addReadsSearchParserArguments(parser)
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
    parser = createArgumentParser("GA4GH reference client")
    addClientGlobalOptions(parser)
    subparsers = parser.add_subparsers(title='subcommands',)
    addHelpParser(subparsers)
    addVariantsSearchParser(subparsers)
    addVariantSetsSearchParser(subparsers)
    addVariantAnnotationSearchParser(subparsers)
    addVariantAnnotationSetsSearchParser(subparsers)
    addVariantSetsGetParser(subparsers)
    addVariantAnnotationSetsGetParser(subparsers)
    addFeaturesSearchParser(subparsers)
    addFeaturesGetParser(subparsers)
    addFeatureSetsGetParser(subparsers)
    addFeatureSetsSearchParser(subparsers)
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
    parser = createArgumentParser((
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
    parser = createArgumentParser("GA4GH SAM conversion tool")
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
        parser = createArgumentParser(
            "GA4GH server configuration validator")
    parser.add_argument(
        "--config", "-c", default='DevelopmentConfig', type=str,
        help="The configuration to use")
    parser.add_argument(
        "--config-file", "-f", type=str, default=None,
        help="The configuration file to use")
    addVersionArgument(parser)

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


def getNameFromPath(filePath):
    """
    Returns the filename of the specified path without its extensions.
    This is usually how we derive the default name for a given object.
    """
    if len(filePath) == 0:
        raise ValueError("Cannot have empty path for name")
    fileName = os.path.split(os.path.normpath(filePath))[1]
    # We need to handle things like .fa.gz, so we can't use
    # os.path.splitext
    ret = fileName.split(".")[0]
    assert ret != ""
    return ret


def getRawInput(display):
    """
    Wrapper around raw_input; put into separate function so that it
    can be easily mocked for tests.
    """
    return raw_input(display)


class RepoManager(object):
    """
    Class that provide command line functionality to manage a
    data repository.
    """
    def __init__(self, args):
        self._args = args
        self._registryPath = args.registryPath
        self._repo = datarepo.SqlDataRepository(self._registryPath)

    def _confirmDelete(self, objectType, name, func):
        if self._args.force:
            func()
        else:
            displayString = (
                "Are you sure you want to delete the {} '{}'? "
                "[y|N] ".format(objectType, name))
            userResponse = getRawInput(displayString)
            if userResponse.strip() == 'y':
                func()
            else:
                print("Aborted")

    def _updateRepo(self, func, *args, **kwargs):
        """
        Runs the specified function that updates the repo with the specified
        arguments. This method ensures that all updates are transactional,
        so that if any part of the update fails no changes are made to the
        repo.
        """
        # TODO how do we make this properly transactional?
        self._repo.open(datarepo.MODE_WRITE)
        try:
            func(*args, **kwargs)
            self._repo.commit()
        finally:
            self._repo.close()

    def _openRepo(self):
        if not self._repo.exists():
            raise exceptions.RepoManagerException(
                "Repo '{}' does not exist. Please create a new repo "
                "using the 'init' command.".format(self._registryPath))
        self._repo.open(datarepo.MODE_READ)

    def _checkSequenceOntology(self, ontology):
        so = ontologies.SEQUENCE_ONTOLOGY_PREFIX
        if ontology.getOntologyPrefix() != so:
            raise exceptions.RepoManagerException(
                "Ontology '{}' does not have ontology prefix '{}'".format(
                    ontology.getName(), so))

    def _getFilePath(self, filePath, useRelativePath):
        return filePath if useRelativePath else os.path.abspath(filePath)

    def init(self):
        forceMessage = (
            "Respository '{}' already exists. Use --force to overwrite")
        if self._repo.exists():
            if self._args.force:
                self._repo.delete()
            else:
                raise exceptions.RepoManagerException(
                    forceMessage.format(self._registryPath))
        self._updateRepo(self._repo.initialise)

    def list(self):
        """
        Lists the contents of this repo.
        """
        self._openRepo()
        # TODO this is _very_ crude. We need much more options and detail here.
        self._repo.printSummary()

    def verify(self):
        """
        Checks that the data pointed to in the repository works and
        we don't have any broken URLs, missing files, etc.
        """
        self._openRepo()
        self._repo.verify()

    def addOntology(self):
        """
        Adds a new Ontology to this repo.
        """
        self._openRepo()
        name = self._args.name
        filePath = self._getFilePath(self._args.filePath,
                                     self._args.relativePath)
        if name is None:
            name = getNameFromPath(filePath)
        ontology = ontologies.Ontology(name)
        ontology.populateFromFile(filePath)
        self._updateRepo(self._repo.insertOntology, ontology)

    def addDataset(self):
        """
        Adds a new dataset into this repo.
        """
        self._openRepo()
        dataset = datasets.Dataset(self._args.datasetName)
        dataset.setDescription(self._args.description)
        self._updateRepo(self._repo.insertDataset, dataset)

    def addReferenceSet(self):
        """
        Adds a new reference set into this repo.
        """
        self._openRepo()
        name = self._args.name
        filePath = self._getFilePath(self._args.filePath,
                                     self._args.relativePath)
        if name is None:
            name = getNameFromPath(self._args.filePath)
        referenceSet = references.HtslibReferenceSet(name)
        referenceSet.populateFromFile(filePath)
        referenceSet.setDescription(self._args.description)
        referenceSet.setNcbiTaxonId(self._args.ncbiTaxonId)
        referenceSet.setIsDerived(self._args.isDerived)
        referenceSet.setAssemblyId(self._args.assemblyId)
        sourceAccessions = []
        if self._args.sourceAccessions is not None:
            sourceAccessions = self._args.sourceAccessions.split(",")
        referenceSet.setSourceAccessions(sourceAccessions)
        referenceSet.setSourceUri(self._args.sourceUri)
        self._updateRepo(self._repo.insertReferenceSet, referenceSet)

    def addReadGroupSet(self):
        """
        Adds a new ReadGroupSet into this repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        dataUrl = self._args.dataFile
        indexFile = self._args.indexFile
        parsed = urlparse.urlparse(dataUrl)
        # TODO, add https support and others when they have been
        # tested.
        if parsed.scheme in ['http', 'ftp']:
            if indexFile is None:
                raise exceptions.MissingIndexException(dataUrl)
        else:
            if indexFile is None:
                indexFile = dataUrl + ".bai"
            dataUrl = self._getFilePath(self._args.dataFile,
                                        self._args.relativePath)
            indexFile = self._getFilePath(indexFile, self._args.relativePath)
        name = self._args.name
        if self._args.name is None:
            name = getNameFromPath(dataUrl)
        readGroupSet = reads.HtslibReadGroupSet(dataset, name)
        readGroupSet.populateFromFile(dataUrl, indexFile)
        referenceSetName = self._args.referenceSetName
        if referenceSetName is None:
            # Try to find a reference set name from the BAM header.
            referenceSetName = readGroupSet.getBamHeaderReferenceSetName()
        referenceSet = self._repo.getReferenceSetByName(referenceSetName)
        readGroupSet.setReferenceSet(referenceSet)
        self._updateRepo(self._repo.insertReadGroupSet, readGroupSet)

    def addVariantSet(self):
        """
        Adds a new VariantSet into this repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        dataUrls = self._args.dataFiles
        name = self._args.name
        if len(dataUrls) == 1:
            if self._args.name is None:
                name = getNameFromPath(dataUrls[0])
            if os.path.isdir(dataUrls[0]):
                # Read in the VCF files from the directory.
                # TODO support uncompressed VCF and BCF files
                vcfDir = dataUrls[0]
                pattern = os.path.join(vcfDir, "*.vcf.gz")
                dataUrls = glob.glob(pattern)
                if len(dataUrls) == 0:
                    raise exceptions.RepoManagerException(
                        "Cannot find any VCF files in the directory "
                        "'{}'.".format(vcfDir))
                dataUrls[0] = self._getFilePath(dataUrls[0],
                                                self._args.relativePath)
        elif self._args.name is None:
            raise exceptions.RepoManagerException(
                "Cannot infer the intended name of the VariantSet when "
                "more than one VCF file is provided. Please provide a "
                "name argument using --name.")
        parsed = urlparse.urlparse(dataUrls[0])
        if parsed.scheme not in ['http', 'ftp']:
            dataUrls = map(lambda url: self._getFilePath(
                url, self._args.relativePath), dataUrls)
        # Now, get the index files for the data files that we've now obtained.
        indexFiles = self._args.indexFiles
        if indexFiles is None:
            # First check if all the paths exist locally, as they must
            # if we are making a default index path.
            for dataUrl in dataUrls:
                if not os.path.exists(dataUrl):
                    raise exceptions.MissingIndexException(
                        "Cannot find file '{}'. All variant files must be "
                        "stored locally if the default index location is "
                        "used. If you are trying to create a VariantSet "
                        "based on remote URLs, please download the index "
                        "files to the local file system and provide them "
                        "with the --indexFiles argument".format(dataUrl))
            # We assume that the indexes are made by adding .tbi
            indexSuffix = ".tbi"
            # TODO support BCF input properly here by adding .csi
            indexFiles = [filename + indexSuffix for filename in dataUrls]
        indexFiles = map(lambda url: self._getFilePath(
            url, self._args.relativePath), indexFiles)
        variantSet = variants.HtslibVariantSet(dataset, name)
        variantSet.populateFromFile(dataUrls, indexFiles)
        # Get the reference set that is associated with the variant set.
        referenceSetName = self._args.referenceSetName
        if referenceSetName is None:
            # Try to find a reference set name from the VCF header.
            referenceSetName = variantSet.getVcfHeaderReferenceSetName()
        if referenceSetName is None:
            raise exceptions.RepoManagerException(
                "Cannot infer the ReferenceSet from the VCF header. Please "
                "specify the ReferenceSet to associate with this "
                "VariantSet using the --referenceSetName option")
        referenceSet = self._repo.getReferenceSetByName(referenceSetName)
        variantSet.setReferenceSet(referenceSet)

        # Now check for annotations
        annotationSets = []
        if variantSet.isAnnotated() and self._args.addAnnotationSets:
            ontologyName = self._args.ontologyName
            if ontologyName is None:
                raise exceptions.RepoManagerException(
                    "A sequence ontology name must be provided")
            ontology = self._repo.getOntologyByName(ontologyName)
            self._checkSequenceOntology(ontology)
            for annotationSet in variantSet.getVariantAnnotationSets():
                annotationSet.setOntology(ontology)
                annotationSets.append(annotationSet)

        # Add the annotation sets and the variant set as an atomic update
        def updateRepo():
            self._repo.insertVariantSet(variantSet)
            for annotationSet in annotationSets:
                self._repo.insertVariantAnnotationSet(annotationSet)
        self._updateRepo(updateRepo)

    def removeReferenceSet(self):
        """
        Removes a referenceSet from the repo.
        """
        self._openRepo()
        referenceSet = self._repo.getReferenceSetByName(
            self._args.referenceSetName)

        def func():
            self._updateRepo(self._repo.removeReferenceSet, referenceSet)
        self._confirmDelete("ReferenceSet", referenceSet.getLocalId(), func)

    def removeReadGroupSet(self):
        """
        Removes a readGroupSet from the repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        readGroupSet = dataset.getReadGroupSetByName(
            self._args.readGroupSetName)

        def func():
            self._updateRepo(self._repo.removeReadGroupSet, readGroupSet)
        self._confirmDelete("ReadGroupSet", readGroupSet.getLocalId(), func)

    def removeVariantSet(self):
        """
        Removes a variantSet from the repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        variantSet = dataset.getVariantSetByName(self._args.variantSetName)

        def func():
            self._updateRepo(self._repo.removeVariantSet, variantSet)
        self._confirmDelete("VariantSet", variantSet.getLocalId(), func)

    def removeDataset(self):
        """
        Removes a dataset from the repo.
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)

        def func():
            self._updateRepo(self._repo.removeDataset, dataset)
        self._confirmDelete("Dataset", dataset.getLocalId(), func)

    def addFeatureSet(self):
        """
        Adds a new feature set into this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        filePath = self._getFilePath(self._args.filePath,
                                     self._args.relativePath)
        name = getNameFromPath(self._args.filePath)
        featureSet = sequenceAnnotations.Gff3DbFeatureSet(
            dataset, name)
        referenceSetName = self._args.referenceSetName
        if referenceSetName is None:
            raise exceptions.RepoManagerException(
                "A reference set name must be provided")
        referenceSet = self._repo.getReferenceSetByName(referenceSetName)
        featureSet.setReferenceSet(referenceSet)
        ontologyName = self._args.ontologyName
        if ontologyName is None:
            raise exceptions.RepoManagerException(
                "A sequence ontology name must be provided")
        ontology = self._repo.getOntologyByName(ontologyName)
        self._checkSequenceOntology(ontology)
        featureSet.setOntology(ontology)
        featureSet.populateFromFile(filePath)
        self._updateRepo(self._repo.insertFeatureSet, featureSet)

    def removeFeatureSet(self):
        """
        Removes a feature set from this repo
        """
        self._openRepo()
        dataset = self._repo.getDatasetByName(self._args.datasetName)
        featureSet = dataset.getFeatureSetByName(self._args.featureSetName)

        def func():
            self._updateRepo(self._repo.removeFeatureSet, featureSet)
        self._confirmDelete("FeatureSet", featureSet.getLocalId(), func)

    def removeOntology(self):
        """
        Removes an ontology from the repo.
        """
        self._openRepo()
        ontology = self._repo.getOntologyByName(self._args.ontologyName)

        def func():
            self._updateRepo(self._repo.removeOntology, ontology)
        self._confirmDelete("Ontology", ontology.getName(), func)

    #
    # Methods to simplify adding common arguments to the parser.
    #

    @classmethod
    def addRepoArgument(cls, subparser):
        subparser.add_argument(
            "registryPath",
            help="the location of the registry database")

    @classmethod
    def addForceOption(cls, subparser):
        subparser.add_argument(
            "-f", "--force", action='store_true',
            default=False, help="do not prompt for confirmation")

    @classmethod
    def addRelativePathOption(cls, subparser):
        subparser.add_argument(
            "-r", "--relativePath", action='store_true',
            default=False, help="store relative path in database")

    @classmethod
    def addDescriptionOption(cls, subparser, objectType):
        subparser.add_argument(
            "-d", "--description", default="",
            help="The human-readable description of the {}.".format(
                objectType))

    @classmethod
    def addDatasetNameArgument(cls, subparser):
        subparser.add_argument(
            "datasetName", help="the name of the dataset")

    @classmethod
    def addReferenceSetNameOption(cls, subparser, objectType):
        helpText = (
            "the name of the reference set to associate with this {}"
        ).format(objectType)
        subparser.add_argument(
            "-R", "--referenceSetName", default=None, help=helpText)

    @classmethod
    def addSequenceOntologyNameOption(cls, subparser, objectType):
        helpText = (
            "the name of the sequence ontology instance used to "
            "translate ontology term names to IDs in this {}"
        ).format(objectType)
        subparser.add_argument(
            "-O", "--ontologyName", default=None, help=helpText)

    @classmethod
    def addOntologyNameArgument(cls, subparser):
        subparser.add_argument(
            "ontologyName",
            help="the name of the ontology")

    @classmethod
    def addReadGroupSetNameArgument(cls, subparser):
        subparser.add_argument(
            "readGroupSetName",
            help="the name of the read group set")

    @classmethod
    def addVariantSetNameArgument(cls, subparser):
        subparser.add_argument(
            "variantSetName",
            help="the name of the feature set")

    @classmethod
    def addFeatureSetNameArgument(cls, subparser):
        subparser.add_argument(
            "featureSetName",
            help="the name of the variant set")

    @classmethod
    def addFilePathArgument(cls, subparser, helpText):
        subparser.add_argument("filePath", help=helpText)

    @classmethod
    def addNameOption(cls, parser, objectType):
        parser.add_argument(
            "-n", "--name", default=None,
            help="The name of the {}".format(objectType))

    @classmethod
    def getParser(cls):
        parser = createArgumentParser(
            "GA4GH data repository management tool")
        subparsers = parser.add_subparsers(title='subcommands',)
        addVersionArgument(parser)

        initParser = addSubparser(
            subparsers, "init", "Initialize a data repository")
        initParser.set_defaults(runner="init")
        cls.addRepoArgument(initParser)
        cls.addForceOption(initParser)

        verifyParser = addSubparser(
            subparsers, "verify",
            "Verifies the repository by examing all data files")
        verifyParser.set_defaults(runner="verify")
        cls.addRepoArgument(verifyParser)

        listParser = addSubparser(
            subparsers, "list", "List the contents of the repo")
        listParser.set_defaults(runner="list")
        cls.addRepoArgument(listParser)

        addDatasetParser = addSubparser(
            subparsers, "add-dataset", "Add a dataset to the data repo")
        addDatasetParser.set_defaults(runner="addDataset")
        cls.addRepoArgument(addDatasetParser)
        cls.addDatasetNameArgument(addDatasetParser)
        cls.addDescriptionOption(addDatasetParser, "dataset")

        removeDatasetParser = addSubparser(
            subparsers, "remove-dataset",
            "Remove a dataset from the data repo")
        removeDatasetParser.set_defaults(runner="removeDataset")
        cls.addRepoArgument(removeDatasetParser)
        cls.addDatasetNameArgument(removeDatasetParser)
        cls.addForceOption(removeDatasetParser)

        objectType = "reference set"
        addReferenceSetParser = addSubparser(
            subparsers, "add-referenceset",
            "Add a reference set to the data repo")
        addReferenceSetParser.set_defaults(runner="addReferenceSet")
        cls.addRepoArgument(addReferenceSetParser)
        cls.addFilePathArgument(
            addReferenceSetParser,
            "The path of the FASTA file to use as a reference set. This "
            "file must be bgzipped and indexed.")
        cls.addRelativePathOption(addReferenceSetParser)
        cls.addNameOption(addReferenceSetParser, objectType)
        cls.addDescriptionOption(addReferenceSetParser, objectType)
        addReferenceSetParser.add_argument(
            "--ncbiTaxonId", default=None, help="The NCBI Taxon Id")
        addReferenceSetParser.add_argument(
            "--isDerived", default=False, type=bool,
            help="Indicates if this reference set is derived from another")
        addReferenceSetParser.add_argument(
            "--assemblyId", default=None,
            help="The assembly id")
        addReferenceSetParser.add_argument(
            "--sourceAccessions", default=None,
            help="The source accessions (pass as comma-separated list)")
        addReferenceSetParser.add_argument(
            "--sourceUri", default=None,
            help="The source URI")

        removeReferenceSetParser = addSubparser(
            subparsers, "remove-referenceset",
            "Remove a reference set from the repo")
        removeReferenceSetParser.set_defaults(runner="removeReferenceSet")
        cls.addRepoArgument(removeReferenceSetParser)
        removeReferenceSetParser.add_argument(
            "referenceSetName",
            help="the name of the reference set")
        cls.addForceOption(removeReferenceSetParser)

        objectType = "ReadGroupSet"
        addReadGroupSetParser = addSubparser(
            subparsers, "add-readgroupset",
            "Add a read group set to the data repo")
        addReadGroupSetParser.set_defaults(runner="addReadGroupSet")
        cls.addRepoArgument(addReadGroupSetParser)
        cls.addDatasetNameArgument(addReadGroupSetParser)
        cls.addNameOption(addReadGroupSetParser, objectType)
        cls.addReferenceSetNameOption(addReadGroupSetParser, "ReadGroupSet")
        cls.addRelativePathOption(addReadGroupSetParser)
        addReadGroupSetParser.add_argument(
            "dataFile",
            help="The file path or URL of the BAM file for this ReadGroupSet")
        addReadGroupSetParser.add_argument(
            "-I", "--indexFile", default=None,
            help=(
                "The file path of the BAM index for this ReadGroupSet. "
                "If the dataFile argument is a local file, this will "
                "be automatically inferred by appending '.bai' to the "
                "file name. If the dataFile is a remote URL the path to "
                "a local file containing the BAM index must be provided"))

        addOntologyParser = addSubparser(
            subparsers, "add-ontology",
            "Adds an ontology in OBO format to the repo. Currently, "
            "a sequence ontology (SO) instance is required to translate "
            "ontology term names held in annotations to ontology IDs. "
            "Sequence ontology files can be found at "
            "https://github.com/The-Sequence-Ontology/SO-Ontologies")
        addOntologyParser.set_defaults(runner="addOntology")
        cls.addRepoArgument(addOntologyParser)
        cls.addFilePathArgument(
            addOntologyParser,
            "The path of the OBO file defining this ontology.")
        cls.addRelativePathOption(addOntologyParser)
        cls.addNameOption(addOntologyParser, "ontology")

        removeOntologyParser = addSubparser(
            subparsers, "remove-ontology",
            "Remove an ontology from the repo")
        removeOntologyParser.set_defaults(runner="removeOntology")
        cls.addRepoArgument(removeOntologyParser)
        cls.addOntologyNameArgument(removeOntologyParser)
        cls.addForceOption(removeOntologyParser)

        removeReadGroupSetParser = addSubparser(
            subparsers, "remove-readgroupset",
            "Remove a read group set from the repo")
        removeReadGroupSetParser.set_defaults(runner="removeReadGroupSet")
        cls.addRepoArgument(removeReadGroupSetParser)
        cls.addDatasetNameArgument(removeReadGroupSetParser)
        cls.addReadGroupSetNameArgument(removeReadGroupSetParser)
        cls.addForceOption(removeReadGroupSetParser)

        objectType = "VariantSet"
        addVariantSetParser = addSubparser(
            subparsers, "add-variantset",
            "Add a variant set to the data repo based on one or "
            "more VCF files. ")
        addVariantSetParser.set_defaults(runner="addVariantSet")
        cls.addRepoArgument(addVariantSetParser)
        cls.addDatasetNameArgument(addVariantSetParser)
        cls.addRelativePathOption(addVariantSetParser)
        addVariantSetParser.add_argument(
            "dataFiles", nargs="+",
            help=(
                "The VCF/BCF files representing the new VariantSet. "
                "These may be specified either one or more paths "
                "to local files or remote URLS, or as a path to "
                "a local directory containing VCF files. Either "
                "a single directory argument may be passed or a "
                "list of file paths/URLS, but not a mixture of "
                "directories and paths.")
            )
        addVariantSetParser.add_argument(
            "-I", "--indexFiles", nargs="+", metavar="indexFiles",
            help=(
                "The index files for the VCF/BCF files provided in "
                "the dataFiles argument. These must be provided in the "
                "same order as the data files."))
        cls.addNameOption(addVariantSetParser, objectType)
        cls.addReferenceSetNameOption(addVariantSetParser, objectType)
        cls.addSequenceOntologyNameOption(addVariantSetParser, objectType)
        addVariantSetParser.add_argument(
            "-a", "--addAnnotationSets", action="store_true",
            help=(
                "If the supplied VCF file contains annotations, create the "
                "corresponding VariantAnnotationSet."))

        removeVariantSetParser = addSubparser(
            subparsers, "remove-variantset",
            "Remove a variant set from the repo")
        removeVariantSetParser.set_defaults(runner="removeVariantSet")
        cls.addRepoArgument(removeVariantSetParser)
        cls.addDatasetNameArgument(removeVariantSetParser)
        cls.addVariantSetNameArgument(removeVariantSetParser)
        cls.addForceOption(removeVariantSetParser)

        addFeatureSetParser = addSubparser(
            subparsers, "add-featureset", "Add a feature set to the data repo")
        addFeatureSetParser.set_defaults(runner="addFeatureSet")
        cls.addRepoArgument(addFeatureSetParser)
        cls.addDatasetNameArgument(addFeatureSetParser)
        cls.addRelativePathOption(addFeatureSetParser)
        cls.addFilePathArgument(
            addFeatureSetParser,
            "The path to the converted SQLite database containing Feature "
            "data")
        cls.addReferenceSetNameOption(addFeatureSetParser, "feature set")
        cls.addSequenceOntologyNameOption(addFeatureSetParser, "feature set")

        removeFeatureSetParser = addSubparser(
            subparsers, "remove-featureset",
            "Remove a feature set from the repo")
        removeFeatureSetParser.set_defaults(runner="removeFeatureSet")
        cls.addRepoArgument(removeFeatureSetParser)
        cls.addDatasetNameArgument(removeFeatureSetParser)
        cls.addFeatureSetNameArgument(removeFeatureSetParser)
        cls.addForceOption(removeFeatureSetParser)

        return parser

    @classmethod
    def runCommand(cls, args):
        parser = cls.getParser()
        parsedArgs = parser.parse_args(args)
        if "runner" not in parsedArgs:
            parser.print_help()
        manager = RepoManager(parsedArgs)
        runMethod = getattr(manager, parsedArgs.runner)
        runMethod()


def getRepoManagerParser():
    """
    Used by sphinx.argparse.
    """
    return RepoManager.getParser()


def repoExitError(message):
    """
    Exits the repo manager with error status.
    """
    wrapper = textwrap.TextWrapper(
        break_on_hyphens=False, break_long_words=False)
    formatted = wrapper.fill("{}: error: {}".format(sys.argv[0], message))
    sys.exit(formatted)


def repo_main(args=None):
    try:
        RepoManager.runCommand(args)
    except exceptions.RepoManagerException as exception:
        # These are exceptions that happen throughout the CLI, and are
        # used to communicate back to the user
        repoExitError(str(exception))
    except exceptions.NotFoundException as exception:
        # We expect NotFoundExceptions to occur when we're looking for
        # datasets, readGroupsets, etc.
        repoExitError(str(exception))
    except exceptions.DataException as exception:
        # We expect DataExceptions to occur when a file open fails,
        # a URL cannot be reached, etc.
        repoExitError(str(exception))
    except Exception as exception:
        # Uncaught exception: this is a bug
        message = """
An internal error has occurred.  Please file a bug report at
https://github.com/ga4gh/server/issues
with all the relevant details, and the following stack trace.
"""
        print("{}: error:".format(sys.argv[0]), file=sys.stderr)
        print(message, file=sys.stderr)
        traceback.print_exception(*sys.exc_info())
        sys.exit(1)
