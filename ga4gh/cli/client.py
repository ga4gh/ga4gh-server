"""
Client cli
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import logging
import requests

import ga4gh.cli as cli
import ga4gh.backend as backend
import ga4gh.client as client
import ga4gh.datarepo as datarepo
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol


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
        self._client.set_page_size(self._pageSize)

    def getAllDatasets(self):
        """
        Returns all datasets on the server.
        """
        return self._client.search_datasets()

    def getAllVariantSets(self):
        """
        Returns all variant sets on the server.
        """
        for dataset in self.getAllDatasets():
            iterator = self._client.search_variant_sets(
                dataset_id=dataset.id)
            for variantSet in iterator:
                yield variantSet

    def getAllFeatureSets(self):
        """
        Returns all feature sets on the server.
        """
        for dataset in self.getAllDatasets():
            iterator = self._client.search_feature_sets(
                dataset_id=dataset.id)
            for featureSet in iterator:
                yield featureSet

    def getAllReadGroupSets(self):
        """
        Returns all readgroup sets on the server.
        """
        for dataset in self.getAllDatasets():
            iterator = self._client.search_read_group_sets(
                dataset_id=dataset.id)
            for readGroupSet in iterator:
                yield readGroupSet

    def getAllReadGroups(self):
        """
        Get all read groups in a read group set
        """
        for dataset in self.getAllDatasets():
            iterator = self._client.search_read_group_sets(
                dataset_id=dataset.id)
            for readGroupSet in iterator:
                readGroupSet = self._client.get_read_group_set(
                    readGroupSet.id)
                for readGroup in readGroupSet.read_groups:
                    yield readGroup.id

    def getAllReferenceSets(self):
        """
        Returns all reference sets on the server.
        """
        return self._client.search_reference_sets()

# Runners for the various search methods


class SearchDatasetsRunner(AbstractSearchRunner):
    """
    Runner class for the datasets/search method
    """
    def __init__(self, args):
        super(SearchDatasetsRunner, self).__init__(args)

    def run(self):
        iterator = self._client.search_datasets()
        self._output(iterator)


class SearchReferenceSetsRunner(AbstractSearchRunner):
    """
    Runner class for the referencesets/search method.
    """
    def __init__(self, args):
        super(SearchReferenceSetsRunner, self).__init__(args)
        self._accession = args.accession
        self._md5checksum = args.md5checksum
        self._assemblyId = args.assemblyId

    def run(self):
        iterator = self._client.search_reference_sets(
            accession=self._accession,
            md5checksum=self._md5checksum,
            assembly_id=self._assemblyId)
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
        iterator = self._client.search_references(
            accession=self._accession, md5checksum=self._md5checksum,
            reference_set_id=referenceSetId)
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
        iterator = self._client.search_variant_sets(dataset_id=datasetId)
        self._output(iterator)

    def run(self):
        if self._datasetId is None:
            for dataset in self.getAllDatasets():
                self._run(dataset.id)
        else:
            self._run(self._datasetId)


class SearchBioSamplesRunner(AbstractSearchRunner):
    """
    Runner class for the biosamples/search method.
    """
    def __init__(self, args):
        super(SearchBioSamplesRunner, self).__init__(args)
        self._datasetId = args.datasetId
        self._individualId = args.individualId
        self._name = args.name

    def _run(self, datasetId):
        iterator = self._client.search_bio_samples(
            datasetId,
            name=self._name,
            individual_id=self._individualId)
        self._output(iterator)

    def run(self):
        if self._datasetId is None:
            for dataset in self.getAllDatasets():
                self._run(dataset.id)
        else:
            self._run(self._datasetId)


class SearchIndividualsRunner(AbstractSearchRunner):
    """
    Runner class for the individuals/search method.
    """
    def __init__(self, args):
        super(SearchIndividualsRunner, self).__init__(args)
        self._datasetId = args.datasetId
        self._name = args.name

    def _run(self, datasetId):
        iterator = self._client.search_bio_samples(
            datasetId,
            name=self._name)
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
        iterator = self._client.search_variant_annotation_sets(
            variant_set_id=variantSetId)
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
        iterator = self._client.search_feature_sets(dataset_id=datasetId)
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
        iterator = self._client.search_read_group_sets(
            dataset_id=datasetId, name=self._name)
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
        iterator = self._client.search_call_sets(
            variant_set_id=variantSetId, name=self._name)
        self._output(iterator)

    def run(self):
        if self._variantSetId is None:
            for variantSet in self.getAllVariantSets():
                self._run(variantSet.id)
        else:
            self._run(self._variantSetId)


class SearchGenotypePhenotypeRunner(AbstractSearchRunner):
    """
    Runner class for the featurephenotypeassociations/search method.
    """
    def __init__(self, args):
        super(SearchGenotypePhenotypeRunner, self).__init__(args)

        # if arg is JSON; parse; else return as string
        def checkJson(value):
            if value is not None:
                try:
                    return json.loads(value)
                except ValueError:
                    return value

        self._feature_ids = None
        self._evidence = None
        self._phenotype_ids = None
        self._phenotype_association_set_id = args.phenotype_association_set_id
        if args.feature_ids:
            self._feature_ids = args.feature_ids.split(",")
        if args.phenotype_ids:
            self._phenotype_ids = args.phenotype_ids.split(",")
        if args.evidence:
            self._evidence = checkJson(args.evidence)

    def run(self):
        iterator = self._client.search_genotype_phenotype(
            phenotype_association_set_id=self._phenotype_association_set_id,
            feature_ids=self._feature_ids,
            phenotype_ids=self._phenotype_ids,
            evidence=self._evidence)
        self._output(iterator)

    def _textOutput(self, gaObjects):
        """
        Prints out the specified FeaturePhenotypeAssociation objects.
        """
        for association in gaObjects:
            print(association.id)


class SearchPhenotypeRunner(AbstractSearchRunner):
    """
    Runner class for the phenotype/search method.
    """
    def __init__(self, args):
        super(SearchPhenotypeRunner, self).__init__(args)

        self._phenotype_association_set_id = args.phenotype_association_set_id
        self._phenotype_id = args.phenotype_id
        self._description = args.description
        self._type = args.type
        self._age_of_onset = args.age_of_onset

    def run(self):
        iterator = self._client.search_phenotype(
            phenotype_association_set_id=self._phenotype_association_set_id,
            phenotype_id=self._phenotype_id,
            description=self._description,
            type_=self._type,
            age_of_onset=self._age_of_onset)
        self._output(iterator)


class SearchPhenotypeAssociationSetsRunner(AbstractSearchRunner):
    """
    Runner class for the phenotypeassociationsets/search method.
    """
    def __init__(self, args):
        super(SearchPhenotypeAssociationSetsRunner, self).__init__(args)
        self._dataset_id = args.datasetId

    def run(self):
        iterator = self._client.search_phenotype_association_sets(
            dataset_id=self._dataset_id)
        self._output(iterator)


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
                val = value.values[0].string_value
                print(key, val, sep="=", end=";")
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
                variantAnnotation.created, sep="\t", end="\t")
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
                feature.id, feature.parent_id, feature.feature_set_id,
                feature.reference_name, feature.start, feature.end,
                feature.strand, "FeatureType:", feature.feature_type.id,
                feature.feature_type.term, sep="\t")


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
        iterator = self._client.search_variants(
            start=self._start, end=self._end,
            reference_name=self._referenceName,
            variant_set_id=variantSetId,
            call_set_ids=self._callSetIds)
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
        iterator = self._client.search_variant_annotations(
            variant_annotation_set_id=variantAnnotationSetId,
            reference_name=self._referenceName,
            reference_id=self._referenceId,
            start=self._start, end=self._end,
            effects=self._effects)
        self._output(iterator)

    def getAllAnnotationSets(self):
        """
        Returns all variant annotation sets on the server.
        """
        for variantSet in self.getAllVariantSets():
            iterator = self._client.search_variant_annotation_sets(
                variant_set_id=variantSet.id)
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
        iterator = self._client.search_features(
            start=self._start, end=self._end,
            reference_name=self._referenceName,
            feature_set_id=featureSetId, parent_id=self._parentId,
            feature_types=self._featureTypes)
        self._output(iterator)

    def run(self):
        if self._featureSetId is None and not self._parentId:
            for featureSet in self.getAllFeatureSets():
                self._run(featureSet.id)
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
            rg = self._client.get_read_group(
                read_group_id=referenceGroupId)
            iterator = self._client.search_references(rg.reference_set_id)
            for reference in iterator:
                self._run(referenceGroupId, reference.id)
        else:
            iterator = self._client.search_reads(
                read_group_ids=[referenceGroupId],
                reference_id=referenceId,
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


class SearchRnaQuantificationSetsRunner(AbstractSearchRunner):
    """
    Runner class for the rnaquantificationsets/search method
    """
    def __init__(self, args):
        super(SearchRnaQuantificationSetsRunner, self).__init__(args)
        self._datasetId = args.datasetId

    def run(self):
        iterator = self._client.search_rna_quantification_sets(
            self._datasetId)
        self._output(iterator)

    def _textOutput(self, rnaQuants):
        for rnaQuant in rnaQuants:
            print(
                rnaQuant.id, rnaQuant.dataset_id, rnaQuant.name,
                sep="\t", end="\t")
            print()


class SearchRnaQuantificationsRunner(AbstractSearchRunner):
    """
    Runner class for the rnaquantifications/search method
    """
    def __init__(self, args):
        super(SearchRnaQuantificationsRunner, self).__init__(args)
        self._rnaQuantificationSetId = args.rnaQuantificationSetId

    def run(self):
        iterator = self._client.search_rna_quantifications(
            self._rnaQuantificationSetId)
        self._output(iterator)

    def _textOutput(self, rnaQuants):
        for rnaQuant in rnaQuants:
            print(
                rnaQuant.id, rnaQuant.description, rnaQuant.name,
                sep="\t", end="\t")
            for featureSet in rnaQuant.featureSetIds:
                print(featureSet, sep=",", end="\t")
            for readGroup in rnaQuant.readGroupIds:
                print(readGroup, sep=",", end="")
            print()


class SearchExpressionLevelsRunner(AbstractSearchRunner):
    """
    Runner class for the expressionlevels/search method
    """
    def __init__(self, args):
        super(SearchExpressionLevelsRunner, self).__init__(args)
        self._rnaQuantificationId = args.rnaQuantificationId
        self._feature_ids = []
        if len(args.featureIds) > 0:
            self._feature_ids = args.featureIds.split(",")
        self.threshold = args.threshold

    def run(self):
        iterator = self._client.search_expression_levels(
            rna_quantification_id=self._rnaQuantificationId,
            feature_ids=self._feature_ids,
            threshold=self.threshold)
        self._output(iterator)

    def _textOutput(self, expressionObjs):
        for expression in expressionObjs:
            print(
                expression.id, expression.expression, expression.name,
                expression.isNormalized, expression.rawReadCount,
                expression.score, expression.units, sep="\t", end="\t")
            print()


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
        sequence = self._client.list_reference_bases(
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
        self._method = self._client.get_reference_set


class GetReferenceRunner(AbstractGetRunner):
    """
    Runner class for the references/{id} method
    """
    def __init__(self, args):
        super(GetReferenceRunner, self).__init__(args)
        self._method = self._client.get_reference


class GetReadGroupSetRunner(AbstractGetRunner):
    """
    Runner class for the readgroupsets/{id} method
    """
    def __init__(self, args):
        super(GetReadGroupSetRunner, self).__init__(args)
        self._method = self._client.get_read_group_set


class GetReadGroupRunner(AbstractGetRunner):
    """
    Runner class for the references/{id} method
    """
    def __init__(self, args):
        super(GetReadGroupRunner, self).__init__(args)
        self._method = self._client.get_read_group


class GetBioSampleRunner(AbstractGetRunner):
    """
    Runner class for the references/{id} method
    """
    def __init__(self, args):
        super(GetBioSampleRunner, self).__init__(args)
        self._method = self._client.getBioSample


class GetIndividualRunner(AbstractGetRunner):
    """
    Runner class for the references/{id} method
    """
    def __init__(self, args):
        super(GetIndividualRunner, self).__init__(args)
        self._method = self._client.get_individual


class GetCallSetRunner(AbstractGetRunner):
    """
    Runner class for the callsets/{id} method
    """
    def __init__(self, args):
        super(GetCallSetRunner, self).__init__(args)
        self._method = self._client.get_call_set


class GetDatasetRunner(AbstractGetRunner):
    """
    Runner class for the datasets/{id} method
    """
    def __init__(self, args):
        super(GetDatasetRunner, self).__init__(args)
        self._method = self._client.get_dataset


class GetVariantRunner(VariantFormatterMixin, AbstractGetRunner):
    """
    Runner class for the variants/{id} method
    """
    def __init__(self, args):
        super(GetVariantRunner, self).__init__(args)
        self._method = self._client.get_variant


class GetVariantSetRunner(AbstractGetRunner):
    """
    Runner class for the variantsets/{id} method
    """
    def __init__(self, args):
        super(GetVariantSetRunner, self).__init__(args)
        self._method = self._client.get_variant_set


class GetVariantAnnotationSetRunner(AbstractGetRunner):
    """
    Runner class for the variantannotationsets/{id} method
    """
    def __init__(self, args):
        super(GetVariantAnnotationSetRunner, self).__init__(args)
        self._method = self._client.get_variant_annotation_set


class GetFeatureRunner(FeatureFormatterMixin, AbstractGetRunner):
    """
    Runner class for the features/{id} method
    """
    def __init__(self, args):
        super(GetFeatureRunner, self).__init__(args)
        self._method = self._client.get_feature


class GetFeatureSetRunner(AbstractGetRunner):
    """
    Runner class for the featuresets/{id} method
    """
    def __init__(self, args):
        super(GetFeatureSetRunner, self).__init__(args)
        self._method = self._client.get_feature_set


class GetRnaQuantificationRunner(AbstractGetRunner):
    """
    Runner class for the rnaquantifications/{id} method
    """
    def __init__(self, args):
        super(GetRnaQuantificationRunner, self).__init__(args)
        self._method = self._client.get_rna_quantification


class GetExpressionLevelRunner(AbstractGetRunner):
    """
    Runner class for the expressionlevels/{id} method
    """
    def __init__(self, args):
        super(GetExpressionLevelRunner, self).__init__(args)
        self._method = self._client.get_expression_level


class GetRnaQuantificationSetRunner(AbstractGetRunner):
    """
    Runner class for the rnaquantificationsets/{id} method
    """
    def __init__(self, args):
        super(GetRnaQuantificationSetRunner, self).__init__(args)
        self._method = self._client.get_rna_quantification_set


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
    addFeaturesReferenceNameArgument(parser)
    addStartArgument(parser)
    addEndArgument(parser)
    addParentFeatureIdArgument(parser)
    addFeatureTypesArgument(parser)


def addGenotypePhenotypeSearchOptions(parser):
    """
    Adds options to a g2p searches command line parser.
    """
    parser.add_argument(
        "--phenotype_association_set_id", "-s", default=None,
        help="Only return associations from this phenotype_association_set.")
    parser.add_argument(
        "--feature_ids", "-f", default=None,
        help="Only return associations for these features.")
    parser.add_argument(
        "--phenotype_ids", "-p", default=None,
        help="Only return associations for these phenotypes.")
    parser.add_argument(
        "--evidence", "-E", default=None,
        help="Only return associations to this evidence.")


def addPhenotypeSearchOptions(parser):
    """
    Adds options to a phenotype searches command line parser.
    """
    parser.add_argument(
        "--phenotype_association_set_id", "-s", default=None,
        help="Only return phenotypes from this phenotype_association_set.")
    parser.add_argument(
        "--phenotype_id", "-p", default=None,
        help="Only return this phenotype.")
    parser.add_argument(
        "--description", "-d", default=None,
        help="Only return phenotypes matching this description.")
    parser.add_argument(
        "--age_of_onset", "-a", default=None,
        help="Only return phenotypes with this age_of_onset.")
    parser.add_argument(
        "--type", "-T", default=None,
        help="Only return phenotypes with this type.")


def addPhenotypeAssociationSetsSearchOptions(parser):
    """
    Adds options to a phenotype assoc. sets searches command line parser.
    """
    addDatasetIdArgument(parser)


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


def addFeaturesReferenceNameArgument(parser):
    parser.add_argument(
        "--referenceName", "-r", default="",
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
        "--parentId", "-p", default="",
        help="Filter features by supplied parent ID")


def addStartArgument(parser):
    parser.add_argument(
        "--start", "-s", default=0, type=int,
        help="The start of the search range (inclusive).")


def addEndArgument(parser, defaultValue=cli.AVRO_LONG_MAX):
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


def addIndividualIdArgument(parser):
    parser.add_argument(
        "--individualId", default=None,
        help="The ID of the individual")


def addBioSampleIdArgument(parser):
    parser.add_argument(
        "--bioSampleId", default=None,
        help="The ID of the biosample")


def addClientGlobalOptions(parser):
    parser.add_argument(
        '--verbose', '-v', action='count', default=0,
        help="Increase verbosity; can be supplied multiple times")
    parser.add_argument(
        "--key", "-k", default='invalid',
        help="Auth Key. Found on server index page.")
    cli.addDisableUrllibWarningsArgument(parser)
    cli.addVersionArgument(parser)


def addHelpParser(subparsers):
    parser = subparsers.add_parser(
        "help", description="ga4gh_client help",
        help="show this help message and exit")
    return parser


def addVariantsSearchParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "variants-search", "Search for variants")
    parser.set_defaults(runner=SearchVariantsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addVariantSearchOptions(parser)
    return parser


def addVariantSetsSearchParser(subparsers):
    parser = cli.addSubparser(
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
    parser = cli.addSubparser(
        subparsers, "variantannotationsets-get", "Get a variantAnnotationSet")
    parser.set_defaults(runner=GetVariantAnnotationSetRunner)
    addGetArguments(parser)


def addVariantSetsGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "variantsets-get", "Get a variantSet")
    parser.set_defaults(runner=GetVariantSetRunner)
    addGetArguments(parser)


def addFeaturesGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "features-get", "Get a feature by ID")
    parser.set_defaults(runner=GetFeatureRunner)
    addGetArguments(parser)


def addFeatureSetsGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "featuresets-get", "Get a featureSet by ID")
    parser.set_defaults(runner=GetFeatureSetRunner)
    addGetArguments(parser)


def addBioSamplesGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "biosamples-get", "Get a biosample by ID")
    parser.set_defaults(runner=GetBioSampleRunner)
    addGetArguments(parser)


def addIndividualsGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "individuals-get", "Get a individual by ID")
    parser.set_defaults(runner=GetIndividualRunner)
    addGetArguments(parser)


def addBioSamplesSearchParser(subparsers):
    parser = subparsers.add_parser(
        "biosamples-search",
        description="Search for biosamples",
        help="Search for biosamples.")
    parser.set_defaults(runner=SearchBioSamplesRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addPageSizeArgument(parser)
    addDatasetIdArgument(parser)
    addNameArgument(parser)
    addIndividualIdArgument(parser)
    return parser


def addIndividualsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "individuals-search",
        description="Search for individuals",
        help="Search for individuals.")
    parser.set_defaults(runner=SearchIndividualsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addDatasetIdArgument(parser)
    addPageSizeArgument(parser)
    addNameArgument(parser)
    return parser


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
    parser = cli.addSubparser(
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
    parser = cli.addSubparser(
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
    parser = cli.addSubparser(
        subparsers, "readgroupsets-search", "Search for readGroupSets")
    parser.set_defaults(runner=SearchReadGroupSetsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addBioSampleIdArgument(parser)
    addPageSizeArgument(parser)
    addDatasetIdArgument(parser)
    addNameArgument(parser)
    return parser


def addCallSetsSearchParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "callsets-search", "Search for callSets")
    parser.set_defaults(runner=SearchCallSetsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addBioSampleIdArgument(parser)
    addPageSizeArgument(parser)
    addNameArgument(parser)
    addVariantSetIdArgument(parser)
    return parser


def addReadsSearchParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "reads-search", "Search for reads")
    parser.set_defaults(runner=SearchReadsRunner)
    addOutputFormatArgument(parser)
    addReadsSearchParserArguments(parser)
    return parser


def addDatasetsGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "datasets-get", "Get a dataset")
    parser.set_defaults(runner=GetDatasetRunner)
    addGetArguments(parser)


def addDatasetsSearchParser(subparsers):
    parser = cli.addSubparser(
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
    parser = cli.addSubparser(
        subparsers, "referencesets-get", "Get a referenceset")
    parser.set_defaults(runner=GetReferenceSetRunner)
    addGetArguments(parser)


def addReferencesGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "references-get", "Get a reference")
    parser.set_defaults(runner=GetReferenceRunner)
    addGetArguments(parser)


def addReadGroupSetsGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "readgroupsets-get", "Get a read group set")
    parser.set_defaults(runner=GetReadGroupSetRunner)
    addGetArguments(parser)


def addReadGroupsGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "readgroups-get", "Get a read group")
    parser.set_defaults(runner=GetReadGroupRunner)
    addGetArguments(parser)


def addCallSetsGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "callsets-get", "Get a callSet")
    parser.set_defaults(runner=GetCallSetRunner)
    addGetArguments(parser)


def addVariantsGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "variants-get", "Get a variant")
    parser.set_defaults(runner=GetVariantRunner)
    addGetArguments(parser)


def addRnaQuantificationSetGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "rnaquantificationsets-get",
        "Get a rna quantification set")
    parser.set_defaults(runner=GetRnaQuantificationSetRunner)
    addGetArguments(parser)


def addRnaQuantificationGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "rnaquantifications-get", "Get a rna quantification")
    parser.set_defaults(runner=GetRnaQuantificationRunner)
    addGetArguments(parser)


def addExpressionLevelGetParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "expressionlevels-get", "Get a expression level")
    parser.set_defaults(runner=GetExpressionLevelRunner)
    addGetArguments(parser)


def addReferencesBasesListParser(subparsers):
    parser = cli.addSubparser(
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


def addRnaQuantificationSetsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "rnaquantificationsets-search",
        description="Search for rna quantification set",
        help="Search for rna quantification set")
    parser.set_defaults(runner=SearchRnaQuantificationSetsRunner)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addDatasetIdArgument(parser)
    addOutputFormatArgument(parser)
    return parser


def addRnaQuantificationsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "rnaquantifications-search",
        description="Search for rna quantification",
        help="Search for rna quantification")
    parser.set_defaults(runner=SearchRnaQuantificationsRunner)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    parser.add_argument(
        "--rnaQuantificationSetId", default=None,
        help="The rnaQuantification set to search over")
    addOutputFormatArgument(parser)
    return parser


def addExpressionLevelsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "expressionlevels-search",
        description="Search for feature expression",
        help="Search for feature expression")
    parser.set_defaults(runner=SearchExpressionLevelsRunner)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addFeatureIdsArgument(parser)
    parser.add_argument(
        "--rnaQuantificationId", default='',
        help="The RNA Quantification Id to search over")
    parser.add_argument(
        "--threshold", default=0.0, type=float,
        help="The minimum value for expression results to report.")
    addOutputFormatArgument(parser)
    return parser


def addGenotypePhenotypeSearchParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "genotypephenotype-search",
        "Search for genotype to phenotype associations")
    parser.set_defaults(runner=SearchGenotypePhenotypeRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addGenotypePhenotypeSearchOptions(parser)
    addPageSizeArgument(parser)
    return parser


def addPhenotypeSearchParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "phenotype-search", "Search for phenotypes")
    parser.set_defaults(runner=SearchPhenotypeRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addPhenotypeSearchOptions(parser)
    addPageSizeArgument(parser)
    return parser


def addPhenotypeAssociationSetsSearchParser(subparsers):
    parser = cli.addSubparser(
        subparsers, "phenotypeassociationsets-search",
        "Search for phenotypeassociationsets")
    parser.set_defaults(runner=SearchPhenotypeAssociationSetsRunner)
    addUrlArgument(parser)
    addOutputFormatArgument(parser)
    addPhenotypeAssociationSetsSearchOptions(parser)
    addPageSizeArgument(parser)


def getClientParser():
    parser = cli.createArgumentParser("GA4GH reference client")
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
    addBioSamplesSearchParser(subparsers)
    addBioSamplesGetParser(subparsers)
    addIndividualsSearchParser(subparsers)
    addIndividualsGetParser(subparsers)
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
    addRnaQuantificationSetGetParser(subparsers)
    addRnaQuantificationGetParser(subparsers)
    addExpressionLevelGetParser(subparsers)
    addReferencesBasesListParser(subparsers)
    addRnaQuantificationSetsSearchParser(subparsers)
    addRnaQuantificationsSearchParser(subparsers)
    addExpressionLevelsSearchParser(subparsers)
    addGenotypePhenotypeSearchParser(subparsers)
    addPhenotypeSearchParser(subparsers)
    addPhenotypeAssociationSetsSearchParser(subparsers)
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
