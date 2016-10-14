"""
Client classes for the GA4GH reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import requests
import posixpath
import logging

import ga4gh.protocol as protocol
import ga4gh.pb as pb
import ga4gh.exceptions as exceptions


class AbstractClient(object):
    """
    The abstract superclass of GA4GH Client objects.
    """

    def __init__(self, log_level=0):
        self._page_size = None
        self._log_level = log_level
        self._protocol_bytes_received = 0
        logging.basicConfig()
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(log_level)

    def _deserialize_response(
            self, json_response_string, protocol_response_class):
        self._protocol_bytes_received += len(json_response_string)
        self._logger.debug("response:{}".format(json_response_string))
        if not json_response_string:
            raise exceptions.EmptyResponseException()
        return protocol.fromJson(json_response_string, protocol_response_class)

    def _run_search_page_request(
            self, protocol_request, object_name, protocol_response_class):
        """
        Runs a complete transaction with the server to obtain a single
        page of search results.
        """
        raise NotImplemented()

    def _run_search_request(
            self, protocol_request, object_name, protocol_response_class):
        """
        Runs the specified request at the specified object_name and
        instantiates an object of the specified class. We yield each object in
        listAttr.  If pages of results are present, repeat this process
        until the pageToken is null.
        """
        not_done = True
        while not_done:
            response_object = self._run_search_page_request(
                protocol_request, object_name, protocol_response_class)
            value_list = getattr(
                response_object,
                protocol.getValueListName(protocol_response_class))
            for extract in value_list:
                yield extract
            not_done = bool(response_object.next_page_token)
            protocol_request.page_token = response_object.next_page_token

    def _run_list_reference_bases_page_request(self, protocol_request):
        """
        Runs a complete transaction with the server to get a single
        page of results for the specified ListReferenceBasesRequest.
        """
        raise NotImplemented()

    def list_reference_bases(self, id_, start=0, end=None):
        """
        Returns an iterator over the bases from the server in the form
        of consecutive strings. This command does not conform to the
        patterns of the other search and get requests, and is implemented
        differently.
        """
        request = protocol.ListReferenceBasesRequest()
        request.start = pb.int(start)
        request.end = pb.int(end)
        request.reference_id = id_
        not_done = True
        # TODO We should probably use a StringIO here to make string buffering
        # a bit more efficient.
        bases_list = []
        while not_done:
            response = self._run_list_reference_bases_page_request(request)
            bases_list.append(response.sequence)
            not_done = bool(response.next_page_token)
            request.page_token = response.next_page_token
        return "".join(bases_list)

    def _run_get_request(self, object_name, protocol_response_class, id_):
        """
        Requests an object from the server and returns the object of
        type protocol_response_class that has id id_.
        Used for requests where a single object is the expected response.
        """
        raise NotImplemented()

    def get_bio_sample(self, bio_sample_id):
        """
        Perform a get request for the given BioSample.

        :param str bio_sample_id: The ID of the BioSample
        :return: The requested BioSample.
        :rtype: :class:`ga4gh.protocol.BioSample`
        """
        return self._run_get_request(
            "biosamples", protocol.BioSample, bio_sample_id)

    def get_individual(self, individual_id):
        """
        Perform a get request for the given Individual.

        :param str individual_id: The ID of the Individual
        :return: The requested Individual.
        :rtype: :class:`ga4gh.protocol.Individual`
        """
        return self._run_get_request(
            "individuals", protocol.Individual, individual_id)

    def get_page_size(self):
        """
        Returns the suggested maximum size of pages of results returned by
        the server.
        """
        return self._page_size

    def set_page_size(self, page_size):
        """
        Sets the requested maximum size of pages of results returned by the
        server to the specified value.
        """
        self._page_size = page_size

    def get_protocol_bytes_received(self):
        """
        Returns the total number of protocol bytes received from the server
        by this client.

        :return: The number of bytes consumed by protocol traffic read from
            the server during the lifetime of this client.
        :rtype: int
        """
        return self._protocol_bytes_received

    def get_dataset(self, dataset_id):
        """
        Returns the Dataset with the specified ID from the server.

        :param str dataset_id: The ID of the Dataset of interest.
        :return: The Dataset of interest.
        :rtype: :class:`ga4gh.protocol.Dataset`
        """
        return self._run_get_request(
            "datasets", protocol.Dataset, dataset_id)

    def get_reference_set(self, reference_set_id):
        """
        Returns the ReferenceSet with the specified ID from the server.

        :param str reference_set_id: The ID of the ReferenceSet of interest.
        :return: The ReferenceSet of interest.
        :rtype: :class:`ga4gh.protocol.ReferenceSet`
        """
        return self._run_get_request(
            "referencesets", protocol.ReferenceSet, reference_set_id)

    def get_reference(self, reference_id):
        """
        Returns the Reference with the specified ID from the server.

        :param str reference_id: The ID of the Reference of interest.
        :return: The Reference of interest.
        :rtype: :class:`ga4gh.protocol.Reference`
        """
        return self._run_get_request(
            "references", protocol.Reference, reference_id)

    def get_read_group_set(self, read_group_set_id):
        """
        Returns the ReadGroupSet with the specified ID from the server.

        :param str read_group_set_id: The ID of the ReadGroupSet of interest.
        :return: The ReadGroupSet of interest.
        :rtype: :class:`ga4gh.protocol.ReadGroupSet`
        """
        return self._run_get_request(
            "readgroupsets", protocol.ReadGroupSet, read_group_set_id)

    def get_read_group(self, read_group_id):
        """
        Returns the ReadGroup with the specified ID from the server.

        :param str read_group_id: The ID of the ReadGroup of interest.
        :return: The ReadGroup of interest.
        :rtype: :class:`ga4gh.protocol.ReadGroup`
        """
        return self._run_get_request(
            "readgroups", protocol.ReadGroup, read_group_id)

    def get_call_set(self, call_set_id):
        """
        Returns the CallSet with the specified ID from the server.

        :param str call_set_id: The ID of the CallSet of interest.
        :return: The CallSet of interest.
        :rtype: :class:`ga4gh.protocol.CallSet`
        """
        return self._run_get_request(
            "callsets", protocol.CallSet, call_set_id)

    def get_variant(self, variant_id):
        """
        Returns the Variant with the specified ID from the server.

        :param str variant_id: The ID of the Variant of interest.
        :return: The Variant of interest.
        :rtype: :class:`ga4gh.protocol.Variant`
        """
        return self._run_get_request(
            "variants", protocol.Variant, variant_id)

    def get_variant_set(self, variant_set_id):
        """
        Returns the VariantSet with the specified ID from the server.

        :param str variant_set_id: The ID of the VariantSet of interest.
        :return: The VariantSet of interest.
        :rtype: :class:`ga4gh.protocol.VariantSet`
        """
        return self._run_get_request(
            "variantsets", protocol.VariantSet, variant_set_id)

    def get_variant_annotation_set(self, variant_annotation_set_id):
        """
        Returns the VariantAnnotationSet with the specified ID from
        the server.

        :param str variant_annotation_set_id: The ID of the
            VariantAnnotationSet of interest.
        :return: The VariantAnnotationSet of interest.
        :rtype: :class:`ga4gh.protocol.VariantAnnotationSet`
        """
        return self._run_get_request(
            "variantannotationsets", protocol.VariantAnnotationSet,
            variant_annotation_set_id)

    def get_feature_set(self, feature_set_id):
        """
        Returns the FeatureSet with the specified ID from the server.

        :param str feature_set_id: The ID of the FeatureSet of interest.
        :return: The FeatureSet of interest.
        :rtype: :class:`ga4gh.protocol.FeatureSet`
        """
        return self._run_get_request(
            "featuresets", protocol.FeatureSet, feature_set_id)

    def get_feature(self, feature_id):
        """
        Returns the feature with the specified ID from the server.

        :param str feature_id: The ID of the requested feature
        :return: The requested ga4gh.protocol.Feature object.
        """
        return self._run_get_request(
            "features", protocol.Feature, feature_id)

    def get_rna_quantification_set(self, rna_quantification_set_id):
        """
        Returns the RnaQuantificationSet with the specified ID from the server.
        :param str rna_quantification_set_id: The ID of the
            RnaQuantificationSet of interest.
        :return: The RnaQuantificationSet of interest.
        :rtype: :class:`ga4gh.protocol.RnaQuantificationSet`
        """
        return self._run_get_request(
            "rnaquantificationsets", protocol.RnaQuantificationSet,
            rna_quantification_set_id)

    def get_rna_quantification(self, rna_quantification_id):
        """
        Returns the RnaQuantification with the specified ID from the server.
        :param str rna_quantification_id: The ID of the RnaQuantification of
            interest.
        :return: The RnaQuantification of interest.
        :rtype: :class:`ga4gh.protocol.RnaQuantification`
        """
        return self._run_get_request(
            "rnaquantifications", protocol.RnaQuantification,
            rna_quantification_id)

    def get_expression_level(self, expression_level_id):
        """
        Returns the ExpressionLevel with the specified ID from the server.
        :param str expression_level_id: The ID of the ExpressionLevel of
            interest.
        :return: The ExpressionLevel of interest.
        :rtype: :class:`ga4gh.protocol.ExpressionLevel`
        """
        return self._run_get_request(
            "expressionlevels", protocol.ExpressionLevel,
            expression_level_id)

    def search_variants(
            self, variant_set_id, start=None, end=None, reference_name=None,
            call_set_ids=None):
        """
        Returns an iterator over the Variants fulfilling the specified
        conditions from the specified VariantSet.

        :param str variant_set_id: The ID of the
            :class:`ga4gh.protocol.VariantSet` of interest.
        :param int start: Required. The beginning of the window (0-based,
            inclusive) for which overlapping variants should be returned.
            Genomic positions are non-negative integers less than reference
            length. Requests spanning the join of circular genomes are
            represented as two requests one on each side of the join
            (position 0).
        :param int end: Required. The end of the window (0-based, exclusive)
            for which overlapping variants should be returned.
        :param str reference_name: The name of the
            :class:`ga4gh.protocol.Reference` we wish to return variants from.
        :param list call_set_ids: Only return variant calls which belong to
            call sets with these IDs. If an empty array, returns variants
            without any call objects. If null, returns all variant calls.

        :return: An iterator over the :class:`ga4gh.protocol.Variant` objects
            defined by the query parameters.
        :rtype: iter
        """
        request = protocol.SearchVariantsRequest()
        request.reference_name = pb.string(reference_name)
        request.start = pb.int(start)
        request.end = pb.int(end)
        request.variant_set_id = variant_set_id
        request.call_set_ids.extend(pb.string(call_set_ids))
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "variants", protocol.SearchVariantsResponse)

    def search_variant_annotations(
            self, variant_annotation_set_id, reference_name="",
            reference_id="", start=0, end=0, effects=[]):
        """
        Returns an iterator over the Variant Annotations fulfilling
        the specified conditions from the specified VariantSet.

        :param str variant_annotation_set_id: The ID of the
            :class:`ga4gh.protocol.VariantAnnotationSet` of interest.
        :param int start: Required. The beginning of the window (0-based,
            inclusive) for which overlapping variants should be returned.
            Genomic positions are non-negative integers less than reference
            length. Requests spanning the join of circular genomes are
            represented as two requests one on each side of the join
            (position 0).
        :param int end: Required. The end of the window (0-based, exclusive)
            for which overlapping variants should be returned.
        :param str reference_name: The name of the
            :class:`ga4gh.protocol.Reference` we wish to return variants from.

        :return: An iterator over the
            :class:`ga4gh.protocol.VariantAnnotation` objects
            defined by the query parameters.
        :rtype: iter
        """
        request = protocol.SearchVariantAnnotationsRequest()
        request.variant_annotation_set_id = variant_annotation_set_id
        request.reference_name = reference_name
        request.reference_id = reference_id
        request.start = start
        request.end = end
        for effect in effects:
            request.effects.add().CopyFrom(protocol.OntologyTerm(**effect))
        for effect in request.effects:
            if not effect.id:
                raise exceptions.BadRequestException(
                    "Each ontology term should have an id set")
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "variantannotations",
            protocol.SearchVariantAnnotationsResponse)

    def search_features(
            self, feature_set_id=None, parent_id="", reference_name="",
            start=0, end=0, feature_types=[], name="", gene_symbol=""):
        """
        Returns the result of running a search_features method
        on a request with the passed-in parameters.

        :param str feature_set_id: ID of the feature Set being searched
        :param str parent_id: ID (optional) of the parent feature
        :param str reference_name: name of the reference to search
            (ex: "chr1")
        :param int start: search start position on reference
        :param int end: end position on reference
        :param feature_types: array of terms to limit search by (ex: "gene")
        :param str name: only return features with this name
        :param str gene_symbol: only return features on this gene
        :return: an iterator over Features as returned in the
            SearchFeaturesResponse object.
        """
        request = protocol.SearchFeaturesRequest()
        request.feature_set_id = feature_set_id
        request.parent_id = parent_id
        request.reference_name = reference_name
        request.name = name
        request.gene_symbol = gene_symbol
        request.start = start
        request.end = end
        request.feature_types.extend(feature_types)
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "features",
            protocol.SearchFeaturesResponse)

    def search_datasets(self):
        """
        Returns an iterator over the Datasets on the server.

        :return: An iterator over the :class:`ga4gh.protocol.Dataset`
            objects on the server.
        """
        request = protocol.SearchDatasetsRequest()
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "datasets", protocol.SearchDatasetsResponse)

    def search_variant_sets(self, dataset_id):
        """
        Returns an iterator over the VariantSets fulfilling the specified
        conditions from the specified Dataset.

        :param str dataset_id: The ID of the :class:`ga4gh.protocol.Dataset`
            of interest.
        :return: An iterator over the :class:`ga4gh.protocol.VariantSet`
            objects defined by the query parameters.
        """
        request = protocol.SearchVariantSetsRequest()
        request.dataset_id = dataset_id
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "variantsets", protocol.SearchVariantSetsResponse)

    def search_variant_annotation_sets(self, variant_set_id):
        """
        Returns an iterator over the Annotation Sets fulfilling the specified
        conditions from the specified variant set.

        :param str variant_set_id: The ID of the
            :class:`ga4gh.protocol.VariantSet` of interest.
        :return: An iterator over the :class:`ga4gh.protocol.AnnotationSet`
            objects defined by the query parameters.
        """
        request = protocol.SearchVariantAnnotationSetsRequest()
        request.variant_set_id = variant_set_id
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "variantannotationsets",
            protocol.SearchVariantAnnotationSetsResponse)

    def search_feature_sets(self, dataset_id):
        """
        Returns an iterator over the FeatureSets fulfilling the specified
        conditions from the specified Dataset.

        :param str dataset_id: The ID of the
            :class:`ga4gh.protocol.Dataset` of interest.
        :return: An iterator over the :class:`ga4gh.protocol.FeatureSet`
            objects defined by the query parameters.
        """
        request = protocol.SearchFeatureSetsRequest()
        request.dataset_id = dataset_id
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "featuresets", protocol.SearchFeatureSetsResponse)

    def search_reference_sets(
            self, accession=None, md5checksum=None, assembly_id=None):
        """
        Returns an iterator over the ReferenceSets fulfilling the specified
        conditions.

        :param str accession: If not null, return the reference sets for which
            the `accession` matches this string (case-sensitive, exact match).
        :param str md5checksum: If not null, return the reference sets for
            which the `md5checksum` matches this string (case-sensitive, exact
            match). See :class:`ga4gh.protocol.ReferenceSet::md5checksum` for
            details.
        :param str assembly_id: If not null, return the reference sets for
            which the `assembly_id` matches this string (case-sensitive,
            exact match).
        :return: An iterator over the :class:`ga4gh.protocol.ReferenceSet`
            objects defined by the query parameters.
        """
        request = protocol.SearchReferenceSetsRequest()
        request.accession = pb.string(accession)
        request.md5checksum = pb.string(md5checksum)
        request.assembly_id = pb.string(assembly_id)
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "referencesets", protocol.SearchReferenceSetsResponse)

    def search_references(
            self, reference_set_id, accession=None, md5checksum=None):
        """
        Returns an iterator over the References fulfilling the specified
        conditions from the specified Dataset.

        :param str reference_set_id: The ReferenceSet to search.
        :param str accession: If not None, return the references for which the
            `accession` matches this string (case-sensitive, exact match).
        :param str md5checksum: If not None, return the references for which
            the `md5checksum` matches this string (case-sensitive, exact
            match).
        :return: An iterator over the :class:`ga4gh.protocol.Reference`
            objects defined by the query parameters.
        """
        request = protocol.SearchReferencesRequest()
        request.reference_set_id = reference_set_id
        request.accession = pb.string(accession)
        request.md5checksum = pb.string(md5checksum)
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "references", protocol.SearchReferencesResponse)

    def search_call_sets(self, variant_set_id, name=None, bio_sample_id=None):
        """
        Returns an iterator over the CallSets fulfilling the specified
        conditions from the specified VariantSet.

        :param str variant_set_id: Find callsets belonging to the
            provided variant set.
        :param str name: Only CallSets matching the specified name will
            be returned.
        :param str bio_sample_id: Only CallSets matching this id will
            be returned.
        :return: An iterator over the :class:`ga4gh.protocol.CallSet`
            objects defined by the query parameters.
        """
        request = protocol.SearchCallSetsRequest()
        request.variant_set_id = variant_set_id
        request.name = pb.string(name)
        request.bio_sample_id = pb.string(bio_sample_id)
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "callsets", protocol.SearchCallSetsResponse)

    def search_bio_samples(self, dataset_id, name=None, individual_id=None):
        """
        Returns an iterator over the BioSamples fulfilling the specified
        conditions.

        :param str dataset_id: The dataset to search within.
        :param str name: Only BioSamples matching the specified name will
            be returned.
        :param str individual_id: Only BioSamples matching matching this
            id will be returned.
        :return: An iterator over the :class:`ga4gh.protocol.BioSample`
            objects defined by the query parameters.
        """
        request = protocol.SearchBioSamplesRequest()
        request.dataset_id = dataset_id
        request.name = pb.string(name)
        request.individual_id = pb.string(individual_id)
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "biosamples", protocol.SearchBioSamplesResponse)

    def search_individuals(self, dataset_id, name=None):
        """
        Returns an iterator over the Individuals fulfilling the specified
        conditions.

        :param str dataset_id: The dataset to search within.
        :param str name: Only Individuals matching the specified name will
            be returned.
        :return: An iterator over the :class:`ga4gh.protocol.BioSample`
            objects defined by the query parameters.
        """
        request = protocol.SearchIndividualsRequest()
        request.dataset_id = dataset_id
        request.name = pb.string(name)
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "individuals", protocol.SearchIndividualsResponse)

    def search_read_group_sets(
            self, dataset_id, name=None, bio_sample_id=None):
        """
        Returns an iterator over the ReadGroupSets fulfilling the specified
        conditions from the specified Dataset.

        :param str name: Only ReadGroupSets matching the specified name
            will be returned.
        :param str bio_sample_id: Only ReadGroups matching the specified
            bioSample will be included in the response.
        :return: An iterator over the :class:`ga4gh.protocol.ReadGroupSet`
            objects defined by the query parameters.
        :rtype: iter
        """
        request = protocol.SearchReadGroupSetsRequest()
        request.dataset_id = dataset_id
        request.name = pb.string(name)
        request.bio_sample_id = pb.string(bio_sample_id)
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "readgroupsets", protocol.SearchReadGroupSetsResponse)

    def search_reads(
            self, read_group_ids, reference_id=None, start=None, end=None):
        """
        Returns an iterator over the Reads fulfilling the specified
        conditions from the specified read_group_ids.

        :param str read_group_ids: The IDs of the
            :class:`ga4gh.protocol.ReadGroup` of interest.
        :param str reference_id: The name of the
            :class:`ga4gh.protocol.Reference` we wish to return reads
            mapped to.
        :param int start: The start position (0-based) of this query. If a
            reference is specified, this defaults to 0. Genomic positions are
            non-negative integers less than reference length. Requests spanning
            the join of circular genomes are represented as two requests one on
            each side of the join (position 0).
        :param int end: The end position (0-based, exclusive) of this query.
            If a reference is specified, this defaults to the reference's
            length.
        :return: An iterator over the
            :class:`ga4gh.protocol.ReadAlignment` objects defined by
            the query parameters.
        :rtype: iter
        """
        request = protocol.SearchReadsRequest()
        request.read_group_ids.extend(read_group_ids)
        request.reference_id = pb.string(reference_id)
        request.start = pb.int(start)
        request.end = pb.int(end)
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "reads", protocol.SearchReadsResponse)

    def search_phenotype_association_sets(self, dataset_id):
        """
        Returns an iterator over the PhenotypeAssociationSets on the server.
        """
        request = protocol.SearchPhenotypeAssociationSetsRequest()
        request.dataset_id = dataset_id
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "phenotypeassociationsets",
            protocol.SearchPhenotypeAssociationSetsResponse)

    def search_genotype_phenotype(
            self, phenotype_association_set_id=None, feature_ids=None,
            phenotype_ids=None, evidence=None):
        """
        Returns an iterator over the GeneotypePhenotype associations from
        the server
        """
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = phenotype_association_set_id
        if feature_ids:
            request.feature_ids.extend(feature_ids)
        if phenotype_ids:
            request.phenotype_ids.extend(phenotype_ids)
        if evidence:
            request.evidence.extend(evidence)
        request.page_size = pb.int(self._page_size)
        self._logger.debug("search_genotype_phenotype {}".format(request))
        return self._run_search_request(
            request, "featurephenotypeassociations",
            protocol.SearchGenotypePhenotypeResponse)

    def search_phenotype(
            self, phenotype_association_set_id=None, phenotype_id=None,
            description=None, type_=None, age_of_onset=None):
        """
        Returns an iterator over the Phenotypes from the server
        """
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = phenotype_association_set_id
        if phenotype_id:
            request.id = phenotype_id
        if description:
            request.description = description
        if type_:
            request.type.mergeFrom(type_)
        if age_of_onset:
            request.age_of_onset = age_of_onset
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "phenotypes",
            protocol.SearchPhenotypesResponse)

    def search_rna_quantification_sets(self, dataset_id):
        """
        Returns an iterator over the RnaQuantificationSet objects from the
        server
        """
        request = protocol.SearchRnaQuantificationSetsRequest()
        request.dataset_id = dataset_id
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "rnaquantificationsets",
            protocol.SearchRnaQuantificationSetsResponse)

    def search_rna_quantifications(
            self, rna_quantification_set_id="", bio_sample_id=""):
        """
        Returns an iterator over the RnaQuantification objects from the server

        :param str rna_quantification_set_id: The ID of the
            :class:`ga4gh.protocol.RnaQuantificationSet` of interest.
        """
        request = protocol.SearchRnaQuantificationsRequest()
        request.rna_quantification_set_id = rna_quantification_set_id
        if bio_sample_id:
            request.bio_sample_id = bio_sample_id
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "rnaquantifications",
            protocol.SearchRnaQuantificationsResponse)

    def search_expression_levels(
            self, rna_quantification_id="", feature_ids=[], threshold=0.0):
        """
        Returns an iterator over the ExpressionLevel objects from the server

        :param str feature_ids: The IDs of the
            :class:`ga4gh.protocol.Feature` of interest.
        :param str rna_quantification_id: The ID of the
            :class:`ga4gh.protocol.RnaQuantification` of interest.
        :param float threshold: Minimum expression of responses to return.
        """
        request = protocol.SearchExpressionLevelsRequest()
        request.rna_quantification_id = rna_quantification_id
        request.feature_ids.extend(feature_ids)
        request.threshold = threshold
        request.page_size = pb.int(self._page_size)
        return self._run_search_request(
            request, "expressionlevels",
            protocol.SearchExpressionLevelsResponse)


class HttpClient(AbstractClient):
    """
    The GA4GH HTTP client. This class provides methods corresponding to the
    GA4GH search and object GET methods.

    .. todo:: Add a better description of the role of this class and include
        links to the high-level API documention.

    :param str urlPrefix: The base URL of the GA4GH server we wish to
        communicate with. This should include the 'http' or 'https' prefix.
    :param int logLevel: The amount of debugging information to log using
        the :mod:`logging` module. This is :data:`logging.WARNING` by default.
    :param str authentication_key: The authentication key provided by the
        server after logging in.
    """

    def __init__(
            self, url_prefix, logLevel=logging.WARNING,
            authentication_key=None):
        super(HttpClient, self).__init__(logLevel)
        self._url_prefix = url_prefix
        self._authentication_key = authentication_key
        self._session = requests.Session()
        self._setup_http_session()
        requests_log = logging.getLogger("requests.packages.urllib3")
        requests_log.setLevel(logLevel)
        requests_log.propagate = True

    def _setup_http_session(self):
        """
        Sets up the common HTTP session parameters used by requests.
        """
        headers = {"Content-type": "application/json"}
        self._session.headers.update(headers)
        # TODO is this unsafe????
        self._session.verify = False

    def _check_response_status(self, response):
        """
        Checks the speficied HTTP response from the requests package and
        raises an exception if a non-200 HTTP code was returned by the
        server.
        """
        if response.status_code != requests.codes.ok:
            self._logger.error("%s %s", response.status_code, response.text)
            raise exceptions.RequestNonSuccessException(
                "Url {0} had status_code {1}".format(
                    response.url, response.status_code))

    def _get_http_parameters(self):
        """
        Returns the basic HTTP parameters we need all requests.
        """
        return {'key': self._authentication_key}

    def _run_search_page_request(
            self, protocol_request, object_name, protocol_response_class):
        url = posixpath.join(self._url_prefix, object_name + '/search')
        data = protocol.toJson(protocol_request)
        self._logger.debug("request:{}".format(data))
        response = self._session.post(
            url, params=self._get_http_parameters(), data=data)
        self._check_response_status(response)
        return self._deserialize_response(
            response.text, protocol_response_class)

    def _run_get_request(self, object_name, protocol_response_class, id_):
        url_suffix = "{object_name}/{id}".format(
            object_name=object_name, id=id_)
        url = posixpath.join(self._url_prefix, url_suffix)
        response = self._session.get(url, params=self._get_http_parameters())
        self._check_response_status(response)
        return self._deserialize_response(
            response.text, protocol_response_class)

    def _run_list_reference_bases_page_request(self, request):
        url_suffix = "listreferencebases"
        url = posixpath.join(self._url_prefix, url_suffix)
        response = self._session.post(
            url, params=self._get_http_parameters(),
            data=protocol.toJson(request))
        self._check_response_status(response)
        return self._deserialize_response(
            response.text, protocol.ListReferenceBasesResponse)


class LocalClient(AbstractClient):

    def __init__(self, backend):
        super(LocalClient, self).__init__()
        self._backend = backend
        self._get_method_map = {
            "callsets": self._backend.runGetCallSet,
            "datasets": self._backend.runGetDataset,
            "referencesets": self._backend.runGetReferenceSet,
            "references": self._backend.runGetReference,
            "variantsets": self._backend.runGetVariantSet,
            "featuresets": self._backend.runGetFeatureSet,
            "variants": self._backend.runGetVariant,
            "features": self._backend.runGetFeature,
            "readgroupsets": self._backend.runGetReadGroupSet,
            "readgroups": self._backend.runGetReadGroup,
            "variantannotationsets": self._backend.runGetVariantAnnotationSet,
            "biosamples": self._backend.runGetBioSample,
            "individuals": self._backend.runGetIndividual,
            "rnaquantificationsets": self._backend.runGetRnaQuantificationSet,
            "rnaquantifications": self._backend.runGetRnaQuantification,
            "expressionlevels": self._backend.runGetExpressionLevel,
        }
        self._search_method_map = {
            "callsets": self._backend.runSearchCallSets,
            "datasets": self._backend.runSearchDatasets,
            "referencesets": self._backend.runSearchReferenceSets,
            "references": self._backend.runSearchReferences,
            "variantsets": self._backend.runSearchVariantSets,
            "featuresets": self._backend.runSearchFeatureSets,
            "variants": self._backend.runSearchVariants,
            "features": self._backend.runSearchFeatures,
            "readgroupsets": self._backend.runSearchReadGroupSets,
            "reads": self._backend.runSearchReads,
            "variantannotations": self._backend.runSearchVariantAnnotations,
            "variantannotationsets":
                self._backend.runSearchVariantAnnotationSets,
            "biosamples": self._backend.runSearchBioSamples,
            "individuals": self._backend.runSearchIndividuals,
            "featurephenotypeassociations":
                self._backend.runSearchGenotypePhenotypes,
            "phenotypes": self._backend.runSearchPhenotypes,
            "phenotypeassociationsets":
                self._backend.runSearchPhenotypeAssociationSets,
            "rnaquantificationsets":
                self._backend.runSearchRnaQuantificationSets,
            "rnaquantifications": self._backend.runSearchRnaQuantifications,
            "expressionlevels": self._backend.runSearchExpressionLevels,
        }

    def _run_get_request(self, object_name, protocol_response_class, id_):
        get_method = self._get_method_map[object_name]
        response_json = get_method(id_)
        return self._deserialize_response(
            response_json, protocol_response_class)

    def _run_search_page_request(
            self, protocol_request, object_name, protocol_response_class):
        search_method = self._search_method_map[object_name]
        response_json = search_method(protocol.toJson(protocol_request))
        return self._deserialize_response(
            response_json, protocol_response_class)

    def _run_list_reference_bases_page_request(self, request):
        response_json = self._backend.runListReferenceBases(
            protocol.toJson(request))
        return self._deserialize_response(
            response_json, protocol.ListReferenceBasesResponse)
