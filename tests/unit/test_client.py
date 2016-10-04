"""
Tests for the client
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import mock
import json

import ga4gh.protocol as protocol
import ga4gh.backend as backend
import ga4gh.client as client
import ga4gh.datarepo as datarepo
import tests.utils as utils
import ga4gh.exceptions as exceptions


class TestSearchMethodsCallRunRequest(unittest.TestCase):
    """
    Test that search methods call lower-level functionality correctly
    """
    def setUp(self):
        self.httpClient = client.HttpClient("http://example.com")
        self.httpClient._run_search_request = mock.Mock()
        self.httpClient._run_get_request = mock.Mock()
        self.objectId = "SomeId"
        self.objectName = "objectName"
        self.datasetId = "datasetId"
        self.variantSetId = "variantSetId"
        self.variantAnnotationSetId = "variantAnnotationSetId"
        self.featureSetId = "featureSetId"
        self.parentId = "parentId"
        self.feature = "feature"
        self.referenceSetId = "referenceSetId"
        self.referenceId = "referenceId"
        self.readGroupIds = ["readGroupId"]
        self.referenceName = "referenceName"
        self.bioSampleId = "bioSampleId"
        self.bioSampleName = "bioSampleName"
        self.individualName = "individualName"
        self.individualId = "individualId"
        self.geneSymbol = "geneSymbol"
        self.start = 100
        self.end = 101
        self.referenceName = "referenceName"
        self.callSetIds = ["id1", "id2"]
        self.pageSize = 1000
        self.httpClient.set_page_size(self.pageSize)
        self.assemblyId = "assemblyId"
        self.accession = "accession"
        self.md5checksum = "md5checksum"
        self.phenotype_association_set_id = "phenotype_association_set_id"
        self.feature_ids = ["id1", "id2"]
        self.phenotype_ids = ["id3", "id4"]
        self.evidence = protocol.EvidenceQuery()
        self.rnaQuantificationSetId = "rnaQuantificationSetId"
        self.rnaQuantificationId = "rnaQuantificationId"
        self.expressionLevelId = "expressionLevelId"
        self.threshold = 0.0

    def testSetPageSize(self):
        testClient = client.AbstractClient()
        # pageSize is None by default
        self.assertIsNone(testClient.get_page_size())
        for pageSize in [1, 10, 100]:
            testClient.set_page_size(pageSize)
            self.assertEqual(testClient.get_page_size(), pageSize)

    def testSearchVariants(self):
        request = protocol.SearchVariantsRequest()
        request.reference_name = self.referenceName
        request.start = self.start
        request.end = self.end
        request.variant_set_id = self.variantSetId
        request.call_set_ids.extend(self.callSetIds)
        request.page_size = self.pageSize
        self.httpClient.search_variants(
            self.variantSetId, start=self.start, end=self.end,
            reference_name=self.referenceName, call_set_ids=self.callSetIds)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "variants", protocol.SearchVariantsResponse)

    def testSearchDatasets(self):
        request = protocol.SearchDatasetsRequest()
        request.page_size = self.pageSize
        self.httpClient.search_datasets()
        self.httpClient._run_search_request.assert_called_once_with(
            request, "datasets", protocol.SearchDatasetsResponse)

    def testSearchVariantSets(self):
        request = protocol.SearchVariantSetsRequest()
        request.dataset_id = self.datasetId
        request.page_size = self.pageSize
        self.httpClient.search_variant_sets(self.datasetId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "variantsets", protocol.SearchVariantSetsResponse)

    def testSearchVariantAnnotationSets(self):
        request = protocol.SearchVariantAnnotationSetsRequest()
        request.variant_set_id = self.variantSetId
        request.page_size = self.pageSize
        self.httpClient.search_variant_annotation_sets(self.variantSetId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "variantannotationsets",
            protocol.SearchVariantAnnotationSetsResponse)

    def testSearchVariantAnnotations(self):
        request = protocol.SearchVariantAnnotationsRequest()
        request.variant_annotation_set_id = self.variantAnnotationSetId
        request.page_size = self.pageSize
        request.reference_name = self.referenceName
        request.reference_id = self.referenceId
        request.start = self.start
        request.end = self.end
        self.httpClient.search_variant_annotations(
            self.variantAnnotationSetId,
            reference_name=self.referenceName,
            start=self.start,
            end=self.end,
            effects=[],
            reference_id=self.referenceId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "variantannotations",
            protocol.SearchVariantAnnotationsResponse)
        with self.assertRaises(exceptions.BadRequestException):
            self.httpClient.search_variant_annotations(
                self.variantAnnotationSetId,
                reference_name=self.referenceName,
                start=self.start,
                end=self.end,
                effects=[{"term": "just a term"}, {"id": "an id"}],
                reference_id=self.referenceId)

    def testSearchFeatures(self):
        request = protocol.SearchFeaturesRequest()
        request.feature_set_id = self.featureSetId
        request.parent_id = self.parentId
        request.page_size = self.pageSize
        request.reference_name = self.referenceName
        request.start = self.start
        request.end = self.end
        request.name = self.objectName
        request.gene_symbol = self.geneSymbol
        request.feature_types.append(self.feature)
        self.httpClient.search_features(
            self.featureSetId, parent_id=self.parentId,
            reference_name=self.referenceName, start=self.start,
            end=self.end, feature_types=[self.feature],
            name=self.objectName, gene_symbol=self.geneSymbol)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "features", protocol.SearchFeaturesResponse)

    def testSearchFeatureSets(self):
        request = protocol.SearchFeatureSetsRequest()
        request.dataset_id = self.datasetId
        request.page_size = self.pageSize
        self.httpClient.search_feature_sets(self.datasetId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "featuresets", protocol.SearchFeatureSetsResponse)

    def testSearchReferenceSets(self):
        request = protocol.SearchReferenceSetsRequest()
        request.page_size = self.pageSize
        request.accession = self.accession
        request.md5checksum = self.md5checksum
        request.assembly_id = self.assemblyId
        self.httpClient.search_reference_sets(
            accession=self.accession, md5checksum=self.md5checksum,
            assembly_id=self.assemblyId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "referencesets", protocol.SearchReferenceSetsResponse)

    def testSearchReferences(self):
        request = protocol.SearchReferencesRequest()
        request.reference_set_id = self.referenceSetId
        request.page_size = self.pageSize
        request.accession = self.accession
        request.md5checksum = self.md5checksum
        self.httpClient.search_references(
            self.referenceSetId, accession=self.accession,
            md5checksum=self.md5checksum)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "references", protocol.SearchReferencesResponse)

    def testSearchReadGroupSets(self):
        request = protocol.SearchReadGroupSetsRequest()
        request.dataset_id = self.datasetId
        request.name = self.objectName
        request.bio_sample_id = self.bioSampleId
        request.page_size = self.pageSize
        self.httpClient.search_read_group_sets(
            self.datasetId,
            name=self.objectName,
            bio_sample_id=self.bioSampleId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "readgroupsets", protocol.SearchReadGroupSetsResponse)

    def testSearchCallSets(self):
        request = protocol.SearchCallSetsRequest()
        request.variant_set_id = self.variantSetId
        request.name = self.objectName
        request.bio_sample_id = self.bioSampleId
        request.page_size = self.pageSize
        self.httpClient.search_call_sets(
            self.variantSetId,
            name=self.objectName,
            bio_sample_id=self.bioSampleId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "callsets", protocol.SearchCallSetsResponse)

    def testSearchReads(self):
        request = protocol.SearchReadsRequest()
        request.read_group_ids.extend(self.readGroupIds)
        request.reference_id = self.referenceId
        request.start = self.start
        request.end = self.end
        request.page_size = self.pageSize
        self.httpClient.search_reads(
            self.readGroupIds, reference_id=self.referenceId,
            start=self.start, end=self.end)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "reads", protocol.SearchReadsResponse)

    def testSearchExpressionLevels(self):
        request = protocol.SearchExpressionLevelsRequest()
        request.rna_quantification_id = self.rnaQuantificationId
        request.threshold = self.threshold
        request.page_size = self.pageSize
        self.httpClient.search_expression_levels(
            self.rnaQuantificationId, threshold=self.threshold)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "expressionlevels",
            protocol.SearchExpressionLevelsResponse)

    def testSearchRnaQuantificationSets(self):
        request = protocol.SearchRnaQuantificationSetsRequest()
        request.dataset_id = self.datasetId
        request.page_size = self.pageSize
        self.httpClient.search_rna_quantification_sets(self.datasetId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "rnaquantificationsets",
            protocol.SearchRnaQuantificationSetsResponse)

    def testSearchRnaQuantifications(self):
        request = protocol.SearchRnaQuantificationsRequest()
        request.rna_quantification_set_id = self.rnaQuantificationSetId
        request.page_size = self.pageSize
        self.httpClient.search_rna_quantifications(
            rna_quantification_set_id=self.rnaQuantificationSetId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "rnaquantifications",
            protocol.SearchRnaQuantificationsResponse)

    def testSearchBioSamples(self):
        request = protocol.SearchBioSamplesRequest()
        request.dataset_id = self.datasetId
        request.name = self.bioSampleName
        request.individual_id = self.individualId
        request.page_size = self.pageSize
        self.httpClient.search_bio_samples(
            self.datasetId, self.bioSampleName, self.individualId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "biosamples", protocol.SearchBioSamplesResponse)

    def testSearchIndividuals(self):
        request = protocol.SearchIndividualsRequest()
        request.dataset_id = self.datasetId
        request.name = self.individualName
        request.page_size = self.pageSize
        self.httpClient.search_individuals(
            self.datasetId, self.individualName)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "individuals", protocol.SearchIndividualsResponse)

    def testGetReferenceSet(self):
        self.httpClient.get_reference_set(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "referencesets", protocol.ReferenceSet, self.objectId)

    def testGetVariantAnnotationSet(self):
        self.httpClient.get_variant_annotation_set(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "variantannotationsets", protocol.VariantAnnotationSet,
            self.objectId)

    def testGetVariantSet(self):
        self.httpClient.get_variant_set(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "variantsets", protocol.VariantSet, self.objectId)

    def testGetReference(self):
        self.httpClient.get_reference(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "references", protocol.Reference, self.objectId)

    def testGetReadGroupSets(self):
        self.httpClient.get_read_group_set(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "readgroupsets", protocol.ReadGroupSet, self.objectId)

    def testGetReadGroup(self):
        self.httpClient.get_read_group(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "readgroups", protocol.ReadGroup, self.objectId)

    def testGetCallSets(self):
        self.httpClient.get_call_set(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "callsets", protocol.CallSet, self.objectId)

    def testGetDatasets(self):
        self.httpClient.get_dataset(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "datasets", protocol.Dataset, self.objectId)

    def testGetVariant(self):
        self.httpClient.get_variant(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "variants", protocol.Variant, self.objectId)

    def testGetBioSample(self):
        self.httpClient.get_bio_sample(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "biosamples", protocol.BioSample, self.objectId)

    def testGetIndividual(self):
        self.httpClient.get_individual(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "individuals", protocol.Individual, self.objectId)

    def testGetRnaQuantificationSet(self):
        self.httpClient.get_rna_quantification_set(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "rnaquantificationsets", protocol.RnaQuantificationSet,
            self.objectId)

    def testGetRnaQuantification(self):
        self.httpClient.get_rna_quantification(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "rnaquantifications", protocol.RnaQuantification, self.objectId)

    def testGetExpressionLevel(self):
        self.httpClient.get_expression_level(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "expressionlevels", protocol.ExpressionLevel, self.objectId)

    def testGetFeatureSet(self):
        self.httpClient.get_feature_set(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "featuresets", protocol.FeatureSet, self.objectId)

    def testGetFeature(self):
        self.httpClient.get_feature(self.objectId)
        self.httpClient._run_get_request.assert_called_once_with(
            "features", protocol.Feature, self.objectId)

    def testSearchGenotypePhenotype(self):
        request = protocol.SearchGenotypePhenotypeRequest()
        request.phenotype_association_set_id = \
            self.phenotype_association_set_id
        request.feature_ids.extend(self.feature_ids)
        request.phenotype_ids.extend(self.phenotype_ids)
        request.evidence.extend([self.evidence])
        request.page_size = self.pageSize
        self.httpClient.search_genotype_phenotype(
            phenotype_association_set_id=self.phenotype_association_set_id,
            feature_ids=self.feature_ids,
            phenotype_ids=self.phenotype_ids,
            evidence=[self.evidence])
        self.httpClient._run_search_request.assert_called_once_with(
            request, "featurephenotypeassociations",
            protocol.SearchGenotypePhenotypeResponse)

    def testSearchPhenotype(self):
        request = protocol.SearchPhenotypesRequest()
        request.phenotype_association_set_id = \
            self.phenotype_association_set_id
        request.id = self.phenotype_ids[0]
        request.page_size = self.pageSize
        self.httpClient.search_phenotype(
            phenotype_association_set_id=self.phenotype_association_set_id,
            phenotype_id=self.phenotype_ids[0])
        self.httpClient._run_search_request.assert_called_once_with(
            request, "phenotypes",
            protocol.SearchPhenotypesResponse)

    def testSearchPhenotypeAssociationSets(self):
        request = protocol.SearchPhenotypeAssociationSetsRequest()
        request.dataset_id = self.datasetId
        request.page_size = self.pageSize
        self.httpClient.search_phenotype_association_sets(
            dataset_id=self.datasetId)
        self.httpClient._run_search_request.assert_called_once_with(
            request, "phenotypeassociationsets",
            protocol.SearchPhenotypeAssociationSetsResponse)


class DatamodelObjectWrapper(object):
    """
    Thin wrapper class that allows us to treat data model objects uniformly.
    We should update the data model interface so that objects are always
    returned so that we always call toProtocolElement on the results.
    Variants and Reads are the exceptions here.
    """
    def __init__(self, gaObject):
        self.gaObject = gaObject

    def toProtocolElement(self):
        return self.gaObject


class DummyResponse(object):
    """
    Stand in for requests Response object;
    """
    def __init__(self, text):
        self.text = text
        self.status_code = 200


class DummyRequestsSession(object):
    """
    Takes the place of a requests session so that we can check that all
    values are sent and received correctly.
    """
    def __init__(self, backend, urlPrefix):
        self._backend = backend
        self._urlPrefix = urlPrefix
        self._getMethodMap = {
            "datasets": self._backend.runGetDataset,
            "referencesets": self._backend.runGetReferenceSet,
            "references": self._backend.runGetReference,
            "variantsets": self._backend.runGetVariantSet,
            "variants": self._backend.runGetVariant,
            "readgroupsets": self._backend.runGetReadGroupSet,
            "readgroups": self._backend.runGetReadGroup,
            "rnaquantifications": self._backend.runGetRnaQuantification,
            "rnaquantificationsets": self._backend.runGetRnaQuantificationSet,
            "expressionlevels": self._backend.runGetExpressionLevel,
        }
        self._searchMethodMap = {
            "datasets": self._backend.runSearchDatasets,
            "referencesets": self._backend.runSearchReferenceSets,
            "references": self._backend.runSearchReferences,
            "variantsets": self._backend.runSearchVariantSets,
            "variants": self._backend.runSearchVariants,
            "readgroupsets": self._backend.runSearchReadGroupSets,
            "reads": self._backend.runSearchReads,
            "rnaquantifications": self._backend.runSearchRnaQuantifications,
            "rnaquantificationsets":
                self._backend.runSearchRnaQuantificationSets,
            "expressionlevels": self._backend.runSearchExpressionLevels,
        }
        self.headers = {}

    def checkSessionParameters(self):
        contentType = "Content-type"
        assert contentType in self.headers
        assert self.headers[contentType] == "application/json"

    def get(self, url, params):
        # TODO add some more checks for params to see if Key is set,
        # and we're not sending any extra stuff.
        self.checkSessionParameters()
        assert url.startswith(self._urlPrefix)
        suffix = url[len(self._urlPrefix):]
        splits = suffix.split("/")
        assert len(splits) == 3
        assert splits[0] == ''
        datatype, id_ = splits[1:]
        assert datatype in self._getMethodMap
        method = self._getMethodMap[datatype]
        result = method(id_)
        return DummyResponse(result)

    def post(self, url, params=None, data=None):
        self.checkSessionParameters()
        assert url.startswith(self._urlPrefix)
        suffix = url[len(self._urlPrefix):]
        searchSuffix = "/search"
        if suffix.endswith(searchSuffix):
            datatype = suffix[1:-len(searchSuffix)]
            assert datatype in self._searchMethodMap
            method = self._searchMethodMap[datatype]
            result = method(data)
        else:
            # ListReferenceBases is an oddball and needs to be treated
            # separately.
            data = json.loads(data)
            args = protocol.ListReferenceBasesRequest()
            args.reference_id = data.get('referenceId', "")
            args.start = int(data.get('start', 0))
            args.end = int(data.get('end', 0))
            args.page_token = data.get('pageToken', "")
            result = self._backend.runListReferenceBases(
                protocol.toJson(args))
        return DummyResponse(result)


class DummyHttpClient(client.HttpClient):
    """
    Client in which we intercept calls to the underlying requests connection.
    """
    def __init__(self, backend):
        self._urlPrefix = "http://example.com"
        super(DummyHttpClient, self).__init__(self._urlPrefix)
        self._session = DummyRequestsSession(backend, self._urlPrefix)
        self._setup_http_session()


class ExhaustiveListingsMixin(object):
    """
    Tests exhaustive listings using the high-level API with a Simulated
    backend.
    """
    @classmethod
    def setUpClass(cls):
        cls.backend = backend.Backend(datarepo.SimulatedDataRepository(
            randomSeed=100, numDatasets=3,
            numVariantSets=3, numCalls=3, variantDensity=0.5,
            numReferenceSets=3, numReferencesPerReferenceSet=3,
            numReadGroupSets=3, numReadGroupsPerReadGroupSet=3,
            numAlignments=3, numRnaQuantSets=3))
        cls.dataRepo = cls.backend.getDataRepository()

    def setUp(self):
        self.client = self.getClient()

    def verifyObjectList(self, gaObjects, datamodelObjects, getMethod):
        """
        Verifies that the specified list of protocol objects corresponds
        to the specified list of datamodel objects.
        """
        for gaObject, datamodelObject in utils.zipLists(
                gaObjects, datamodelObjects):
            self.assertEqual(gaObject, datamodelObject.toProtocolElement())
            otherGaObject = getMethod(gaObject.id)
            self.assertEqual(gaObject, otherGaObject)

    def testAllDatasets(self):
        datasets = list(self.client.search_datasets())
        self.verifyObjectList(
            datasets, self.dataRepo.getDatasets(), self.client.get_dataset)

    def testAllReferenceSets(self):
        referenceSets = list(self.client.search_reference_sets())
        self.verifyObjectList(
            referenceSets, self.dataRepo.getReferenceSets(),
            self.client.get_reference_set)

    def testAllReferences(self):
        for referenceSet in self.client.search_reference_sets():
            references = list(self.client.search_references(referenceSet.id))
            datamodelReferences = self.dataRepo.getReferenceSet(
                referenceSet.id).getReferences()
            self.verifyObjectList(
                references, datamodelReferences, self.client.get_reference)
            for datamodelReference in datamodelReferences:
                bases = self.client.list_reference_bases(
                    datamodelReference.getId())
                otherBases = datamodelReference.getBases(
                    0, datamodelReference.getLength())
                self.assertEqual(bases, otherBases)

    def testAllVariantSets(self):
        for dataset in self.client.search_datasets():
            variantSets = list(self.client.search_variant_sets(dataset.id))
            datamodelVariantSets = self.dataRepo.getDataset(
                dataset.id).getVariantSets()
            self.verifyObjectList(
                variantSets, datamodelVariantSets, self.client.get_variant_set)

    def testAllVariants(self):
        for datamodelDataset in self.dataRepo.getDatasets():
            for datamodelVariantSet in datamodelDataset.getVariantSets():
                # TODO the values should be derived from the datamodel
                # variant set object.
                start = 0
                end = 20
                referenceName = "fixme"
                variants = list(self.client.search_variants(
                    datamodelVariantSet.getId(), start=start, end=end,
                    reference_name=referenceName))
                datamodelVariants = [
                    DatamodelObjectWrapper(variant) for variant in
                    datamodelVariantSet.getVariants(
                        referenceName, start, end)]
                self.verifyObjectList(
                    variants, datamodelVariants, self.client.get_variant)

    def testAllReadGroupSets(self):
        for dataset in self.client.search_datasets():
            readGroupSets = list(
                self.client.search_read_group_sets(dataset.id))
            datamodelReadGroupSets = self.dataRepo.getDataset(
                dataset.id).getReadGroupSets()
            self.verifyObjectList(
                readGroupSets, datamodelReadGroupSets,
                self.client.get_read_group_set)
            # Check the readGroups.
            for readGroupSet, datamodelReadGroupSet in zip(
                    readGroupSets, datamodelReadGroupSets):
                datamodelReadGroups = datamodelReadGroupSet.getReadGroups()
                self.verifyObjectList(
                    readGroupSet.read_groups, datamodelReadGroups,
                    self.client.get_read_group)

    def testAllReads(self):
        for dmDataset in self.dataRepo.getDatasets():
            for dmReadGroupSet in dmDataset.getReadGroupSets():
                dmReferenceSet = dmReadGroupSet.getReferenceSet()
                for dmReadGroup in dmReadGroupSet.getReadGroups():
                    for dmReference in dmReferenceSet.getReferences():
                        # TODO fix these coordinates.
                        start = 0
                        end = 10
                        dmReads = list(dmReadGroup.getReadAlignments(
                            dmReference, start, end))
                        reads = list(self.client.search_reads(
                            [dmReadGroup.getId()], dmReference.getId(),
                            start, end))
                        self.assertGreater(len(reads), 0)
                        for dmRead, read in utils.zipLists(dmReads, reads):
                            self.assertEqual(dmRead, read)

    def testAllRnaQuantificationSets(self):
        for dataset in self.client.search_datasets():
            rnaQuantificationSets = \
                list(self.client.search_rna_quantification_sets(dataset.id))
            datamodelRnaQuantificationSets = self.dataRepo.getDataset(
                dataset.id).getRnaQuantificationSets()
            self.verifyObjectList(
                rnaQuantificationSets, datamodelRnaQuantificationSets,
                self.client.get_rna_quantification_set)


class TestExhaustiveListingsHttp(ExhaustiveListingsMixin, unittest.TestCase):
    """
    Tests the exhaustive listings using the HTTP client.
    """

    def getClient(self):
        return DummyHttpClient(self.backend)


class TestExhaustiveListingsLocal(ExhaustiveListingsMixin, unittest.TestCase):
    """
    Tests the exhaustive listings using the local client.
    """

    def getClient(self):
        return client.LocalClient(self.backend)


class PagingMixin(object):
    """
    Tests the paging code using a simulated backend.
    """
    @classmethod
    def setUpClass(cls):
        cls.numReferences = 25
        cls.backend = backend.Backend(datarepo.SimulatedDataRepository(
            randomSeed=100, numDatasets=0,
            numReferenceSets=1,
            numReferencesPerReferenceSet=cls.numReferences))
        cls.dataRepo = cls.backend.getDataRepository()

    def setUp(self):
        self.client = self.getClient()
        self.datamodelReferenceSet = self.dataRepo.getReferenceSetByIndex(0)
        self.datamodelReferences = self.datamodelReferenceSet.getReferences()
        self.references = [
            dmReference.toProtocolElement()
            for dmReference in self.datamodelReferences]
        self.assertEqual(len(self.references), self.numReferences)

    def verifyAllReferences(self):
        """
        Verifies that we correctly return all references.
        """
        references = list(self.client.search_references(
            self.datamodelReferenceSet.getId()))
        self.assertEqual(references, self.references)

    def testDefaultPageSize(self):
        self.verifyAllReferences()

    def verifyPageSize(self, pageSize):
        self.client.set_page_size(pageSize)
        self.assertEqual(pageSize, self.client.get_page_size())
        self.verifyAllReferences()

    def testPageSize1(self):
        self.verifyPageSize(1)

    def testPageSize2(self):
        self.verifyPageSize(2)

    def testPageSize3(self):
        self.verifyPageSize(3)

    def testPageSizeAlmostListLength(self):
        self.verifyPageSize(self.numReferences - 1)

    def testPageSizeListLength(self):
        self.verifyPageSize(self.numReferences)


class TestPagingLocal(PagingMixin, unittest.TestCase):
    """
    Tests paging using the local client.
    """

    def getClient(self):
        return client.LocalClient(self.backend)


class TestPagingHttp(PagingMixin, unittest.TestCase):
    """
    Tests paging using the HTTP client.
    """

    def getClient(self):
        return DummyHttpClient(self.backend)
