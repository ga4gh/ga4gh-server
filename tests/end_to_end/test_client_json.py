"""
Performs an end-to-end test where we verify that the data
output by the client command line interface is equal to
the values we expect using the test dataset.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import json
import unittest
import freezegun

import ga4gh.client.client as client
import ga4gh.server.backend as backend
import ga4gh.client.cli as cli_client
import ga4gh.server.datarepo as datarepo
import ga4gh.common.utils as utils
import tests.paths as paths
import server

import ga4gh.schemas.protocol as protocol


def setUpModule():
    global moduleTestServer
    moduleTestServer = server.Ga4ghServerForTestingDataSource(
        paths.testDataRepo)
    moduleTestServer.start()


def tearDownModule():
    moduleTestServer.shutdown()


@freezegun.freeze_time(datetime.datetime.now())
class TestClientOutput(unittest.TestCase):
    """
    Base class for client output tests
    """
    def setUp(self):
        self._maxDiff = None
        repoPath = paths.testDataRepo
        self._dataUrl = moduleTestServer.getUrl()
        dataRepository = datarepo.SqlDataRepository(repoPath)
        dataRepository.open(datarepo.MODE_READ)
        self._backend = backend.Backend(dataRepository)
        self._client = client.LocalClient(self._backend)
        # TODO probably could use a cache of objects, so we don't
        # need to keep re-fetching them

    def captureCliOutput(self, command, arguments, outputFormat):
        clientCommand = "{} {} {} -O {}".format(
            command, self._dataUrl, arguments, outputFormat)
        stdout, stderr = utils.captureOutput(
            cli_client.client_main, clientCommand.split())
        self.assertEqual(len(stderr), 0)
        return stdout


class TestClientFasta(TestClientOutput):
    """
    Tests client FASTA output
    """
    def captureFastaOutput(self, command, arguments=""):
        stdout = self.captureCliOutput(command, arguments, "fasta")
        lines = stdout.split()
        return lines

    def testListReferenceBases(self):
        referenceSetIterator = self._client.search_reference_sets()
        referenceSet = next(referenceSetIterator)
        referencesIterator = self._client.search_references(referenceSet.id)
        reference = next(referencesIterator)
        start = 1
        end = 5
        lines = self.captureFastaOutput(
            "references-list-bases --start {} --end {}".format(start, end),
            reference.id)
        self.assertEqual(
            lines[0], ">{}:{}-{}".format(reference.id, start, end))
        cliBases = ''.join(lines[1:])
        bases = self._client.list_reference_bases(reference.id, start, end)
        self.assertEqual(cliBases, bases)


class TestClientJson(TestClientOutput):
    """
    Tests that the JSON output by the client on the command line for
    various options is equal to the values we find using the Python
    client API.
    """
    def captureJsonOutput(self, command, arguments=""):
        """
        Runs the specified command add the JSON output option and
        returns the result as a list of JSON parsed dictionaries.
        """
        stdout = self.captureCliOutput(command, arguments, "json")
        cliOutput = []
        for line in stdout.splitlines():
            try:
                cliOutput.append(json.loads(line))
            except ValueError, e:
                raise Exception((e, line, stdout, command, arguments))
        return cliOutput

    def verifyParsedOutputsEqual(
            self, clientIterator, cliCommand, cliArguments="",
            scrubFunc=None):
        """
        Verify that the parsed JSON of all the objects in the specified
        client iterator are equal to the parsed JSON from the specified
        CLI command.
        """
        cliOutput = self.captureJsonOutput(cliCommand, cliArguments)
        clientOutput = [
            protocol.toJsonDict(gObj) for gObj in clientIterator]
        if scrubFunc is not None:
            # this is a hack to deal with the issue of the timestamps on
            # objects being different
            scrubFunc(cliOutput, clientOutput)
        self.assertEqual(clientOutput, cliOutput)
        return len(clientOutput)

    def _scrubReadGroup(self, cliReadGroup, clientReadGroup):
        cliReadGroup['created'] = clientReadGroup['created']
        cliReadGroup['updated'] = clientReadGroup['updated']
        cliReadGroup['experiment']['messageUpdateTime'] = \
            clientReadGroup['experiment']['messageUpdateTime']
        cliReadGroup['experiment']['messageCreateTime'] = \
            clientReadGroup['experiment']['messageCreateTime']

    def scrubReadGroup(self, cliOutput, clientOutput):
        for cliReadGroup, clientReadGroup in utils.zipLists(
                cliOutput, clientOutput):
            self._scrubReadGroup(cliReadGroup, clientReadGroup)

    def scrubReadGroupSet(self, cliOutput, clientOutput):
        for cliReadGroupSet, clientReadGroupSet in utils.zipLists(
                cliOutput, clientOutput):
            for cliReadGroup, clientReadGroup in utils.zipLists(
                    cliReadGroupSet['readGroups'],
                    clientReadGroupSet['readGroups']):
                self._scrubReadGroup(cliReadGroup, clientReadGroup)

    def verifyParsedOutputsEqualReadGroup(
            self, clientIterator, cliCommand, cliArguments=""):
        cliOutput = self.captureJsonOutput(cliCommand, cliArguments)
        clientOutput = [
            protocol.toJsonDict(gObj) for gObj in clientIterator]
        self.assertEqual(clientOutput, cliOutput)
        return len(clientOutput)

    def testGetCallSet(self):
        for dataset in self._client.search_datasets():
            for variantSet in self._client.search_variant_sets(dataset.id):
                for callSet in self._client.search_call_sets(variantSet.id):
                    self.verifyParsedOutputsEqual(
                        [callSet], "callsets-get", callSet.id)

    def testGetDataset(self):
        for dataset in self._client.search_datasets():
            self.verifyParsedOutputsEqual(
                [dataset], "datasets-get", dataset.id)

    def testGetReadGroup(self):
        for dataset in self._client.search_datasets():
            for readGroupSet in self._client.search_read_group_sets(
                    dataset.id):
                for readGroup in readGroupSet.read_groups:
                    self.verifyParsedOutputsEqual(
                        [readGroup], "readgroups-get", readGroup.id,
                        scrubFunc=self.scrubReadGroup)

    def testGetReadGroupSet(self):
        for dataset in self._client.search_datasets():
            for readGroupSet in self._client.search_read_group_sets(
                    dataset.id):
                self.verifyParsedOutputsEqual(
                    [readGroupSet], "readgroupsets-get", readGroupSet.id,
                    scrubFunc=self.scrubReadGroupSet)

    def testGetReference(self):
        for referenceSet in self._client.search_reference_sets():
            for reference in self._client.search_references(
                    referenceSet.id):
                self.verifyParsedOutputsEqual(
                    [reference], "references-get", reference.id)

    def testGetReferenceSet(self):
        for referenceSet in self._client.search_reference_sets():
            self.verifyParsedOutputsEqual(
                [referenceSet], "referencesets-get", referenceSet.id)

    @unittest.skip("TODO: clarify semantics of callsets and fix")
    def testGetVariant(self):
        test_executed = 0
        start = 0
        end = 1000
        referenceName = "1"
        for dataset in self._client.search_datasets():
            for variantSet in self._client.search_variant_sets(dataset.id):
                variants = self._client.search_variants(
                    variantSet.id, start=start, end=end,
                    reference_name=referenceName)
                for variant in variants:
                    test_executed += self.verifyParsedOutputsEqual(
                        [variant], "variants-get", variant.id)
        self.assertGreater(test_executed, 0)

    def testGetVariantAnnotationSet(self):
        # TODO this doesn't actually test get_variant_annotation_set
        test_executed = 0
        for dataset in self._client.search_datasets():
            for variantSet in self._client.search_variant_sets(dataset.id):
                for annSet in self._client.search_variant_annotation_sets(
                        variantSet.id):
                    test_executed += self.verifyParsedOutputsEqual(
                        [annSet], "variantannotationsets-get", annSet.id)
        self.assertGreater(test_executed, 0)

    def testGetVariantSet(self):
        for dataset in self._client.search_datasets():
            for variantSet in self._client.search_variant_sets(dataset.id):
                self.verifyParsedOutputsEqual(
                    [variantSet], "variantsets-get", variantSet.id)

    def testSearchCallSets(self):
        for dataset in self._client.search_datasets():
            for variantSet in self._client.search_variant_sets(dataset.id):
                iterator = self._client.search_call_sets(variantSet.id)
                args = "--variantSetId {}".format(variantSet.id)
                self.verifyParsedOutputsEqual(
                    iterator, "callsets-search", args)

    def testGetBiosamples(self):
        for dataset in self._client.search_datasets():
            for biosample in self._client.search_biosamples(dataset.id):
                self.verifyParsedOutputsEqual(
                    [biosample], "biosamples-get", biosample.id)

    def testGetIndividuals(self):
        for dataset in self._client.search_datasets():
            for individual in self._client.search_individuals(dataset.id):
                self.verifyParsedOutputsEqual(
                    [individual], "individuals-get", individual.id)

    def testSearchDatasets(self):
        iterator = self._client.search_datasets()
        self.verifyParsedOutputsEqual(iterator, "datasets-search")

    def testSearchReadGroupSets(self):
        for dataset in self._client.search_datasets():
            iterator = self._client.search_read_group_sets(dataset.id)
            self.verifyParsedOutputsEqual(
                iterator, "readgroupsets-search",
                "--datasetId {}".format(dataset.id),
                scrubFunc=self.scrubReadGroupSet)
            for readGroupSet in iterator:
                bioIterator = self._client.search_read_group_sets(
                    dataset.id,
                    name=readGroupSet.name,
                    biosample_id=readGroupSet.biosampleId)
                self.verifyParsedOutputsEqual(
                    bioIterator, "readgroupsets-search",
                    "--datasetId {} --name {} --biosampleId {}".format(
                        dataset.id,
                        readGroupSet.name,
                        readGroupSet.biosampleId),
                    scrubFunc=self.scrubReadGroupSet)

    def testSearchReads(self):
        test_executed = 0
        start = 0
        end = 1000000
        for dataset in self._client.search_datasets():
            for readGroupSet in self._client.search_read_group_sets(
                    dataset.id):
                for readGroup in readGroupSet.read_groups:
                    reference = self._client.search_references(
                        reference_set_id=readGroup.reference_set_id).next()
                    referenceId = reference.id
                    iterator = self._client.search_reads(
                        [readGroup.id], reference_id=referenceId,
                        start=start, end=end)
                    args = "--start {} --end {} --readGroupIds {}\
                    --referenceId {}".format(
                        start, end, readGroup.id, referenceId)
                    test_executed += self.verifyParsedOutputsEqual(
                        iterator, "reads-search", args)
        self.assertGreater(test_executed, 0)

    def testSearchReferenceSets(self):
        iterator = self._client.search_reference_sets()
        self.verifyParsedOutputsEqual(iterator, "referencesets-search")

    def testSearchReferences(self):
        for referenceSet in self._client.search_reference_sets():
            iterator = self._client.search_references(
                reference_set_id=referenceSet.id)
            args = "--referenceSetId={}".format(referenceSet.id)
            self.verifyParsedOutputsEqual(iterator, "references-search", args)

    def testSearchVariantSets(self):
        for dataset in self._client.search_datasets():
            iterator = self._client.search_variant_sets(dataset.id)
            self.verifyParsedOutputsEqual(iterator, "variantsets-search")

    def testSearchVariants(self):
        test_executed = 0
        start = 0
        end = 1000
        referenceName = "1"
        for dataset in self._client.search_datasets():
            for variantSet in self._client.search_variant_sets(dataset.id):
                iterator = self._client.search_variants(
                    variantSet.id, start=start, end=end,
                    reference_name=referenceName, call_set_ids=[])
                args = "--variantSetId {} --start {} --end {} -r {}".format(
                    variantSet.id, start, end, referenceName)
                test_executed += self.verifyParsedOutputsEqual(
                    iterator, "variants-search", args)
        self.assertGreater(test_executed, 0)

    def testSearchBiosamples(self):
        for dataset in self._client.search_datasets():
            iterator = self._client.search_biosamples(dataset.id)
            self.verifyParsedOutputsEqual(iterator, "biosamples-search")
        for dataset in self._client.search_datasets():
            for biosample in self._client.search_biosamples(dataset.id):
                iterator = self._client.search_biosamples(
                    dataset.id, biosample.name)
                self.verifyParsedOutputsEqual(
                    iterator,
                    "biosamples-search",
                    "--datasetId {} --name {}".format(
                        dataset.id,
                        biosample.name))
                if biosample.individualId:
                    iterator = self._client.search_biosamples(
                        dataset.id, individual_id=biosample.individualId)
                    self.verifyParsedOutputsEqual(
                        iterator,
                        "biosamples-search",
                        "--datasetId {} --individualId {}".format(
                            dataset.id,
                            biosample.individualId))

    def testSearchIndividuals(self):
        for dataset in self._client.search_datasets():
            iterator = self._client.search_individuals(dataset.id)
            self.verifyParsedOutputsEqual(iterator, "individuals-search")
        for dataset in self._client.search_datasets():
            for individual in self._client.search_individuals(dataset.id):
                iterator = self._client.search_individuals(
                    dataset.id, individual.name)
                self.verifyParsedOutputsEqual(
                    iterator, "individuals-search", "--name {}".format(
                        individual.name))

    def testSearchVariantAnnotationSets(self):
        for dataset in self._client.search_datasets():
            for variantSet in self._client.search_variant_sets(dataset.id):
                iterator = self._client.search_variant_annotation_sets(
                    variantSet.id)
                args = "{}".format(variantSet.id)
                self.verifyParsedOutputsEqual(
                    iterator, "variantannotationsets-search", args)

    def testSearchVariantAnnotations(self):
        test_executed = 0
        start = 0
        end = 10000000
        referenceName = "1"
        for dataset in self._client.search_datasets():
            for variantSet in self._client.search_variant_sets(dataset.id):
                searchIterator = self._client.search_variant_annotation_sets(
                    variantSet.id)
                for variantAnnotationSet in searchIterator:
                    iterator = self._client.search_variant_annotations(
                        variantAnnotationSet.id,
                        start=start,
                        end=end,
                        reference_name=referenceName)
                    args = ("--variantAnnotationSetId {}"
                            " --start {} --end {} -r {}").format(
                        variantAnnotationSet.id, start, end, referenceName)
                    test_executed += self.verifyParsedOutputsEqual(
                        iterator, "variantannotations-search", args)
        self.assertGreater(test_executed, 0)

    def testGetFeatures(self):
        for dataset in self._client.search_datasets():
            datasetId = dataset.id
            for featureSet in self._client.search_feature_sets(datasetId):
                for feature in self._client.search_features(
                        featureSet.id):
                    self.verifyParsedOutputsEqual(
                        [feature], "features-get", feature.id)
                    break  # this test takes too long otherwise

    def testSearchFeatures(self):
        for dataset in self._client.search_datasets():
            datasetId = dataset.id
            for featureSet in self._client.search_feature_sets(datasetId):
                iterator = self._client.search_features(featureSet.id)
                self.verifyParsedOutputsEqual(
                    iterator, "features-search",
                    "--featureSetId {}".format(featureSet.id))

    def testGetFeatureSets(self):
        for dataset in self._client.search_datasets():
            datasetId = dataset.id
            for featureSet in self._client.search_feature_sets(datasetId):
                self.verifyParsedOutputsEqual(
                    [featureSet], "featuresets-get", featureSet.id)

    def testSearchFeatureSets(self):
        for dataset in self._client.search_datasets():
            iterator = self._client.search_feature_sets(dataset.id)
            self.verifyParsedOutputsEqual(
                iterator, "featuresets-search",
                "--datasetId {}".format(dataset.id))

    def testSearchContinuous(self):
        for dataset in self._client.search_datasets():
            datasetId = dataset.id
            for continuousSet in self._client.search_continuous_sets(
                                                datasetId):
                iterator = self._client.search_continuous(
                                continuousSet.id, 'chr19', 49305897, 49306090)
                self.verifyParsedOutputsEqual(
                    iterator, "continuous-search",
                    "--continuousSetId {} --referenceName {}"
                    " --start {} --end {}".format(
                        continuousSet.id, 'chr19', 49305897, 49306090))

    def testGetContinuousSets(self):
        for dataset in self._client.search_datasets():
            datasetId = dataset.id
            for continuousSet in self._client.search_continuous_sets(
                                                datasetId):
                self.verifyParsedOutputsEqual(
                    [continuousSet], "continuoussets-get", continuousSet.id)

    def testSearchContinuousSets(self):
        for dataset in self._client.search_datasets():
            iterator = self._client.search_continuous_sets(dataset.id)
            self.verifyParsedOutputsEqual(
                iterator, "continuoussets-search",
                "--datasetId {}".format(dataset.id))

    def testSearchGenotypePhenotype(self):
        phenotype_id = "http://ohsu.edu/cgd/87795e43"
        test_executed = 0
        for dataset in self._client.search_datasets():
            # pas = phenotype_association_set
            for pas in \
              self._client.search_phenotype_association_sets(dataset.id):
                    iterator = self._client.search_genotype_phenotype(
                        phenotype_association_set_id=pas.id,
                        phenotype_ids=[phenotype_id])
                    args = (
                        "--phenotype_association_set_id {}"
                        " --phenotype_ids {} ").format(
                        pas.id, phenotype_id)
                    test_executed += self.verifyParsedOutputsEqual(
                        iterator, "genotypephenotype-search", args)
        self.assertGreater(test_executed, 0)

    def testSearchPhenotype(self):
        phenotype_id = "http://ohsu.edu/cgd/87795e43"
        test_executed = 0
        for dataset in self._client.search_datasets():
            # pas = phenotype_association_set
            for pas in \
              self._client.search_phenotype_association_sets(dataset.id):
                    iterator = self._client.search_phenotype(
                        phenotype_association_set_id=pas.id,
                        phenotype_id=phenotype_id)
                    args = (
                        "--phenotype_association_set_id {}"
                        " --phenotype_id {} ").format(
                        pas.id, phenotype_id)
                    test_executed += self.verifyParsedOutputsEqual(
                        iterator, "phenotype-search", args)
        self.assertGreater(test_executed, 0)

    def testSearchSearchPhenotypeAssociationSets(self):
        test_executed = 0
        for dataset in self._client.search_datasets():
            iterator = self._client.search_phenotype_association_sets(
                dataset_id=dataset.id)
            args = "--datasetId {}".format(dataset.id)
            test_executed += self.verifyParsedOutputsEqual(
                iterator, "phenotypeassociationsets-search", args)
        self.assertGreater(test_executed, 0)

    def testSearchExpressionLevels(self):
        for dataset in self._client.search_datasets():
            for rnaQuantificationSet in \
                    self._client.search_rna_quantification_sets(dataset.id):
                for rnaQuantification in \
                        self._client.search_rna_quantifications(
                                rnaQuantificationSet.id):
                    iterator = self._client.search_expression_levels(
                        rnaQuantification.id)
                    cliString = (
                        "expressionlevels-search "
                        "--rnaQuantificationId {}".format(
                            rnaQuantification.id))
                    self.verifyParsedOutputsEqual(iterator, cliString)

    def testSearchRnaQuantifications(self):
        for dataset in self._client.search_datasets():
            for rnaQuantificationSet in \
                    self._client.search_rna_quantification_sets(dataset.id):
                iterator = self._client.search_rna_quantifications(
                    rnaQuantificationSet.id)
                cliString = (
                    "rnaquantifications-search "
                    "--rnaQuantificationSetId {}".format(
                        rnaQuantificationSet.id))
                self.verifyParsedOutputsEqual(iterator, cliString)

    def testSearchRnaQuantificationSets(self):
        for dataset in self._client.search_datasets():
            iterator = self._client.search_rna_quantification_sets(dataset.id)
            cliString = (
                "rnaquantificationsets-search --datasetId {}".format(
                    dataset.id))
            self.verifyParsedOutputsEqual(iterator, cliString)

    def testGetExpressionLevel(self):
        for dataset in self._client.search_datasets():
            for rnaQuantificationSet in \
                    self._client.search_rna_quantification_sets(dataset.id):
                for rnaQuantification in \
                        self._client.search_rna_quantifications(
                            rnaQuantificationSet.id):
                            for expressionLevel in \
                                    self._client.search_expression_levels(
                                        rnaQuantification.id):
                                self.verifyParsedOutputsEqual(
                                    [expressionLevel],
                                    "expressionlevels-get",
                                    expressionLevel.id)

    def testGetRnaQuantification(self):
        for dataset in self._client.search_datasets():
            for rnaQuantificationSet in \
                    self._client.search_rna_quantification_sets(dataset.id):
                for rnaQuantification in \
                        self._client.search_rna_quantifications(
                            rnaQuantificationSet.id):
                    self.verifyParsedOutputsEqual(
                        [rnaQuantification],
                        "rnaquantifications-get",
                        rnaQuantification.id)

    def testGetRnaQuantificationSet(self):
        for dataset in self._client.search_datasets():
            for rnaQuantificationSet in \
                    self._client.search_rna_quantification_sets(dataset.id):
                self.verifyParsedOutputsEqual(
                    [rnaQuantificationSet],
                    "rnaquantificationsets-get",
                    rnaQuantificationSet.id)

    def testListPeers(self):
        iterator = self._client.list_peers()
        cliString = "list-peers"
        self.verifyParsedOutputsEqual(iterator, cliString)

    def testInfo(self):
        info = self._client.get_info()
        cliString = "get-info"
        self.verifyParsedOutputsEqual([info], cliString)

    def testAnnounce(self):
        url = "http://1kgenomes.ga4gh.org"
        response = self._client.announce(url)
        cliString = "announce"
        self.verifyParsedOutputsEqual([response], cliString, url)
