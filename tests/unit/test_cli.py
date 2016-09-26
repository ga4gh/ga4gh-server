"""
Tests the cli
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import mock
import unittest
import shlex

import ga4gh.cli.server as cli_server
import ga4gh.cli.client as cli_client
import ga4gh.cli.repomanager as cli_repomanager
import ga4gh.cli.ga2vcf as cli_ga2vcf
import ga4gh.cli.ga2sam as cli_ga2sam
import ga4gh.protocol as protocol
import google.protobuf.descriptor as descriptor
import google.protobuf.internal.python_message as python_message
import tests.utils as utils


class TestServerArguments(unittest.TestCase):
    """
    Tests that the server can parse expected arguments
    """
    def testParseArguments(self):
        cliInput = """--port 7777 --host 123.4.5.6 --config MockConfigName
        --config-file /path/to/config --tls --dont-use-reloader"""
        parser = cli_server.getServerParser()
        args = parser.parse_args(cliInput.split())
        self.assertEqual(args.port, 7777)
        self.assertEqual(args.host, "123.4.5.6")
        self.assertEqual(args.config, "MockConfigName")
        self.assertEqual(args.config_file, "/path/to/config")
        self.assertTrue(args.tls)
        self.assertTrue(args.dont_use_reloader)


class TestGa2VcfArguments(unittest.TestCase):
    """
    Tests the ga2vcf cli can parse all arguments it is supposed to
    """
    def testParseArguments(self):
        cliInput = """--key KEY -O vcf --outputFile /dev/null
        --referenceName REFERENCENAME --callSetIds CALL,SET,IDS --start 0
        --end 1 --pageSize 2 BASEURL VARIANTSETID"""
        parser = cli_ga2vcf.getGa2VcfParser()
        args = parser.parse_args(cliInput.split())
        self.assertEqual(args.key, "KEY")
        self.assertEqual(args.outputFormat, "vcf")
        self.assertEqual(args.outputFile, "/dev/null")
        self.assertEqual(args.referenceName, "REFERENCENAME")
        self.assertEqual(args.callSetIds, "CALL,SET,IDS")
        self.assertEqual(args.start, 0)
        self.assertEqual(args.end, 1)
        self.assertEqual(args.pageSize, 2)
        self.assertEquals(args.baseUrl, "BASEURL")
        self.assertEquals(args.variantSetId, "VARIANTSETID")


class TestGa2SamArguments(unittest.TestCase):
    """
    Tests the ga2sam cli can parse all arguments it is supposed to
    """
    def testParseArguments(self):
        cliInput = """--key KEY --outputFormat sam
        --pageSize 1 --start 2 --end 3 --outputFile OUT.SAM
        --referenceId REFERENCEID BASEURL READGROUPID"""
        parser = cli_ga2sam.getGa2SamParser()
        args = parser.parse_args(cliInput.split())
        self.assertEqual(args.key, "KEY")
        self.assertEqual(args.outputFormat, "sam")
        self.assertEqual(args.outputFile, "OUT.SAM")
        self.assertEqual(args.referenceId, "REFERENCEID")
        self.assertEqual(args.start, 2)
        self.assertEqual(args.end, 3)
        self.assertEqual(args.pageSize, 1)
        self.assertEquals(args.baseUrl, "BASEURL")
        self.assertEquals(args.readGroupId, "READGROUPID")


class TestClientArguments(unittest.TestCase):
    """
    Tests the client cli can parse all arguments it is supposed to
    and can initialize the runner in preparation for a request
    """
    class ParseFailureException(Exception):
        pass

    def _raiseParseFailureException(self, exitVal):
        raise self.ParseFailureException(exitVal)

    def setUp(self):
        self.parser = cli_client.getClientParser()

    def testParseFailure(self):
        cliInput = "invalidCommand"
        with utils.suppressOutput():
            with mock.patch('sys.exit', self._raiseParseFailureException):
                with self.assertRaises(self.ParseFailureException):
                    self.parser.parse_args(cliInput.split())

    def testOutputFormat(self):
        # Most of the commands support the output format option.
        cliInput = "variants-search BASEURL --outputFormat=json"
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.outputFormat, "json")
        cliInput = "variants-search BASEURL -O text"
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.outputFormat, "text")

    def testVariantsSearchArguments(self):
        cliInput = (
            "variants-search --referenceName REFERENCENAME "
            "--callSetIds CALL,SET,IDS --start 0 "
            "--end 1 --pageSize 2 --variantSetId VARIANTSETID BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.start, 0)
        self.assertEqual(args.end, 1)
        self.assertEqual(args.referenceName, "REFERENCENAME")
        self.assertEqual(args.callSetIds, "CALL,SET,IDS")
        self.assertEqual(args.pageSize, 2)
        self.assertEqual(args.variantSetId, "VARIANTSETID")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEqual(args.runner, cli_client.SearchVariantsRunner)

    def testVariantSetsSearchArguments(self):
        cliInput = (
            "variantsets-search --pageSize 1 --datasetId DATASETID BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 1)
        self.assertEqual(args.datasetId, "DATASETID")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchVariantSetsRunner)

    def testReferenceSetsSearchArguments(self):
        cliInput = (
            "referencesets-search --pageSize 1 --accession ACCESSION "
            "--md5checksum MD5CHECKSUM --assemblyId ASSEMBLYID "
            "BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 1)
        self.assertEqual(args.md5checksum, "MD5CHECKSUM")
        self.assertEqual(args.assemblyId, "ASSEMBLYID")
        self.assertEqual(args.accession, "ACCESSION")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchReferenceSetsRunner)

    def testReferencesSearchArguments(self):
        cliInput = (
            "references-search --pageSize 10 --accession ACCESSION "
            "--md5checksum MD5CHECKSUM BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 10)
        self.assertEqual(args.md5checksum, "MD5CHECKSUM")
        self.assertEqual(args.accession, "ACCESSION")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchReferencesRunner)

    def testReadGroupSetsSearchArguments(self):
        cliInput = (
            "readgroupsets-search --pageSize 1 --datasetId DATASETID "
            "--name NAME BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 1)
        self.assertEqual(args.datasetId, "DATASETID")
        self.assertEqual(args.name, "NAME")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchReadGroupSetsRunner)

    def testCallSetsSearchArguments(self):
        cliInput = (
            "callsets-search --pageSize 1 --name NAME "
            "--variantSetId VARIANTSETID BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 1)
        self.assertEqual(args.variantSetId, "VARIANTSETID")
        self.assertEqual(args.name, "NAME")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchCallSetsRunner)

    def testReadsSearchArguments(self):
        cliInput = (
            "reads-search --pageSize 2 --start 5 --end 10 "
            "--readGroupIds READ,GROUP,IDS --referenceId REFERENCEID "
            "BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 2)
        self.assertEqual(args.start, 5)
        self.assertEqual(args.end, 10)
        self.assertEqual(args.readGroupIds, "READ,GROUP,IDS")
        self.assertEqual(args.referenceId, "REFERENCEID")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchReadsRunner)

    def testBioSamplesSearchArguments(self):
        cliInput = (
            "biosamples-search --pageSize 2 --name BIOSAMPLENAME "
            "--datasetId DATASETID "
            "BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 2)
        self.assertEqual(args.name, "BIOSAMPLENAME")
        self.assertEqual(args.datasetId, "DATASETID")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchBioSamplesRunner)

    def testIndividualsSearchArguments(self):
        cliInput = (
            "individuals-search --pageSize 2 --name INDIVIDUALNAME "
            "--datasetId DATASETID "
            "BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 2)
        self.assertEqual(args.name, "INDIVIDUALNAME")
        self.assertEqual(args.datasetId, "DATASETID")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchIndividualsRunner)

    def testDatasetsSearchArguments(self):
        cliInput = "datasets-search BASEURL"
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchDatasetsRunner)

    def testGenotypePhenotypeSearchArguments(self):
        cliInput = (
            "genotypephenotype-search --phenotype_association_set_id SET_ID "
            "--feature_ids A,B,C --phenotype_ids D,E,F --pageSize 1 -E E1 "
            "BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 1)
        self.assertEqual(args.phenotype_association_set_id, "SET_ID")
        self.assertEqual(args.feature_ids, "A,B,C")
        self.assertEqual(args.phenotype_ids, "D,E,F")
        self.assertEqual(args.evidence, "E1")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(
            args.runner, cli_client.SearchGenotypePhenotypeRunner)

    def testPhenotypeSearchArguments(self):
        cliInput = (
            "phenotype-search --phenotype_association_set_id SET_ID "
            "--phenotype_id ID1 --description FOO --type T --age_of_onset 2 "
            "--pageSize 1 BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 1)
        self.assertEqual(args.phenotype_association_set_id, "SET_ID")
        self.assertEqual(args.phenotype_id, "ID1")
        self.assertEqual(args.description, "FOO")
        self.assertEqual(args.type, "T")
        self.assertEqual(args.age_of_onset, "2")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(args.runner, cli_client.SearchPhenotypeRunner)

    def testPhenotypeAssociationSetsSearchArguments(self):
        cliInput = (
            "phenotypeassociationsets-search --datasetId SET_ID "
            "--pageSize 1 BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 1)
        self.assertEqual(args.datasetId, "SET_ID")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(
            args.runner, cli_client.SearchPhenotypeAssociationSetsRunner)

    def verifyGetArguments(self, command, runnerClass):
        cliInput = "{} BASEURL ID".format(command)
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEqual(args.id, "ID")
        self.assertEquals(args.runner, runnerClass)

    def testReferenceSetGetArguments(self):
        self.verifyGetArguments(
            "referencesets-get", cli_client.GetReferenceSetRunner)

    def testBioSamplesGetArguments(self):
        self.verifyGetArguments(
            "biosamples-get", cli_client.GetBioSampleRunner)

    def testIndividualsGetArguments(self):
        self.verifyGetArguments(
            "individuals-get", cli_client.GetIndividualRunner)

    def testReferenceGetArguments(self):
        self.verifyGetArguments(
            "references-get", cli_client.GetReferenceRunner)

    def testReadGroupSetGetArguments(self):
        self.verifyGetArguments(
            "readgroupsets-get", cli_client.GetReadGroupSetRunner)

    def testReadGroupGetArguments(self):
        self.verifyGetArguments(
            "readgroups-get", cli_client.GetReadGroupRunner)

    def testCallSetGetArguments(self):
        self.verifyGetArguments(
            "callsets-get", cli_client.GetCallSetRunner)

    def testDatasetsGetArguments(self):
        self.verifyGetArguments(
            "datasets-get", cli_client.GetDatasetRunner)

    def testVariantGetArguments(self):
        self.verifyGetArguments(
            "variants-get", cli_client.GetVariantRunner)

    def testRnaQuantificationGetArguments(self):
        self.verifyGetArguments(
            "rnaquantifications-get", cli_client.GetRnaQuantificationRunner)

    def testVariantSetsGet(self):
        self.verifyGetArguments(
            "variantsets-get", cli_client.GetVariantSetRunner)

    def testFeaturesGet(self):
        self.verifyGetArguments(
            "features-get", cli_client.GetFeatureRunner)

    def testFeatureSetsGet(self):
        self.verifyGetArguments(
            "featuresets-get", cli_client.GetFeatureSetRunner)

    def testExpressionLevelsGet(self):
        self.verifyGetArguments(
            "expressionlevels-get", cli_client.GetExpressionLevelRunner)

    def testRnaQuantificationSetsGet(self):
        self.verifyGetArguments(
            "rnaquantificationsets-get",
            cli_client.GetRnaQuantificationSetRunner)

    def testReferenceBasesListArguments(self):
        cliInput = (
            "references-list-bases BASEURL ID --start 1 --end 2 "
            "--outputFormat fasta")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEqual(args.id, "ID")
        self.assertEqual(args.start, 1)
        self.assertEqual(args.end, 2)
        self.assertEquals(args.outputFormat, "fasta")
        self.assertEquals(args.runner, cli_client.ListReferenceBasesRunner)

    def testVariantAnnotationsSearch(self):
        cliInput = (
            "variantannotations-search "
            "--variantAnnotationSetId VARIANTANNOTATIONSETID "
            "--referenceName REFERENCENAME --start 1 "
            "--end 2 --effects EFFECTS "
            "--pageSize 3 BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(
            args.variantAnnotationSetId, "VARIANTANNOTATIONSETID")
        self.assertEqual(args.referenceName, "REFERENCENAME")
        self.assertEqual(args.start, 1)
        self.assertEqual(args.end, 2)
        self.assertEqual(args.effects, "EFFECTS")
        self.assertEqual(args.pageSize, 3)
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEquals(
            args.runner, cli_client.SearchVariantAnnotationsRunner)

    def testVariationAnnotationSetsSearch(self):
        cliInput = (
            "variantannotationsets-search "
            "--pageSize 3 BASEURL VARIANTSETID")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 3)
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEqual(args.variantSetId, "VARIANTSETID")
        self.assertEquals(
            args.runner, cli_client.SearchVariantAnnotationSetsRunner)

    def testVariationAnnotationSetsGet(self):
        cliInput = (
            "variantannotationsets-get "
            "BASEURL VARIANTANNOTATIONSETID")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEqual(args.id, "VARIANTANNOTATIONSETID")
        self.assertEquals(
            args.runner, cli_client.GetVariantAnnotationSetRunner)

    def testRnaQuantificationSearchArguments(self):
        cliInput = (
            "rnaquantifications-search --rnaQuantificationSetId ID BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.rnaQuantificationSetId, "ID")
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEqual(
            args.runner, cli_client.SearchRnaQuantificationsRunner)

    def testExpressionLevelSearchArguments(self):
        cliInput = (
            "expressionlevels-search --rnaQuantificationId rID "
            "--threshold 0.0 BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.rnaQuantificationId, "rID")
        self.assertEqual(args.threshold, 0.0)
        self.assertEqual(args.baseUrl, "BASEURL")
        self.assertEqual(args.runner, cli_client.SearchExpressionLevelsRunner)

    def testFeaturesSearch(self):
        cliInput = (
            "features-search "
            "--pageSize 3 "
            "--featureSetId FEATURESETID "
            "--start 1 "
            "--end 2 "
            "--parentId PARENTID "
            "--featureTypes FEATURE,TYPES "
            "BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 3)
        self.assertEqual(args.featureSetId, "FEATURESETID")
        self.assertEqual(args.start, 1)
        self.assertEqual(args.end, 2)
        self.assertEqual(args.parentId, "PARENTID")
        self.assertEqual(args.featureTypes, "FEATURE,TYPES")
        self.assertEqual(args.baseUrl, "BASEURL")

    def testFeatureSetsSearch(self):
        cliInput = (
            "featuresets-search "
            "--pageSize 3 "
            "--datasetId DATASETID "
            "BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 3)
        self.assertEqual(args.datasetId, "DATASETID")
        self.assertEqual(args.baseUrl, "BASEURL")

    def testRnaQuantificationSetsSearch(self):
        cliInput = (
            "rnaquantificationsets-search "
            "--pageSize 3 "
            "--datasetId DATASETID "
            "BASEURL")
        args = self.parser.parse_args(cliInput.split())
        self.assertEqual(args.pageSize, 3)
        self.assertEqual(args.datasetId, "DATASETID")
        self.assertEqual(args.baseUrl, "BASEURL")


class TestRepoManagerCli(unittest.TestCase):

    def setUp(self):
        self.parser = cli_repomanager.RepoManager.getParser()
        self.registryPath = 'a/repo/path'
        self.datasetName = "datasetName"
        self.filePath = 'a/file/path'
        self.dirPath = 'a/dir/path/'
        self.individualName = "test"
        self.bioSampleName = "test"
        self.individual = protocol.toJson(
            protocol.Individual(
                name="test",
                created="2016-05-19T21:00:19Z",
                updated="2016-05-19T21:00:19Z"))
        self.bioSample = protocol.toJson(
            protocol.BioSample(
                name="test",
                created="2016-05-19T21:00:19Z",
                updated="2016-05-19T21:00:19Z"))

    def testInit(self):
        cliInput = "init {}".format(self.registryPath)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.runner, "init")

    def testVerify(self):
        cliInput = "verify {}".format(self.registryPath)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.runner, "verify")

    def testList(self):
        cliInput = "list {}".format(self.registryPath)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.runner, "list")

    def testAddDataset(self):
        cliInput = "add-dataset {} {}".format(
            self.registryPath, self.datasetName)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.runner, "addDataset")

    def testRemoveDataset(self):
        cliInput = "remove-dataset {} {} -f".format(
            self.registryPath, self.datasetName)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.runner, "removeDataset")
        self.assertEquals(args.force, True)

    def testAddReferenceSet(self):
        description = "description"
        cliInput = (
            "add-referenceset {} {} --description={} "
            "--ncbiTaxonId NCBITAXONID "
            "--isDerived True "
            "--assemblyId ASSEMBLYID "
            "--sourceAccessions SOURCEACCESSIONS "
            "--sourceUri SOURCEURI ").format(
            self.registryPath, self.filePath, description)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.filePath, self.filePath)
        self.assertEquals(args.description, description)
        self.assertEquals(args.ncbiTaxonId, "NCBITAXONID")
        self.assertEquals(args.isDerived, True)
        self.assertEquals(args.assemblyId, "ASSEMBLYID")
        self.assertEquals(args.sourceAccessions, "SOURCEACCESSIONS")
        self.assertEquals(args.sourceUri, "SOURCEURI")
        self.assertEquals(args.runner, "addReferenceSet")

    def testRemoveReferenceSet(self):
        referenceSetName = "referenceSetName"
        cliInput = "remove-referenceset {} {} -f".format(
            self.registryPath, referenceSetName)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.referenceSetName, referenceSetName)
        self.assertEquals(args.runner, "removeReferenceSet")
        self.assertEquals(args.force, True)

    def testAddReadGroupSet(self):
        cliInput = "add-readgroupset {} {} {} ".format(
            self.registryPath, self.datasetName, self.filePath)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.dataFile, self.filePath)
        self.assertEquals(args.indexFile, None)
        self.assertEquals(args.runner, "addReadGroupSet")

    def testAddReadGroupSetWithIndexFile(self):
        indexPath = self.filePath + ".bai"
        cliInput = "add-readgroupset {} {} {} -I {}".format(
            self.registryPath, self.datasetName, self.filePath,
            indexPath)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.dataFile, self.filePath)
        self.assertEquals(args.indexFile, indexPath)
        self.assertEquals(args.runner, "addReadGroupSet")

    def testRemoveReadGroupSet(self):
        readGroupSetName = "readGroupSetName"
        cliInput = "remove-readgroupset {} {} {} -f".format(
            self.registryPath, self.datasetName, readGroupSetName)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.readGroupSetName, readGroupSetName)
        self.assertEquals(args.runner, "removeReadGroupSet")
        self.assertEquals(args.force, True)

    def testAddVariantSet(self):
        cliInput = "add-variantset {} {} {} ".format(
            self.registryPath, self.datasetName, self.filePath)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.dataFiles, [self.filePath])
        self.assertEquals(args.indexFiles, None)
        self.assertEquals(args.runner, "addVariantSet")

    def testAddVariantSetWithIndexFiles(self):
        file1 = "file1"
        file2 = "file2"
        indexFile1 = file1 + ".tbi"
        indexFile2 = file2 + ".tbi"
        cliInput = "add-variantset {} {} {} {} -I {} {}".format(
            self.registryPath, self.datasetName, file1, file2,
            indexFile1, indexFile2)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.dataFiles, [file1, file2])
        self.assertEquals(args.indexFiles, [indexFile1, indexFile2])
        self.assertEquals(args.runner, "addVariantSet")

    def testRemoveVariantSet(self):
        variantSetName = "variantSetName"
        cliInput = "remove-variantset {} {} {}".format(
            self.registryPath, self.datasetName, variantSetName)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.variantSetName, variantSetName)
        self.assertEquals(args.runner, "removeVariantSet")
        self.assertEquals(args.force, False)

    def testAddOntology(self):
        cliInput = "add-ontology {} {}".format(
            self.registryPath, self.filePath)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.filePath, self.filePath)
        self.assertEquals(args.runner, "addOntology")

    def testRemoveOntology(self):
        ontologyName = "the_ontology_name"
        cliInput = "remove-ontology {} {}".format(
            self.registryPath, ontologyName)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.ontologyName, ontologyName)
        self.assertEquals(args.runner, "removeOntology")
        self.assertEquals(args.force, False)

    def testAddBioSample(self):
        cliInput = "add-biosample {} {} {} '{}'".format(
            self.registryPath,
            self.datasetName,
            self.bioSampleName,
            self.bioSample)
        args = self.parser.parse_args(shlex.split(cliInput))
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.bioSampleName, self.bioSampleName)
        self.assertEquals(args.bioSample, self.bioSample)
        self.assertEquals(args.runner, "addBioSample")

    def testRemoveBioSample(self):
        cliInput = "remove-biosample {} {} {}".format(
            self.registryPath,
            self.datasetName,
            self.bioSampleName)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.bioSampleName, self.bioSampleName)
        self.assertEquals(args.runner, "removeBioSample")
        self.assertEquals(args.force, False)

    def testAddIndividual(self):
        cliInput = "add-individual {} {} {} '{}'".format(
            self.registryPath,
            self.datasetName,
            self.individualName,
            self.individual)
        args = self.parser.parse_args(shlex.split(cliInput))
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.individualName, self.individualName)
        self.assertEquals(args.individual, self.individual)
        self.assertEquals(args.runner, "addIndividual")

    def testRemoveIndividual(self):
        cliInput = "remove-individual {} {} {}".format(
            self.registryPath,
            self.datasetName,
            self.individualName)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.individualName, self.individualName)
        self.assertEquals(args.runner, "removeIndividual")
        self.assertEquals(args.force, False)

    def testAddPhenotypeAssociationSet(self):
        cliInput = "add-phenotypeassociationset {} {} {} -n NAME".format(
            self.registryPath,
            self.datasetName,
            self.dirPath)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.dirPath, self.dirPath)
        self.assertEquals(args.name, "NAME")

    def testRemovePhenotypeAssociationSet(self):
        cliInput = "remove-phenotypeassociationset {} {} NAME".format(
            self.registryPath, self.datasetName)
        args = self.parser.parse_args(cliInput.split())
        self.assertEquals(args.registryPath, self.registryPath)
        self.assertEquals(args.datasetName, self.datasetName)
        self.assertEquals(args.name, "NAME")


class TestOutputFormats(unittest.TestCase):
    """
    Tests the different output formats of the cli
    """
    class FakeArgs(object):
        def __init__(self, outputFormat='text'):
            self.outputFormat = outputFormat
            self.id = 'id'
            self.key = 'key'
            self.baseUrl = 'baseUrl'
            self.verbose = 'verbose'

    class FakeObject(protocol.message.Message):
        __metaclass__ = python_message.GeneratedProtocolMessageType

        FILE = descriptor.FileDescriptor(__file__, "test", "")
        DESCRIPTOR = descriptor.Descriptor(
            "FakeObject",
            "test.FakeObject",
            filename=__file__,
            file=FILE,
            containing_type=None,
            fields=[
                descriptor.FieldDescriptor(
                    name="name",
                    full_name="test.FakeObject.name",
                    index=0,
                    number=1,
                    type=descriptor.FieldDescriptor.TYPE_STRING,
                    cpp_type=descriptor.FieldDescriptor.CPPTYPE_STRING,
                    label=descriptor.FieldDescriptor.LABEL_REQUIRED,
                    default_value="",
                    message_type=None,
                    enum_type=None,
                    containing_type=None,
                    is_extension=False,
                    extension_scope=None
                ),
                descriptor.FieldDescriptor(
                    name="id",
                    full_name="test.FakeObject.id",
                    index=1,
                    number=2,
                    type=descriptor.FieldDescriptor.TYPE_STRING,
                    cpp_type=descriptor.FieldDescriptor.CPPTYPE_STRING,
                    label=descriptor.FieldDescriptor.LABEL_REQUIRED,
                    default_value="",
                    message_type=None,
                    enum_type=None,
                    containing_type=None,
                    is_extension=False,
                    extension_scope=None
                )
            ], nested_types=[], enum_types=[], extensions=[])

    def makeFakeObject(self):
        returnObj = self.FakeObject()
        returnObj.id = 'id'
        returnObj.name = 'name'
        return returnObj

    def _getRunPrintMethodCalls(self, runner):
        printCalls = []
        with mock.patch('__builtin__.print') as printMethod:
            printMethod.side_effect = \
                lambda *args, **kwargs: printCalls.append((args, kwargs))
            runner.run()
        return printCalls

    def testListReferenceBasesFasta(self):
        args = self.FakeArgs('fasta')
        args.start = 1
        args.end = 100
        returnVal = 'AGCT' * 100  # 400 bases
        runner = cli_client.ListReferenceBasesRunner(args)
        runner._client.list_reference_bases = mock.Mock(
            return_value=returnVal)
        printCalls = self._getRunPrintMethodCalls(runner)
        self.assertEqual(printCalls[0][0][0], '>id:1-100')
        self.assertEqual(len(printCalls), 7)
        self.assertEqual(
            printCalls[-1][0][0],
            returnVal[-50:])  # 50 = 400 % 70

    def testTextOutput(self):
        returnObj = self.makeFakeObject()
        args = self.FakeArgs()
        runner = cli_client.AbstractGetRunner(args)
        runner._method = mock.Mock(return_value=returnObj)
        printCalls = self._getRunPrintMethodCalls(runner)
        self.assertEqual(printCalls, [((u'id', u'name'), {'sep': u'\t'})])

    def testJsonOutput(self):
        returnObj = self.makeFakeObject()
        args = self.FakeArgs('json')
        runner = cli_client.AbstractGetRunner(args)
        runner._method = mock.Mock(return_value=returnObj)
        printCalls = self._getRunPrintMethodCalls(runner)
        self.assertEqual(json.loads(printCalls[0][0][0])['name'], 'name')
