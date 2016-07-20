"""
An end to end test which tests:
- client cmd line parsing
- client operation
- client logging
- server cmd line parsing
- server operation
- simulated variantSet backend
- server logging
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import server_test
import client

import ga4gh.registry as registry


class TestGestalt(server_test.ServerTest):
    """
    An end-to-end test of the client and server
    """
    def testEndToEnd(self):
        # TODO this should be broken in to several test methods, with the
        # expensive setup done using setUpClass
        db_url = self.server.getDbUrl()
        registry_db = registry.RegistryDb(db_url)
        registry_db.open()

        dataset = registry_db.get_datasets()[0]
        datasetId = str(dataset.id)
        variantSet = dataset.variant_sets[0]
        variantSetId = str(variantSet.id)
        readGroupSet = dataset.read_group_sets[0]
        readGroupId = str(readGroupSet.read_groups[0].id)
        referenceSet = registry_db.get_reference_sets()[0]
        referenceSetId = str(referenceSet.id)
        referenceId = str(referenceSet.references[0].id)
        variantAnnotationSetId = str(variantSet.variant_annotation_sets[0].id)

        self.simulatedDatasetId = datasetId
        self.simulatedVariantSetId = variantSetId
        self.simulatedReadGroupId = readGroupId
        self.simulatedReferenceSetId = referenceSetId
        self.simulatedReferenceId = referenceId
        self.simulatedVariantAnnotationSetId = variantAnnotationSetId  # noqa
        self.client = client.ClientForTesting(self.server.getUrl())
        self.runVariantsRequest()
        self.assertLogsWritten()
        self.runReadsRequest()
        self.runReferencesRequest()
        self.runVariantSetsRequestDatasetTwo()
        self.runVariantAnnotationsRequest()
        self.runGetVariantAnnotationSetsRequest()
        self.client.cleanup()

    def assertLogsWritten(self):
        serverOutLines = self.server.getOutLines()
        serverErrLines = self.server.getErrLines()
        clientOutLines = self.client.getOutLines()
        clientErrLines = self.client.getErrLines()

        # nothing should be written to server stdout
        self.assertEqual(
            [], serverOutLines,
            "Server stdout log not empty")

        # server stderr should log at least one response success
        responseFound = False
        for line in serverErrLines:
            if ' 200 ' in line:
                responseFound = True
                break
        self.assertTrue(
            responseFound,
            "No successful server response logged to stderr")

        # client stdout should not be empty
        self.assertNotEqual(
            [], clientOutLines,
            "Client stdout log is empty")

        # number of variants to expect
        expectedNumClientOutLines = 2
        self.assertEqual(len(clientOutLines), expectedNumClientOutLines)

        # client stderr should log at least one post
        requestFound = False
        for line in clientErrLines:
            if 'POST' in line:
                requestFound = True
                break
        self.assertTrue(
            requestFound,
            "No request logged from the client to stderr")

    def runVariantsRequest(self):
        self.runClientCmd(
            self.client,
            "variants-search",
            "-s 0 -e 2 -V {}".format(self.simulatedVariantSetId))

    def runVariantAnnotationsRequest(self):
        self.runClientCmd(
            self.client,
            "variantannotations-search",
            "--variantAnnotationSetId {} -s 0 -e 2".format(
                self.simulatedVariantAnnotationSetId))

    def runGetVariantAnnotationSetsRequest(self):
        self.runClientCmd(
            self.client,
            "variantannotationsets-get",
            "{}".format(self.simulatedVariantAnnotationSetId))

    def runReadsRequest(self):
        args = "--readGroupIds {} --referenceId {} --end 10".format(
            self.simulatedReadGroupId, self.simulatedReferenceId)
        self.runClientCmd(self.client, "reads-search", args)

    def runReferencesRequest(self):
        referenceSetId = self.simulatedReferenceSetId
        referenceId = self.simulatedReferenceId
        cmd = "referencesets-search"
        self.runClientCmd(self.client, cmd)
        cmd = "references-search"
        args = "--referenceSetId={}".format(referenceSetId)
        self.runClientCmd(self.client, cmd, args)
        cmd = "referencesets-get"
        args = "{}".format(referenceSetId)
        self.runClientCmd(self.client, cmd, args)
        cmd = "references-get"
        args = "{}".format(referenceId)
        self.runClientCmd(self.client, cmd, args)
        cmd = "references-list-bases"
        args = "{}".format(referenceId)
        self.runClientCmd(self.client, cmd, args)

    def runVariantSetsRequestDatasetTwo(self):
        cmd = "variantsets-search"
        args = "--datasetId {}".format(self.simulatedDatasetId)
        self.runClientCmd(self.client, cmd, args)
