"""
End to end test that invokes the repo manager
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import mock
import os
import tempfile
import unittest

import ga4gh.cli.repomanager as cli_repomanager
import tests.paths as paths

import ga4gh.protocol as protocol


class RepoManagerEndToEndTest(unittest.TestCase):

    datasetName = 'datasetOne'
    metadata = {'description': 'aDescription'}
    individualName = "test"
    bioSampleName = "test"
    individual = protocol.toJson(protocol.Individual(
        name="test",
        created="2016-05-19T21:00:19Z",
        updated="2016-05-19T21:00:19Z"))
    bioSample = protocol.toJson(protocol.BioSample(
        name="test",
        created="2016-05-19T21:00:19Z",
        updated="2016-05-19T21:00:19Z"))

    def setUp(self):
        _, self.repoFile = tempfile.mkstemp(
            prefix='ga4gh_repo_manager_end2end_test')
        os.unlink(self.repoFile)

    def tearDown(self):
        if os.path.exists(self.repoFile):
            os.unlink(self.repoFile)

    def _runCmd(self, cmd, *args):
        command = [cmd, self.repoFile] + list(args)
        cli_repomanager.repo_main(command)

    def testEndToEnd(self):
        self._runCmd("init")
        self._runCmd("add-ontology", paths.ontologyPath)
        self._runCmd(
            "add-referenceset", paths.faPath,
            '-n', paths.referenceSetName)
        self._runCmd("add-dataset", self.datasetName)
        self._runCmd("add-biosample",
                     self.datasetName,
                     self.bioSampleName,
                     self.bioSample)
        self._runCmd("add-individual",
                     self.datasetName,
                     self.individualName,
                     self.individual)
        self._runCmd(
            "add-readgroupset", self.datasetName, paths.bamPath,
            '-R', paths.referenceSetName, '-n', paths.readGroupSetName)
        self._runCmd(
            "add-featureset", self.datasetName, paths.featuresPath,
            '-R', paths.referenceSetName, '-O', paths.ontologyName)
        # ensure we can handle trailing slashes
        vcfPath = paths.vcfDirPath + '/'
        self._runCmd(
            "add-variantset", self.datasetName,
            vcfPath, '-R', paths.referenceSetName)
        variantAnnotationSetName = "vas"
        self._runCmd(
            "add-variantset", self.datasetName,
            paths.annotatedVcfPath, '-R', paths.referenceSetName,
            "-aO", paths.ontologyName, "-n", variantAnnotationSetName)
        phenotypeAssociationSetName = "paSet"
        self._runCmd(
            "add-phenotypeassociationset",
            self.datasetName,
            paths.phenotypeAssociationSetPath,
            "-n",
            phenotypeAssociationSetName)

        self._runCmd("verify")
        self._runCmd("list")

        self._runCmd(
            "remove-phenotypeassociationset",
            self.datasetName, phenotypeAssociationSetName, "-f")
        self._runCmd(
            "remove-variantset", self.datasetName, paths.variantSetName, "-f")
        self._runCmd(
            "remove-variantset", self.datasetName, variantAnnotationSetName,
            "-f")
        self._runCmd(
            "remove-readgroupset", self.datasetName,
            paths.readGroupSetName, "-f")
        self._runCmd(
            "remove-featureset", self.datasetName, paths.featureSetName,
            "-f")
        self._runCmd(
            "remove-dataset", self.datasetName, "-f")
        self._runCmd(
            "remove-referenceset", paths.referenceSetName, "-f")
        self._runCmd(
            "remove-ontology", paths.ontologyName, "-f")

    def testForce(self):
        datasetName = 'dataset1'
        self._runCmd("init")
        self._runCmd("add-dataset", datasetName)
        with mock.patch('ga4gh.cli.repomanager.getRawInput', lambda x: 'N'):
            self._runCmd("remove-dataset", datasetName)
        with mock.patch('ga4gh.cli.repomanager.getRawInput', lambda x: 'y'):
            self._runCmd("remove-dataset", datasetName)
            with self.assertRaises(SystemExit):
                self._runCmd("remove-dataset", datasetName)
