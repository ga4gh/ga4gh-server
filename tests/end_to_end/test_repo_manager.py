"""
End to end test that invokes the repo manager
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import tempfile
import unittest

import ga4gh.cli as cli
import tests.paths as paths


class RepoManagerEndToEndTest(unittest.TestCase):

    datasetName = 'datasetOne'
    metadata = {'description': 'aDescription'}

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        os.rmdir(self.tempdir)

    def tearDown(self):
        if os.path.exists(self.tempdir):
            os.remove(self.tempdir)

    def _runCmd(self, cmd, *args):
        command = [cmd, self.tempdir] + list(args)
        print(command)
        cli.repo_main(command)

    def testEndToEnd(self):
        self._runCmd("init")
        self._runCmd("add-ontology", paths.ontologyPath)
        self._runCmd("add-referenceset", paths.ncbi37FaPath)
        self._runCmd("add-dataset", self.datasetName)
        self._runCmd("add-readgroupset", self.datasetName, paths.bamPath)
        self._runCmd(
            "add-featureset",
            self.datasetName,
            paths.featuresPath,
            "-R{}".format(paths.ncbi37FaName))
        # ensure we can handle trailing slashes by adding '/'
        self._runCmd(
            "add-variantset", self.datasetName,
            paths.vcfDirPath + '/', "-R{}".format(
                paths.ncbi37FaName))
        self._runCmd("verify")
        self._runCmd("list")
        self._runCmd(
            "remove-variantset", self.datasetName, paths.variantSetName,
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
            "remove-referenceset", paths.ncbi37FaName, "-f")
        self._runCmd(
            "remove-ontology", paths.ontologyName, "-f")
