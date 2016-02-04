"""
End to end test that invokes the repo manager
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import shutil
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
            shutil.rmtree(self.tempdir)

    def _runCmd(self, cmd, *args):
        command = [cmd] + [self.tempdir] + list(args)
        cli.repo_main(command)

    def testEndToEnd(self):
        self._runCmd("init")
        self._runCmd("add-referenceset", paths.faPath)
        self._runCmd("add-dataset", self.datasetName)
        self._runCmd("add-readgroupset", self.datasetName, paths.bamPath)
        self._runCmd("add-variantset", self.datasetName, paths.vcfDirPath)
        self._runCmd("check")
        self._runCmd("list")
        self._runCmd(
            "remove-variantset", self.datasetName, paths.variantSetName)
        self._runCmd(
            "remove-readgroupset", self.datasetName, paths.readGroupSetName)
        self._runCmd("remove-dataset", self.datasetName)
        self._runCmd("remove-referenceset", paths.referenceSetName)
        self._runCmd("destroy")
