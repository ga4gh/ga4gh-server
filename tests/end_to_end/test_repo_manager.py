"""
End to end test that invokes the repo manager
"""
# TODO make this faster
# eliminate repo_dev.py invocatations by calling the cli module
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import shlex
import shutil
import subprocess
import tempfile
import unittest


class RepoManagerEndToEndTest(unittest.TestCase):

    datasetName = 'datasetOne'
    metadata = {'description': 'aDescription'}
    variantSetName = '1kgPhase1'
    readGroupSetName = 'chr17.1-250'
    referenceSetName = 'chr17'
    vcfDirPath = os.path.abspath(
        'tests/data/datasets/dataset1/variants/1kgPhase1')
    faPath = os.path.abspath(
        'tests/data/referenceSets/Default/chr17.fa.gz')
    bamPath = os.path.abspath(
        'tests/data/datasets/dataset1/reads/chr17.1-250.bam')

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        os.rmdir(self.tempdir)

    def tearDown(self):
        if os.path.exists(self.tempdir):
            shutil.rmtree(self.tempdir)

    def _runCmd(self, cmd, *args):
        argString = ' '.join(args)
        command = "python repo_dev.py {} {} {}".format(
            cmd, self.tempdir, argString)
        splits = shlex.split(command)
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(splits, stdout=devnull, stderr=devnull)

    def testEndToEnd(self):
        self._runCmd("init")
        self._runCmd("add-referenceset", self.faPath)
        self._runCmd("add-dataset", self.datasetName)
        self._runCmd("add-readgroupset", self.datasetName, self.bamPath)
        self._runCmd("add-variantset", self.datasetName, self.vcfDirPath)
        self._runCmd("check")
        self._runCmd("list")
        self._runCmd(
            "remove-variantset", self.datasetName, self.variantSetName)
        self._runCmd(
            "remove-readgroupset", self.datasetName, self.readGroupSetName)
        self._runCmd("remove-dataset", self.datasetName)
        self._runCmd("remove-referenceset", self.referenceSetName)
        self._runCmd("destroy")
