"""
Tests the server_benchmark command line tool
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import server_benchmark
import tests.utils as utils


class TestServerBenchmark(unittest.TestCase):
    """
    Tests heap profile
    """
    def testHeapProfile(self):
        # capture std out
        stdout, stderr = utils.captureOutput(
            server_benchmark.main, ['--profile', 'heap', 'tests/data'])
        # store stdout lines in array
        output = stdout.split('\n')
        # we should have a string that has heap partition
        self.assertTrue("Partition of a set of" in output[0])

    """
    Tests cpu profile
    """
    def testCPUProfile(self):
        # capture std out
        stdout, stderr = utils.captureOutput(
            server_benchmark.main, ['--profile', 'cpu', 'tests/data'])
        # store stdout lines in array
        output = stdout.split('\n')
        # we should have a single floating point number that represents time
        self.assertGreater(float(output[0]), 0)
        # we should have a string that has number of function calls
        self.assertTrue("function calls" in output[1])

    """
    Tests no profile
    """
    def testNoProfile(self):
        # capture std out
        stdout, stderr = utils.captureOutput(
            server_benchmark.main, ['--profile', 'none', 'tests/data'])
        # store stdout lines in array
        output = stdout.split('\n')
        # we should have a single floating point number that represents time
        self.assertGreater(float(output[0]), 0)
