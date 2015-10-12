"""
Tests the server_benchmark command line tool
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import server_benchmark

from StringIO import StringIO


class TestServerBenchmark(unittest.TestCase):
    """
    Tests heap profile
    """
    def testHeapProfile(self):
        # capture std out
        saved_stdout, server_benchmark.sys.stdout = \
            server_benchmark.sys.stdout, StringIO('_')
        self.assertEquals(
            server_benchmark.main(['--profile', 'heap', 'tests/data']), 0)
        # store stdout in a variable
        output = server_benchmark.sys.stdout.getvalue().split('\n')
        # restore std out
        server_benchmark.sys.stdout = saved_stdout
        # we should have a string that has heap partition
        self.assertTrue("Partition of a set of" in output[0])

    """
    Tests cpu profile
    """
    def testCPUProfile(self):
        saved_stdout, server_benchmark.sys.stdout = \
            server_benchmark.sys.stdout, StringIO('_')
        self.assertEquals(
            server_benchmark.main(['--profile', 'cpu', 'tests/data']), 0)
        output = server_benchmark.sys.stdout.getvalue().split('\n')
        server_benchmark.sys.stdout = saved_stdout
        # we should have a single floating point number that represents time
        self.assertGreater(float(output[0]), 0)
        # we should have a string that has number of function calls
        self.assertTrue("function calls" in output[1])

    """
    Tests no profile
    """
    def testNoProfile(self):
        saved_stdout, server_benchmark.sys.stdout = \
            server_benchmark.sys.stdout, StringIO('_')
        self.assertEquals(
            server_benchmark.main(['--profile', 'none', 'tests/data']), 0)
        output = server_benchmark.sys.stdout.getvalue().split('\n')
        server_benchmark.sys.stdout = saved_stdout
        # we should have a single floating point number that represents time
        self.assertGreater(float(output[0]), 0)
