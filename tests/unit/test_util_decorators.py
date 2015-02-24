"""
Unit tests for utility decorators
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import sys
import time

import tests.utils as utils


class TestRepeat(unittest.TestCase):
    """
    Test the Repeat decorator
    """
    # lowest possible positive number
    sleepSeconds = sys.float_info.min

    @utils.Repeat(sleepSeconds)
    def repeating(self):
        self.i += 1
        return self.i < self.endIterationAt

    def testRepeat(self):
        self.i = 0
        self.endIterationAt = 3
        self.repeating()
        self.assertEqual(self.i, self.endIterationAt)


class TestTimeout(unittest.TestCase):
    """
    Test the Timeout decorator
    """
    # lowest positive value signal.alarm allows
    timeoutSecondsLow = 1

    @utils.Timeout(timeoutSecondsLow)
    def timingOut(self):
        time.sleep(60)  # much longer than the timeout

    @utils.Timeout()
    def notTimingOut(self):
        return 0

    def testTimeoutException(self):
        with self.assertRaises(utils.TimeoutException):
            self.timingOut()

    def testTimeoutNoException(self):
        self.assertEquals(self.notTimingOut(), 0)
        # no exception thrown
