"""
Unit tests for test utility functions.
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


class TestCaptureOutput(unittest.TestCase):
    """
    Test that the captureOutput correctly returns the value of stdout
    and stderr for a function.
    """

    def testCapture(self):
        stdoutValue = "stdout"
        stderrValue = "stderr"

        def func():
            print(stdoutValue, file=sys.stdout, end="")
            print(stderrValue, file=sys.stderr, end="")
        stdout, stderr = utils.captureOutput(func)
        self.assertEqual(stdout, stdoutValue)
        self.assertEqual(stderr, stderrValue)

        # Empty stdout
        def func():
            print(stderrValue, file=sys.stderr, end="")
        stdout, stderr = utils.captureOutput(func)
        self.assertEqual(stdout, "")
        self.assertEqual(stderr, stderrValue)

        # Empty stderr
        def func():
            print(stdoutValue, file=sys.stdout, end="")
        stdout, stderr = utils.captureOutput(func)
        self.assertEqual(stdout, stdoutValue)
        self.assertEqual(stderr, "")

    def testArgs(self):
        def func(one, two, three, keywordOne=None, keywordTwo=None):
            print(one, two, three, keywordOne, keywordTwo, file=sys.stdout)
        stdout, stderr = utils.captureOutput(
            func, "1", "2", "3", keywordTwo="5", keywordOne="4")
        self.assertEqual(stdout, "1 2 3 4 5\n")
        self.assertEqual(stderr, "")
