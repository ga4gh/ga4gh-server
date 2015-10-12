"""
Functions and utility classes for testing.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import StringIO
import functools
import humanize
import itertools
import os
import signal
import sys
import time


packageName = 'ga4gh'


def captureOutput(func, *args, **kwargs):
    """
    Runs the specified function and arguments, and returns the
    tuple (stdout, stderr) as strings.
    """
    stdout = sys.stdout
    sys.stdout = StringIO.StringIO()
    stderr = sys.stderr
    sys.stderr = StringIO.StringIO()
    try:
        func(*args, **kwargs)
        stdoutOutput = sys.stdout.getvalue()
        stderrOutput = sys.stderr.getvalue()
    finally:
        sys.stdout.close()
        sys.stdout = stdout
        sys.stderr.close()
        sys.stderr = stderr
    return stdoutOutput, stderrOutput


def zipLists(*lists):
    """
    Checks to see if all of the lists are the same length, and throws
    an AssertionError otherwise.  Returns the zipped lists.
    """
    length = len(lists[0])
    for list_ in lists[1:]:
        assert len(list_) == length
    return zip(*lists)


def getLinesFromLogFile(stream):
    stream.flush()
    stream.seek(0)
    lines = stream.readlines()
    return lines


def getProjectRootFilePath():
    # assumes we're in a directory one level below the project root
    return os.path.dirname(os.path.dirname(__file__))


def getGa4ghFilePath():
    return os.path.join(getProjectRootFilePath(), packageName)


def powerset(iterable, maxSets=None):
    """
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)

    See https://docs.python.org/2/library/itertools.html#recipes
    """
    s = list(iterable)
    return itertools.islice(itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(len(s) + 1)),
        0, maxSets)


# ---------------- Decorators ----------------


class TimeoutException(Exception):
    """
    A process has taken too long to execute
    """


class Timed(object):
    """
    Decorator that times a method, reporting runtime at finish
    """
    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            self.start = time.time()
            result = func(*args, **kwargs)
            self.end = time.time()
            self._report()
            return result
        return wrapper

    def _report(self):
        delta = self.end - self.start
        timeString = humanize.time.naturaldelta(delta)
        print("Finished in {} ({} seconds)".format(timeString, delta))


class Repeat(object):
    """
    A decorator to use for repeating a tagged function.
    The tagged function should return true if it wants to run again,
    and false if it wants to stop repeating.
    """
    defaultSleepSeconds = 0.1

    def __init__(self, sleepSeconds=defaultSleepSeconds):
        self.sleepSeconds = sleepSeconds

    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            while func(*args, **kwargs):
                time.sleep(self.sleepSeconds)
        return wrapper


class Timeout(object):
    """
    A decorator to use for only allowing a function to run
    for a limited amount of time
    """
    defaultTimeoutSeconds = 60

    def __init__(self, timeoutSeconds=defaultTimeoutSeconds):
        self.timeoutSeconds = timeoutSeconds

    def __call__(self, func):

        def _handle_timeout(signum, frame):
            raise TimeoutException()

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                # set the alarm and execute func
                signal.signal(signal.SIGALRM, _handle_timeout)
                signal.alarm(self.timeoutSeconds)
                result = func(*args, **kwargs)
            finally:
                # clear the alarm
                signal.alarm(0)
            return result
        return wrapper
