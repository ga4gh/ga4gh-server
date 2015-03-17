"""
The GA4GH data model. Defines all the methods required to translate
data in existing formats into GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import contextlib
import os
import sys


class DataModelObject(object):
    """
    Base class for all data model object types
    """
    def __init__(self, lowLevelOutputSuppression=True):
        self.lowLevelOutputSuppression = lowLevelOutputSuppression

    def setLowLevelOutputSuppression(self, lowLevelOutputSuppression):
        self.lowLevelOutputSuppression = lowLevelOutputSuppression

    @contextlib.contextmanager
    def suppressLowLevelOutput(self):
        # I would like to use sys.stdout.fileno() and sys.stderr.fileno()
        # here instead of literal fd numbers, but nose does something like
        # sys.stdout = StringIO.StringIO() when the -s flag is not enabled
        # (to capture test output so it doesn't get entangled with nose's
        # display) so the sys.stdout and sys.stderr objects are not able to
        # be accessed here.
        if self.lowLevelOutputSuppression:
            try:
                # redirect stdout and stderr to /dev/null
                devnull = open(os.devnull, 'w')
                origStdoutFd = 1
                origStderrFd = 2
                origStdout = os.dup(origStdoutFd)
                origStderr = os.dup(origStderrFd)
                sys.stdout.flush()
                sys.stderr.flush()
                os.dup2(devnull.fileno(), origStdoutFd)
                os.dup2(devnull.fileno(), origStderrFd)
                # enter the wrapped code
                yield
            finally:
                # restore original file streams
                devnull.flush()
                os.dup2(origStdout, origStdoutFd)
                os.dup2(origStderr, origStderrFd)
                # clean up
                os.close(origStdout)
                os.close(origStderr)
                devnull.close()
        else:
            yield
