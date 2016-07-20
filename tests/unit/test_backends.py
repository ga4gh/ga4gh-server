"""
Tests for the backend objects. We instantiate local copies of
the backends and invoke the entry points for the protocol methods.
We do not set up any server processes or communicate over sockets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.exceptions as exceptions
import ga4gh.backend as backend


# TODO most of the tests that were in here no longer make sense. The backend
# class is not responsible for handling parsed requests, and we need tests
# for this functionality.


class TestPrivateBackendMethods(unittest.TestCase):
    """
    keep tests of private backend methods here and not in one of the
    subclasses of TestAbstractBackend, otherwise the tests will needlessly
    be run more than once

    (they could be put in TestAbstractBackend, but I think it's a clearer
    separation to put them in their own test class)
    """

    def testParseIntegerArgument(self):
        good = {"one": "1", "minusone": "-1"}
        expected = {"one": 1, "minusone": -1}
        bad = {"string": "A", "float": "0.98"}
        self.assertEqual(backend._parseIntegerArgument({}, "missing", 0), 0)
        for key in good:
            self.assertEqual(
                backend._parseIntegerArgument(good, key, 0), expected[key])
        for key in bad:
            with self.assertRaises(exceptions.BadRequestIntegerException):
                backend._parseIntegerArgument(bad, key, 0)
