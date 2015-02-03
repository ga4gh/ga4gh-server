"""
Tests related to exceptions
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import inspect

import ga4gh.frontend_exceptions as frontendExceptions


class TestFrontendExceptionConsistency(unittest.TestCase):
    """
    Ensure invariants of frontend exceptions:
    - every frontend exception has a non-None error code
        - except FrontendException, which does
    - every frontend exception has a unique error code
    """

    def _getFrontendExceptionClasses(self):

        def isClassAndExceptionSubclass(clazz):
            return inspect.isclass(clazz) and issubclass(clazz, Exception)

        classes = inspect.getmembers(
            frontendExceptions, isClassAndExceptionSubclass)
        return [clazz for _, clazz in classes]

    def testCodeInvariants(self):
        codes = set()
        for clazz in self._getFrontendExceptionClasses():
            instance = clazz()
            assert instance.code not in codes
            codes.add(instance.code)
            if clazz == frontendExceptions.FrontendException:
                assert instance.code is None
            else:
                assert instance.code is not None
