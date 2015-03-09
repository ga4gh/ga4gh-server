"""
Data-driven tests for the GA4GH reference implementation. A data
driven test applies a given test method to a data file, and
each applicaton is an independent test case under nose.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import inspect


def _wrapTestMethod(method):
    """
    Takes an instance method of a DataDrivenTest subclass, and returns a
    function that can be used in a test generator in nose. This entails
    updating the description attribute so that we can see details of
    the tests being run in nose.
    """
    instance = method.im_self
    cls = instance.__class__

    def testFunction():
        method()
    testFunction.description = "{}.{}.{}:{}".format(
        method.__module__, cls.__name__, method.__name__,
        instance.getSetId())
    return testFunction


def makeTests(testDataDir, testClass):
    """
    Top-level entry point for data driven tests. For every subdirectory
    in testDataDir, create an instance of testClass and then yield
    each of its testMethods in a format suitable for use with nose
    test generators.
    """
    for testSetId in os.listdir(testDataDir):
        tester = testClass(testSetId, testDataDir)
        for name, _ in inspect.getmembers(testClass):
            if name.startswith("test"):
                yield _wrapTestMethod(getattr(tester, name))


class DataDrivenTest(object):
    """
    Superclass of all data driven tests for GA4GH datamodel objects.
    A data driven test class is instantiated with a set of data files
    that represent a some aggregation of data (for example, a
    ReferenceSet or VariantSet). We allocate a GA4GH datamodel object
    corresponding to this, and then test that these objects have the
    properties that we expect.
    """
    def __init__(self, setId, baseDir):
        self._setId = setId
        self._dataDir = os.path.join(baseDir, setId)
        self._gaObject = self.getDataModelClass()(
            self._setId, self._dataDir)

    def getSetId(self):
        """
        Return the ID of this GA4GH datamodel object we are testing.
        """
        return self._setId

    def getDataModelClass(self):
        """
        Returns the GA4GH datamodel class that this data driven test
        is exercising.
        """
        raise NotImplementedError()

    def getProtocolClass(self):
        """
        Returns the GA4GH protocol class that this data driven test
        is exercising.
        """
        raise NotImplementedError()

    def testProtocolElementValid(self):
        protocolElement = self._gaObject.toProtocolElement()
        jsonDict = protocolElement.toJsonDict()
        assert self.getProtocolClass().validate(jsonDict)
