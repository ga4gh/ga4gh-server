"""
Tests the avrotools module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.avrotools as avrotools
import ga4gh.protocol as protocol
import tests.utils as utils


class TestValidationTool(unittest.TestCase):
    """
    Tests the ValidationTool
    """
    def testNonProtocolElement(self):
        # Throws an exception when class_ is not a subclass of ProtocolElement
        with self.assertRaises(avrotools.ValidationToolException):
            avrotools.ValidationTool.getInvalidFields(object, {})

    def testLessFields(self):
        # Throws an exception when there are fields missing from the jsonDict
        for class_ in protocol.getProtocolClasses():
            with self.assertRaises(avrotools.ValidationToolException):
                avrotools.ValidationTool.getInvalidFields(class_, {})

    def testMoreFields(self):
        # Throws an exception when there are extra fields in the jsonDict
        for class_ in protocol.getProtocolClasses():
            jsonDict = class_().toJsonDict()
            jsonDict['extra'] = 'extra'
            with self.assertRaises(avrotools.ValidationToolException):
                avrotools.ValidationTool.getInvalidFields(
                    class_, jsonDict)

    def testGeneratedObjects(self):
        # Test that generated objects pass validation
        instanceGenerator = utils.InstanceGenerator()
        for class_ in protocol.getProtocolClasses():
            generatedInstance = instanceGenerator.generateInstance(class_)
            jsonDict = generatedInstance.toJsonDict()
            returnValue = avrotools.ValidationTool.getInvalidFields(
                class_, jsonDict)
            self.assertEqual(returnValue, {})
