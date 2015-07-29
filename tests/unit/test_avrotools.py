"""
Tests the avrotools module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.avrotools as avrotools
import ga4gh.protocol as protocol


class TestAvrotools(unittest.TestCase):
    """
    Tests avrotools creator and validator
    """
    def testNonProtocolElement(self):
        # Throws an exception when class_ is not a subclass of ProtocolElement
        with self.assertRaises(avrotools.AvrotoolsException):
            avrotools.Validator(object).getInvalidFields({})

    def testLessFields(self):
        # Throws an exception when there are fields missing from the jsonDict
        for class_ in protocol.getProtocolClasses():
            validator = avrotools.Validator(class_)
            with self.assertRaises(avrotools.AvrotoolsException):
                validator.getInvalidFields({})

    def testMoreFields(self):
        # Throws an exception when there are extra fields in the jsonDict
        for class_ in protocol.getProtocolClasses():
            jsonDict = class_().toJsonDict()
            jsonDict['extra'] = 'extra'
            validator = avrotools.Validator(class_)
            with self.assertRaises(avrotools.AvrotoolsException):
                validator.getInvalidFields(jsonDict)

    def testGeneratedObjects(self):
        # Test that generated objects pass validation
        for class_ in protocol.getProtocolClasses():
            creator = avrotools.Creator(class_)
            validator = avrotools.Validator(class_)
            generatedInstance = creator.getTypicalInstance()
            jsonDict = generatedInstance.toJsonDict()
            returnValue = validator.getInvalidFields(jsonDict)
            self.assertEqual(returnValue, {})
