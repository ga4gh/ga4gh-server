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

    def testOptionalFields(self):
        # test that omission of optional fields does not throw
        # an exception, but omission of a required field returns
        # a dict only mentioning that field
        for class_ in protocol.getProtocolClasses():
            validator = avrotools.Validator(class_)
            creator = avrotools.Creator(class_)
            generatedInstance = creator.getTypicalInstance()
            jsonDict = generatedInstance.toJsonDict()
            requiredFields = class_.requiredFields
            for key in jsonDict.keys():
                if key not in requiredFields:
                    del jsonDict[key]
            returnValue = validator.getInvalidFields(jsonDict)
            self.assertEqual(returnValue, {})
            if len(requiredFields) != 0:
                requiredKey = list(requiredFields)[0]
                del jsonDict[requiredKey]
                returnValue = validator.getInvalidFields(jsonDict)
                self.assertEqual(
                    returnValue[requiredKey],
                    avrotools.SchemaValidator.missingValue)
                self.assertEqual(len(returnValue), 1)

    def testLessFields(self):
        # Returns a bogus field indicator
        # when there are fields missing from the jsonDict
        for class_ in protocol.getProtocolClasses():
            validator = avrotools.Validator(class_)
            invalidFields = validator.getInvalidFields({})
            for key, value in invalidFields.items():
                self.assertEqual(
                    value, avrotools.SchemaValidator.missingValue)

    def testMoreFields(self):
        # Returns a bogus field indicator
        # when there are extra fields in the jsonDict
        key = 'extra'
        for class_ in protocol.getProtocolClasses():
            jsonDict = class_().toJsonDict()
            jsonDict[key] = None
            validator = avrotools.Validator(class_)
            invalidFields = validator.getInvalidFields(jsonDict)
            self.assertIn(key, invalidFields)
            self.assertEqual(
                invalidFields[key], avrotools.SchemaValidator.extraValue)

    def testGeneratedObjects(self):
        # Test that generated objects pass validation
        for class_ in protocol.getProtocolClasses():
            creator = avrotools.Creator(class_)
            validator = avrotools.Validator(class_)
            generatedInstance = creator.getTypicalInstance()
            jsonDict = generatedInstance.toJsonDict()
            returnValue = validator.getInvalidFields(jsonDict)
            self.assertEqual(returnValue, {})
