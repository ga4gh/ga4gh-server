"""
Tests related to exceptions
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import mock
import unittest
import inspect

import ga4gh.server.exceptions as exceptions
import ga4gh.server.frontend as frontend

import ga4gh.schemas.protocol as protocol


class TestExceptionHandler(unittest.TestCase):
    """
    Test that caught exceptions are handled correctly
    """
    @classmethod
    def setUpClass(cls):
        frontend.reset()
        frontend.configure(baseConfig="TestConfig")
        frontend.app.log_exception = mock.Mock()

    class UnknownException(Exception):
        pass

    def getGa4ghException(self, data):
        return protocol.fromJson(data, protocol.GAException)

    def testObjectNotFoundException(self):
        exception = exceptions.ObjectNotFoundException()
        response = frontend.handleException(exception)
        self.assertEquals(response.status_code, 404)
        self.assertFalse(frontend.app.log_exception.called)

    def testCallSetNotInVariantSetException(self):
        exception = exceptions.CallSetNotInVariantSetException(
            'csId', 'vsId')
        response = frontend.handleException(exception)
        self.assertEquals(response.status_code, 404)
        gaException = self.getGa4ghException(response.data)
        self.assertGreater(len(gaException.message), 0)
        self.assertFalse(frontend.app.log_exception.called)

    def testUnknownExceptionBecomesServerError(self):
        exception = self.UnknownException()
        response = frontend.handleException(exception)
        self.assertEquals(response.status_code, 500)
        gaException = self.getGa4ghException(response.data)
        self.assertEquals(gaException.message, exceptions.ServerError.message)
        self.assertTrue(frontend.app.log_exception.called)

    def testNotImplementedException(self):
        message = "A string unlikely to occur at random."
        exception = exceptions.NotImplementedException(message)
        response = frontend.handleException(exception)
        self.assertEquals(response.status_code, 501)
        gaException = self.getGa4ghException(response.data)
        self.assertEquals(gaException.message, message)
        self.assertFalse(frontend.app.log_exception.called)


def isClassAndBaseServerExceptionSubclass(class_):
    return (inspect.isclass(class_) and
            issubclass(class_, exceptions.BaseServerException))


class TestExceptionConsistency(unittest.TestCase):
    """
    Ensure invariants of exceptions:
    - every exception has a non-None error code
    - every exception has a unique error code
    - every exception has an error code in the range [0, 2**31-1]
    - every exception can be instantiated successfully
    """
    def _getExceptionClasses(self):
        classes = inspect.getmembers(
            exceptions, isClassAndBaseServerExceptionSubclass)
        return [class_ for _, class_ in classes]

    def testCodeInvariants(self):
        codes = set()
        for class_ in self._getExceptionClasses():
            code = class_.getErrorCode()
            self.assertNotIn(code, codes)
            codes.add(code)
            self.assertIsNotNone(code)
            self.assertGreaterEqual(code, 0)
            self.assertLessEqual(code, 2**31 - 1)

    def testInstantiation(self):
        for class_ in self._getExceptionClasses():
            # some exceptions are becoming too complicated to instantiate
            # like the rest of the exceptions; just do them manually
            if class_ == exceptions.RequestValidationFailureException:
                objClass = protocol.SearchReadsRequest
                obj = objClass()
                obj.start = -1
                jsonDict = protocol.toJsonDict(obj)
                args = (jsonDict, objClass)
            else:
                numInitArgs = len(inspect.getargspec(
                    class_.__init__).args) - 1
                args = ['arg' for _ in range(numInitArgs)]
            instance = class_(*args)
            self.assertIsInstance(instance, exceptions.BaseServerException)
            message = instance.getMessage()
            self.assertIsInstance(message, basestring)
            self.assertGreater(len(message), 0)
            self.assertEqual(instance.getErrorCode(), class_.getErrorCode())

    def testGetExceptionClass(self):
        for class_ in self._getExceptionClasses():
            code = class_.getErrorCode()
            self.assertEqual(class_, exceptions.getExceptionClass(code))


@unittest.skip("Protobuf already does validation")
class TestValidationExceptions(unittest.TestCase):
    """
    Tests for exceptions that occur when validation fails.
    """

    def testValidationFailureExceptionMessages(self):
        # RequestValidationFailureException
        wrongString = "thisIsWrong"
        objClass = protocol.SearchReadsRequest
        obj = objClass()
        obj.start = wrongString
        jsonDict = obj.toJsonDict()
        instance = exceptions.RequestValidationFailureException(
            jsonDict, objClass)
        self.assertIn("invalid fields:", instance.message)
        self.assertIn("u'start': u'thisIsWrong'", instance.message)
        self.assertEqual(instance.message.count(wrongString), 2)
