"""
Exceptions for the backend
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


class BackendException(Exception):
    """
    Base class of backend exceptions
    """
    def __init__(self):
        super(BackendException, self).__init__()


class CallSetNotInVariantSetException(BackendException):
    """
    Indicates a request was made for a callSet not in the actual variantSet
    """
    exceptionMessage = "callSet '{0}' not in variantSet '{1}'"

    def __init__(self, callSetId, variantSetId):
        super(CallSetNotInVariantSetException, self).__init__()
        self.callSetId = callSetId
        self.variantSetId = variantSetId

    def toErrorMessage(self):
        return self.exceptionMessage.format(
            self.callSetId, self.variantSetId)


class RequestValidationFailureException(BackendException):
    """
    A validation of the request data failed
    """


class ResponseValidationFailureException(BackendException):
    """
    A validation of the response data failed
    """
