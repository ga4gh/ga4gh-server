"""
Exceptions for the GA4GH server. Each exception that can occur in the server
is given a unique error code that is derived from its name.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


class BaseClientException(Exception):
    """
    The base class for client exceptions
    """


class EmptyResponseException(BaseClientException):
    """
    The client received an empty response from the server
    """


class RequestNonSuccessException(BaseClientException):
    """
    The client received a 4xx or 5xx error code from the server
    """


class ErrantRequestException(BaseClientException):
    """
    An error was encountered in the process of creating the request
    """
