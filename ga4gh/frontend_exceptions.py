"""
Error objects we return to clients
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import flask.ext.api as api

import ga4gh.backend_exceptions as backendExceptions
import ga4gh.protocol as protocol


class FrontendException(Exception):

    def __init__(self):
        super(FrontendException, self).__init__()
        self.message = ""
        self.code = None
        self.description = ""
        self.httpStatus = None

    def toGaException(self):
        error = protocol.GAException()
        error.errorCode = self.code
        error.message = self.message
        return error


class BadRequestException(FrontendException):

    def __init__(self):
        super(BadRequestException, self).__init__()
        self.httpStatus = 400
        self.message = "Bad request"
        self.code = 1


class BadPageSizeException(BadRequestException):

    def __init__(self):
        super(BadPageSizeException, self).__init__()
        self.message = "Request page size invalid"
        self.code = 2


class BadPageTokenException(BadRequestException):

    def __init__(self):
        super(BadPageTokenException, self).__init__()
        self.message = "Request page token invalid"
        self.code = 3


class NotFoundException(FrontendException):

    def __init__(self):
        super(NotFoundException, self).__init__()
        self.httpStatus = 404
        self.message = "A resource was not found"
        self.code = 4


class PathNotFoundException(NotFoundException):

    def __init__(self):
        super(PathNotFoundException, self).__init__()
        self.message = "The request path was not found"
        self.code = 5


class ObjectNotFoundException(NotFoundException):

    def __init__(self):
        super(ObjectNotFoundException, self).__init__()
        self.message = "The requested object was not found"
        self.code = 6


class ServerException(FrontendException):

    def __init__(self):
        super(ServerException, self).__init__()
        self.httpStatus = 500
        self.message = "Internal server error"
        self.code = 7


class UnsupportedMediaTypeException(FrontendException):

    def __init__(self):
        super(FrontendException, self).__init__()
        self.httpStatus = 415
        self.message = "Unsupported media type"
        self.code = 8


class VersionNotSupportedException(NotFoundException):

    def __init__(self):
        super(VersionNotSupportedException, self).__init__()
        self.message = "API version not supported"
        self.code = 9


class MethodNotAllowedException(FrontendException):

    def __init__(self):
        super(FrontendException, self).__init__()
        self.httpStatus = 405
        self.message = "Method not allowed"
        self.code = 10


class NotImplementedException(ServerException):

    def __init__(self):
        super(NotImplementedException, self).__init__()
        self.message = "Not implemented"
        self.code = 11

# exceptions thrown by the underlying system that we want to
# translate to exceptions that we define before they are
# serialized and returned to the client
exceptionMap = {
    api.exceptions.UnsupportedMediaType: UnsupportedMediaTypeException,
    backendExceptions.BackendException: ServerException,
    backendExceptions.CallSetNotInVariantSetException: NotFoundException,
    NotImplementedError: NotImplementedException,
    backendExceptions.RequestValidationFailureException: BadRequestException,
    backendExceptions.ResponseValidationFailureException: ServerException,
}
