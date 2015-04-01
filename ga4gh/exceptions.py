"""
Exceptions for the GA4GH server. Each exception that can occur in the server
is given a unique error code that is derived from its name.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import zlib
import inspect

import ga4gh.protocol as protocol


def getExceptionClass(errorCode):
    """
    Converts the specified error code into the corresponding class object.
    Raises a KeyError if the errorCode is not found.
    """
    classMap = {}
    for name, class_ in inspect.getmembers(sys.modules[__name__]):
        if inspect.isclass(class_) and issubclass(class_, Exception):
            classMap[class_.getErrorCode()] = class_
    return classMap[errorCode]


def getServerError(exception):
    """
    Converts the specified exception that is not a subclass of
    BaseServerException into a ServerError so that it can be correctly
    serialised and communicated back to the client.
    """
    assert not isinstance(exception, BaseServerException)
    serverException = ServerError()
    return serverException


class BaseServerException(Exception):
    """
    Superclass of all exceptions that can occur in the GA4GH reference
    server.
    """
    message = "Error code not set in exception; this is a bug."

    def __init__(self, *args, **kwargs):
        # We have this constructor so that we can always create an
        # instance of our base classes using inspect. This is useful
        # for testing.
        super(BaseServerException, self).__init__(*args, **kwargs)

    def getMessage(self):
        """
        Returns the message that we map into the GA4GH protocol message.
        For many exceptions we can simply set the message statically
        for a given class, and this message will be returned by
        default. For more elaborate messages that must be contructed
        at run time, we can override this method and use state held
        in instance variables in the exception.
        """
        return self.message

    @classmethod
    def getErrorCode(cls):
        """
        Returns an integer that can be used to identify this class.
        """
        # We use the CRC32 digest of the class name as a unique code.
        # We follow the recommendation of the Python docs to ensure
        # that this value is signed 32 bit integer.
        code = zlib.crc32(cls.__name__) & 0xffffffff
        return code


#####################################################################
#
# Exceptions that occur in the normal operation of the server
#
#####################################################################


class RuntimeException(BaseServerException):
    """
    exceptions that can occur during the processing of a client request
    """
    httpStatus = -1
    message = "Error code not set in exception; this is a bug."

    def toProtocolElement(self):
        """
        Converts this exception into the GA4GH protocol type so that
        it can be communicated back to the client.
        """
        error = protocol.GAException()
        error.errorCode = self.getErrorCode()
        error.message = self.getMessage()
        return error


class BadRequestException(RuntimeException):
    """
    A request that we don't like was sent to the server.
    """
    httpStatus = 400
    message = "Bad request"


class BadPageSizeException(BadRequestException):
    def __init__(self, pageSize):
        self.message = "Request page size '{}' is invalid".format(pageSize)


class BadPageTokenException(BadRequestException):
    message = "Request page token invalid"


class InvalidJsonException(BadRequestException):
    def __init__(self, jsonString):
        self.message = "Cannot parse JSON: '{}'".format(jsonString)


class RequestValidationFailureException(BadRequestException):
    """
    A validation of the request data failed
    """
    def __init__(self, jsonDict, requestClass):
        self.message = (
            "Request '{}' is not a valid instance of {}".format(
                jsonDict, requestClass))


class BadReadsSearchRequestBothRefs(BadRequestException):
    message = "only one of referenceId and referenceName can be specified"


class DatamodelValidationException(BadRequestException):
    """
    Some bad data was passed to us by the client that made no sense
    in the context of a particular datamodel object
    """


class NotFoundException(RuntimeException):
    """
    The superclass of all exceptions in which some resource was not
    found.
    """
    httpStatus = 404
    message = "A resource was not found"


class PathNotFoundException(NotFoundException):
    message = "The request path was not found"


class ObjectNotFoundException(NotFoundException):
    message = "The requested object was not found"


class VariantSetNotFound(NotFoundException):
    def __init__(self, variantSetId):
        self.message = "The requested VariantSet '{}' was not found".format(
            variantSetId)


class ReadGroupNotFoundException(ObjectNotFoundException):
    def __init__(self, readGroupId):
        self.message = "readGroupId '{}' not found".format(readGroupId)


class UnsupportedMediaTypeException(RuntimeException):
    httpStatus = 415
    message = "Unsupported media type"


class VersionNotSupportedException(NotFoundException):
    message = "API version not supported"


class MethodNotAllowedException(RuntimeException):
    httpStatus = 405
    message = "Method not allowed"


class NotImplementedException(RuntimeException):
    """
    Exception raised when a part of the API has not been implemented.
    """
    httpStatus = 501

    def __init__(self, message):
        self.message = message


class CallSetNotInVariantSetException(NotFoundException):
    """
    Indicates a request was made for a callSet not in the actual variantSet
    """
    def __init__(self, callSetId, variantSetId):
        self.message = "callSet '{0}' not in variantSet '{1}'".format(
            callSetId, variantSetId)


class DataException(BaseServerException):
    """
    Exceptions thrown during the server startup, and processing faulty VCFs
    """
    message = "Faulty data found or data file is missing."


class EmptyDirException(DataException):
    message = "Diretory empty, no VCF/BCF file was found"


class MalformedException(DataException):
    """
    A base exception class for exceptions thrown when faulty VCF file
    was found, which stops the server to processed variant search.
    """
    message = "Faulty entry found in data file."


class NotIndexedException(MalformedException):
    """
    VCF/BCF files must be indexed, with the same prefix.
    """
    def __init__(self, fileName):
        self.message = (
            "Indexing file is missing, File {} must be"
            " indexed.".format(fileName))


class OverlappingVcfException(MalformedException):
    """
    Exception thrown when two VCF files within a VariantSet directory
    contain records for the same contig.
    """
    def __init__(self, fileName, contig):
        self.message = (
            "VCF file '{}' contains records for contig '{}'. Other files"
            " in this VariantSet have records for this contig, and overlapping"
            " VCFs are not permitted.".format(fileName, contig))


class InconsistentMetaDataException(MalformedException):
    """
    Exception thrown when two VCF files within a VariantSet direcotry
    contain different metadata entries.
    """
    def __init__(self, fileName):
        self.fileName = fileName
        self.message = (
            "Metadata of {} is not consistent with other files"
            " in this VariantSet, and inconsistent metadata is not"
            " permitted".format(self.fileName))


class DuplicateCallSetIdException(MalformedException):
    """
    Exception thrown when the same sample ID occurs multiple times in
    one VCF file
    """
    def __init__(self, fileName, callSetId):
        self.callSetId = callSetId
        self.message = (
            "CallSetId: {} appeared multiple times in file {}. Duplicated"
            " sample ID is not permitted".format(self.callSetId, fileName))


class InconsistentCallSetIdException(MalformedException):
    """
    Exception thrown when two VCF files contain different sample IDs
    within a VariantSet directory.
    """
    def __init__(self, fileName):
        self.message = (
            "Inconsistent sample names found in {}. Sample IDs must be"
            " consistent within the same VariantSet"
            " directory.".format(fileName))


###############################################################
#
# Internal errors. These are exceptions that we regard as bugs.
#
###############################################################


class ServerError(RuntimeException):
    """
    Superclass of all exceptions that indicate a bug has occured.
    """
    httpStatus = 500
    message = "Internal Server Error"


class ResponseValidationFailureException(ServerError):
    """
    A validation of the response data failed
    """
    def __init__(self, jsonDict, requestClass):
        self.message = (
            "Response '{}' is not a valid instance of {}. "
            "Please file a bug report.".format(
                jsonDict, requestClass))
