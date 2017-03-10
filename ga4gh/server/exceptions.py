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
import json

import ga4gh.schemas.protocol as protocol


def getExceptionClass(errorCode):
    """
    Converts the specified error code into the corresponding class object.
    Raises a KeyError if the errorCode is not found.
    """
    classMap = {}
    for name, class_ in inspect.getmembers(sys.modules[__name__]):
        if inspect.isclass(class_) and issubclass(class_, BaseServerException):
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
        # that this value is signed 32 bit integer, and then mod it
        # to ensure non-negativity
        code = (zlib.crc32(cls.__name__) & 0xffffffff) % 2**31
        return code

    def __str__(self):
        return self.getMessage()


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
        error.error_code = self.getErrorCode()
        error.message = self.getMessage()
        return error


class BadRequestException(RuntimeException):
    """
    A request that we don't like was sent to the server.
    """
    httpStatus = 400
    message = "Bad request"


class BadRequestIntegerException(BadRequestException):
    def __init__(self, attrName, intString):
        self.message = \
            "{} argument '{}' could not be parsed as an integer".format(
                attrName, intString)


class BadPageSizeException(BadRequestException):
    def __init__(self, pageSize):
        self.message = "Request page size '{}' is invalid".format(pageSize)


class BadPageTokenException(BadRequestException):
    message = "Request page token invalid"


class BadIdentifierException(BadRequestException):
    def __init__(self, id_, msg=None):
        self.message = "The identifier provided is invalid: '{}' ".format(id_)
        if msg is not None:
            self.message += msg


class BadIdentifierNotStringException(BadIdentifierException):
    def __init__(self, localId):
        self.message = (
            "The identifier provided is not a string: '{}'".format(
                localId))


class InvalidJsonException(BadRequestException):
    def __init__(self, jsonString):
        self.message = "Cannot parse JSON: '{}'".format(jsonString)


class Validator(object):
    """
    Check that a JSON dictionary is a valid representation of a protocol
    buffer class
    """

    def __init__(self, class_):
        self.class_ = class_

    def getInvalidFields(self, jsonDict):
        # FIXME: get the proper list of fields
        protocol.fromJson(json.dumps(jsonDict), self.class_)
        return []


class RequestValidationFailureException(BadRequestException):
    """
    A validation of the request data failed
    """
    def __init__(self, jsonDict, requestClass):
        messageString = (
            "Request '{}' is not a valid instance of {}; "
            "invalid fields: {}")
        validator = Validator(requestClass)
        self.message = messageString.format(
            jsonDict, requestClass,
            validator.getInvalidFields(jsonDict)
        )


class BadFeatureSetSearchRequestRegularExpression(BadRequestException):
    message = "Malformed regular expression"
    httpStatus = 400


class DatamodelValidationException(BadRequestException):
    """
    Some bad data was passed to us by the client that made no sense
    in the context of a particular datamodel object
    """


class BadUrlException(RuntimeException):
    def __init__(self, url):
        self.message = "The URL: '{}' was malformed".format(url)
        self.httpStatus = 400


class ReadGroupSetNotMappedToReferenceSetException(BadRequestException):

    def __init__(self, readGroupSetId):
        self.message = (
            "ReadGroupSet '{}' is not mapped to any referenceSet".format(
                readGroupSetId))


class NotFoundException(RuntimeException):
    """
    The superclass of all exceptions in which some resource was not
    found.
    """
    httpStatus = 404
    message = "A resource was not found"


class PeerNotFoundException(NotFoundException):
    def __init__(self, url):
        self.message = "A peer with url '{}' was not found".format(
            url)


class PathNotFoundException(NotFoundException):
    message = "The request path was not found"


class ConfigurationException(BaseServerException):
    pass


class ObjectNotFoundException(NotFoundException):
    message = "The requested object was not found"


class VariantSetNotFoundException(NotFoundException):
    def __init__(self, variantSetId):
        self.message = "The requested VariantSet '{}' was not found".format(
            variantSetId)


class PhenotypeAssociationSetNotFoundException(NotFoundException):
    def __init__(self, paSetId):
        self.message = (
            "The requested PhenotypeAssociationSet '{}' "
            "was not found".format(
                paSetId))


class BiosampleNotFoundException(NotFoundException):
    def __init__(self, biosampleId):
        self.message = "The requested Biosample '{}' was not found".format(
            biosampleId)


class IndividualNotFoundException(NotFoundException):
    def __init__(self, individualId):
        self.message = "The requested Individual '{}' was not found".format(
            individualId)


class AnnotationSetNotFoundException(NotFoundException):
    def __init__(self, variantAnnotationSetId):
        self.message = "The requested VariantAnnotationSet '{}'" \
            "was not found".format(variantAnnotationSetId)


class CallSetNotFoundException(NotFoundException):
    def __init__(self, callSetId):
        self.message = "The requested CallSet '{}' was not found".format(
            callSetId)


class DatasetNotFoundException(NotFoundException):
    def __init__(self, datasetId):
        self.message = "The requested dataset '{}' was not found".format(
            datasetId)


class ReadGroupNotFoundException(ObjectNotFoundException):
    def __init__(self, readGroupId):
        self.message = "readGroupId '{}' not found".format(readGroupId)


class ReferenceSetNotFoundException(ObjectNotFoundException):
    def __init__(self, referenceSetId):
        self.message = "referenceSetId '{}' not found".format(referenceSetId)


class ReferenceNotFoundException(ObjectNotFoundException):
    def __init__(self, referenceId):
        self.message = "referenceId '{}' not found".format(referenceId)


class OntologyNotFoundException(ObjectNotFoundException):
    def __init__(self, ontologyId):
        self.message = "ontologyId '{}' not found".format(ontologyId)


class ObjectWithIdNotFoundException(ObjectNotFoundException):
    def __init__(self, objectId):
        self.message = "No object of this type exists with id '{}'".format(
            objectId)


class UnsupportedMediaTypeException(RuntimeException):
    httpStatus = 415
    message = "Unsupported media type"


class RangeErrorException(RuntimeException):
    """
    The superclass of all exceptions for which a query range error occured.
    This raises a HTTP Error 416 "Requested Range not satisfiable".
    """
    httpStatus = 416
    message = "Requested Range not satisfiable"


class ReferenceRangeErrorException(RangeErrorException):
    """
    Exception raised when the client attempts to access coordinates
    outside of the reference.
    """
    def __init__(self, referenceId, start, end):
        self.message = (
            "Query ({}, {}) outside of range for reference {}".format(
                start, end, referenceId))


class MethodNotAllowedException(RuntimeException):
    httpStatus = 405
    message = "Method not allowed"


class NotAuthenticatedException(RuntimeException):
    httpStatus = 403

    def __init__(self, message=None):
        if message is None:
            self.message = "Not authenticated. Use the " \
                           "key on the server index page."
        else:
            self.message = message


class NotAuthorizedException(RuntimeException):
    httpStatus = 401

    def __init__(self, message=None):
        if message is None:
            self.message = "Not authenticated. Use the " \
                           "key on the server index page."
        else:
            self.message = message


class NotImplementedException(RuntimeException):
    """
    Exception raised when a part of the API has not been implemented.
    """
    httpStatus = 501

    def __init__(self, message=None):
        if message is None:
            self.message = "Path not implemented"
        else:
            self.message = message


class UnmappedReadsNotSupported(NotImplementedException):
    def __init__(self):
        self.message = (
            "Unmapped reads are not yet supported; "
            "please specify a reference")


class CallSetNotInVariantSetException(NotFoundException):
    """
    Indicates a request was made for a callSet not in the actual variantSet
    """
    def __init__(self, callSetId, variantSetId):
        self.message = "callSet '{0}' not in variantSet '{1}'".format(
            callSetId, variantSetId)


class CallSetNameNotFoundException(NotFoundException):
    """
    Indicates a request was made for a callSet with a name that
    does not exist.
    """
    def __init__(self, name):
        self.message = "CallSet with name '{0}' not found".format(name)


class VariantSetNameNotFoundException(NotFoundException):
    """
    Indicates a request was made for a VariantSet with a name that
    does not exist.
    """
    def __init__(self, name):
        self.message = "VariantSet with name '{0}' not found".format(name)


class ReadGroupSetNameNotFoundException(NotFoundException):
    """
    Indicates a request was made for a ReadGroupSet with a name that
    does not exist.
    """
    def __init__(self, name):
        self.message = "ReadGroupSet with name '{0}' not found".format(name)


class ReferenceNameNotFoundException(NotFoundException):
    """
    Indicates a request was made for a Reference with a name that
    does not exist.
    """
    def __init__(self, name):
        self.message = "Reference with name '{0}' not found".format(name)


class BiosampleNameNotFoundException(NotFoundException):
    """
    Indicates a request was made for a Biosample with a name that
    does not exist.
    """
    def __init__(self, name):
        self.message = "Biosample with name '{0}' not found".format(name)


class IndividualNameNotFoundException(NotFoundException):
    """
    Indicates a request was made for an Individual with a name that
    does not exist.
    """
    def __init__(self, name):
        self.message = "Individual with name '{0}' not found".format(name)


class ReferenceSetNameNotFoundException(NotFoundException):
    """
    Indicates a request was made for a ReferenceSetSet with a name that
    does not exist.
    """
    def __init__(self, name):
        self.message = "ReferenceSet with name '{0}' not found".format(name)


class DatasetNameNotFoundException(NotFoundException):
    """
    Indicates a request was made for a Dataset with a name that
    does not exist.
    """
    def __init__(self, name):
        self.message = "Dataset with name '{0}' not found".format(name)


class OntologyNameNotFoundException(NotFoundException):
    """
    Indicates a request was made for an Ontology with a name that
    does not exist.
    """
    def __init__(self, name):
        self.message = "Ontology with name '{0}' not found".format(name)


class FeatureSetNotFoundException(NotFoundException):
    def __init__(self, featureSetId):
        self.message = (
            "FeatureSet with id '{0}' not found".format(featureSetId))


class FeatureSetNameNotFoundException(NotFoundException):
    def __init__(self, featureSetId):
        self.message = (
            "FeatureSet with name '{0}' not found".format(featureSetId))


class ParentIncompatibleWithFeatureSet(BadRequestException):
    def __init__(self):
        self.message = (
            "Parent feature incompatible with requested Feature Set."
        )


class FeatureSetNotSpecifiedException(BadRequestException):
    def __init__(self):
        self.message = (
            "One of featureSetId or parentId must be supplied."
        )


class ContinuousSetNotFoundException(NotFoundException):
    def __init__(self, continuousSetId):
        self.message = (
            "ContinuousSet with id '{0}' not found".format(continuousSetId))


class ContinuousSetNameNotFoundException(NotFoundException):
    def __init__(self, continuousSetId):
        self.message = (
            "ContinuousSet with name '{0}' not found".format(continuousSetId))


class RnaQuantificationSetNotFoundException(NotFoundException):
    def __init__(self, name):
        self.message = (
            "RnaQuantificationSet with name '{0}' not found".format(name))


class RnaQuantificationSetNameNotFoundException(NotFoundException):
    def __init__(self, name):
        self.message = (
            "RnaQuantificationSet with name '{0}' not found".format(name))


class RnaQuantificationNotFoundException(NotFoundException):
    def __init__(self, name):
        self.message = (
            "RnaQuantification with name '{0}' not found".format(name))


class ExpressionLevelNotFoundException(NotFoundException):
    def __init__(self, name):
        self.message = (
            "ExpressionLevel with name '{0}' not found".format(name))


class DataException(BaseServerException):
    """
    Exceptions thrown during data ingestion.
    """
    message = "xs"


class RepoNotFoundException(DataException):

    def __init__(self, repoPath):
        self.message = (
            "Database file '{}' does not exist.  Either "
            "change the configuration to point to a valid file or "
            "create one using the repo manager.".format(repoPath))


class RepoInvalidDatabaseException(DataException):

    def __init__(self, repoPath):
        self.message = (
            "Database file '{}' is malformed.  Either "
            "change the configuration to point to a valid file or "
            "create one using the repo manager.".format(repoPath))


class RepoSchemaVersionMismatchException(DataException):

    def __init__(self, repoVersion, expectedVersion):
        self.message = (
            "Repository version '{}' does not match expected "
            "version '{}'".format(repoVersion, expectedVersion))


class FileOpenFailedException(DataException):

    def __init__(self, filename):
        self.message = "Failed to open file '{}'".format(filename)


class EmptyDirException(DataException):

    def __init__(self, dirname, filetype):
        self.message = "Directory '{}' empty, no {} file was found".format(
            dirname, filetype)


class OntologyFileFormatException(DataException):
    """
    Exception thrown when an error occurs processing an OBO ontology
    file.
    """
    def __init__(self, filename, message):
        self.message = "Error processing ontology OBO file '{}': {}".format(
            filename, message)


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


class MultipleReferenceSetsInReadGroupSet(MalformedException):
    """
    A BAM file contains reference sequences from multiple reference
    sets.
    """
    def __init__(self, fileName, referenceSetName, otherReferenceSetName):
        self.message = (
            "The BAM file '{}' contains the referenceSets '{}' and "
            "'{}'; at most one referenceSet per file is allowed.".format(
                fileName, referenceSetName, otherReferenceSetName))


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


#####################################################################
#
# Repo manager exceptions
#
#####################################################################


class RepoManagerException(Exception):
    """
    Signals something went wrong inside the repo manager
    """


class DuplicateNameException(RepoManagerException):
    """
    The user has tried to insert an object with the same name as
    an existing object within a given container.
    """
    def __init__(self, objectName, containerName=None):
        msg = "Name '{}' already exists".format(objectName)
        if containerName is not None:
            msg += " in '{}'".format(containerName)
        super(DuplicateNameException, self).__init__(msg)


class MissingIndexException(RepoManagerException):
    """
    An index file for a remote BAM/VCF has not been provided.
    """
    def __init__(self, dataUrl):
        msg = "An index file must be provided for remote file '{}'".format(
            dataUrl)
        super(MissingIndexException, self).__init__(msg)


class UnsupportedFormatException(RepoManagerException):
    """
    The user has specified a data format which is not supported.
    """
    def __init__(self, dataFormat):
        msg = "Unsupported format: {}".format(dataFormat)
        super(UnsupportedFormatException, self).__init__(msg)
