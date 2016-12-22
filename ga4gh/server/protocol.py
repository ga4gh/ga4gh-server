"""
Definitions of the GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import json
import inspect
from sys import modules

import google.protobuf.json_format as json_format
import google.protobuf.message as message
import google.protobuf.struct_pb2 as struct_pb2

import ga4gh.schemas.pb as pb

from ga4gh.schemas._protocol_version import version  # noqa
from ga4gh.schemas.ga4gh.common_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.assay_metadata_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.metadata_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.metadata_service_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.read_service_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.reads_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.reference_service_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.references_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.variant_service_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.variants_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.allele_annotations_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.allele_annotation_service_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.sequence_annotations_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.sequence_annotation_service_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.bio_metadata_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.bio_metadata_service_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.genotype_phenotype_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.genotype_phenotype_service_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.rna_quantification_pb2 import *  # noqa
from ga4gh.schemas.ga4gh.rna_quantification_service_pb2 import *  # noqa


# A map of response objects to the name of the attribute used to
# store the values returned.
_valueListNameMap = {
    SearchVariantSetsResponse: "variant_sets",  # noqa
    SearchVariantsResponse: "variants",  # noqa
    SearchDatasetsResponse: "datasets",  # noqa
    SearchReferenceSetsResponse: "reference_sets",  # noqa
    SearchReferencesResponse: "references",  # noqa
    SearchReadGroupSetsResponse: "read_group_sets",  # noqa
    SearchReadsResponse: "alignments",  # noqa
    SearchCallSetsResponse: "call_sets",  # noqa
    SearchVariantAnnotationSetsResponse: "variant_annotation_sets",  # noqa
    SearchVariantAnnotationsResponse: "variant_annotations",  # noqa
    SearchFeatureSetsResponse: "feature_sets",  # noqa
    SearchFeaturesResponse: "features",  # noqa
    SearchBiosamplesResponse: "biosamples",  # noqa
    SearchIndividualsResponse: "individuals",  # noqa
    SearchPhenotypeAssociationSetsResponse: "phenotype_association_sets",  # noqa
    SearchPhenotypesResponse: "phenotypes",  # noqa
    SearchGenotypePhenotypeResponse: "associations",  # noqa
    SearchRnaQuantificationSetsResponse: "rna_quantification_sets",  # noqa
    SearchRnaQuantificationsResponse: "rna_quantifications",  # noqa
    SearchExpressionLevelsResponse: "expression_levels",  # noqa
}


def getValueListName(protocolResponseClass):
    """
    Returns the name of the attribute in the specified protocol class
    that is used to hold the values in a search response.
    """
    return _valueListNameMap[protocolResponseClass]


def convertDatetime(t):
    """
    Converts the specified datetime object into its appropriate protocol
    value. This is the number of milliseconds from the epoch.
    """
    epoch = datetime.datetime.utcfromtimestamp(0)
    delta = t - epoch
    millis = delta.total_seconds() * 1000
    return int(millis)


def getValueFromValue(value):
    """
    Extract the currently set field from a Value structure
    """
    if type(value) != struct_pb2.Value:
        raise TypeError("Expected a Value, but got {}".format(type(value)))
    if value.WhichOneof("kind") is None:
        raise AttributeError("Nothing set for {}".format(value))
    return getattr(value, value.WhichOneof("kind"))


def toJson(protoObject, indent=None):
    """
    Serialises a protobuf object as json
    """
    # Using the internal method because this way we can reformat the JSON
    js = json_format.MessageToDict(protoObject, True)
    return json.dumps(js, indent=indent)


def toJsonDict(protoObject):
    """
    Converts a protobuf object to the raw attributes
    i.e. a key/value dictionary
    """
    return json.loads(toJson(protoObject))


def fromJson(json, protoClass):
    """
    Deserialise json into an instance of protobuf class
    """
    return json_format.Parse(json, protoClass())


def validate(json, protoClass):
    """
    Check that json represents data that could be used to make
    a given protobuf class
    """
    try:
        fromJson(json, protoClass)
        # The json conversion automatically validates
        return True
    except Exception:
        return False


class SearchResponseBuilder(object):
    """
    A class to allow sequential building of SearchResponse objects.
    """
    def __init__(self, responseClass, pageSize, maxBufferSize):
        """
        Allocates a new SearchResponseBuilder for the specified
        responseClass, user-requested pageSize and the system mandated
        maxBufferSize (in bytes). The maxBufferSize is an
        approximate limit on the overall length of the serialised
        response.
        """
        self._responseClass = responseClass
        self._pageSize = pageSize
        self._maxBufferSize = maxBufferSize
        self._numElements = 0
        self._nextPageToken = None
        self._protoObject = responseClass()
        self._valueListName = getValueListName(responseClass)
        self._bufferSize = self._protoObject.ByteSize()

    def getPageSize(self):
        """
        Returns the page size for this SearchResponseBuilder. This is the
        user-requested maximum size for the number of elements in the
        value list.
        """
        return self._pageSize

    def getMaxBufferSize(self):
        """
        Returns the maximum internal buffer size for responses, which
        corresponds to total length (in bytes) of the serialised protobuf
        objects. This will always be less than the size of JSON output.
        """
        return self._maxBufferSize

    def getNextPageToken(self):
        """
        Returns the value of the nextPageToken for this
        SearchResponseBuilder.
        """
        return self._nextPageToken

    def setNextPageToken(self, nextPageToken):
        """
        Sets the nextPageToken to the specified value.
        """
        self._nextPageToken = nextPageToken

    def addValue(self, protocolElement):
        """
        Appends the specified protocolElement to the value list for this
        response.
        """
        self._numElements += 1
        self._bufferSize += protocolElement.ByteSize()
        attr = getattr(self._protoObject, self._valueListName)
        obj = attr.add()
        obj.CopyFrom(protocolElement)

    def isFull(self):
        """
        Returns True if the response buffer is full, and False otherwise.
        The buffer is full if either (1) the number of items in the value
        list is >= pageSize or (2) the total length of the serialised
        elements in the page is >= maxBufferSize.

        If page_size or max_response_length were not set in the request
        then they're not checked.
        """
        return (
            (self._pageSize > 0 and self._numElements >= self._pageSize) or
            (self._bufferSize >= self._maxBufferSize)
        )

    def getSerializedResponse(self):
        """
        Returns a string version of the SearchResponse that has
        been built by this SearchResponseBuilder.
        """
        self._protoObject.next_page_token = pb.string(self._nextPageToken)
        s = toJson(self._protoObject)
        return s


def getProtocolClasses(superclass=message.Message):
    """
    Returns all the protocol classes that are subclasses of the
    specified superclass. Only 'leaf' classes are returned,
    corresponding directly to the classes defined in the protocol.
    """
    # We keep a manual list of the superclasses that we define here
    # so we can filter them out when we're getting the protocol
    # classes.
    superclasses = set([message.Message])
    thisModule = modules[__name__]
    subclasses = []
    for name, class_ in inspect.getmembers(thisModule):
        if ((inspect.isclass(class_) and
                issubclass(class_, superclass) and
                class_ not in superclasses)):
            subclasses.append(class_)
    return subclasses


postMethods = \
    [('/callsets/search',
      SearchCallSetsRequest,  # noqa
      SearchCallSetsResponse),  # noqa
     ('/datasets/search',
      SearchDatasetsRequest,  # noqa
      SearchDatasetsResponse),  # noqa
     ('/readgroupsets/search',
      SearchReadGroupSetsRequest,  # noqa
      SearchReadGroupSetsResponse),  # noqa
     ('/reads/search',
      SearchReadsRequest,  # noqa
      SearchReadsResponse),  # noqa
     ('/references/search',
      SearchReferencesRequest,  # noqa
      SearchReferencesResponse),  # noqa
     ('/referencesets/search',
      SearchReferenceSetsRequest,  # noqa
      SearchReferenceSetsResponse),  # noqa
     ('/variants/search',
      SearchVariantsRequest,  # noqa
      SearchVariantsResponse),  # noqa
     ('/datasets/search',
      SearchDatasetsRequest,  # noqa
      SearchDatasetsResponse),  # noqa
     ('/callsets/search',
      SearchCallSetsRequest,  # noqa
      SearchCallSetsResponse),  # noqa
     ('/featuresets/search',
      SearchFeatureSetsRequest,  # noqa
      SearchFeatureSetsResponse),  # noqa
     ('/features/search',
      SearchFeaturesRequest,  # noqa
      SearchFeaturesResponse),  # noqa
     ('/variantsets/search',
      SearchVariantSetsRequest,  # noqa
      SearchVariantSetsResponse),  # noqa
     ('/variantannotations/search',
      SearchVariantAnnotationsRequest,  # noqa
      SearchVariantAnnotationSetsResponse),  # noqa
     ('/variantannotationsets/search',
      SearchVariantAnnotationSetsRequest,  # noqa
      SearchVariantAnnotationSetsResponse),  # noqa
     ('/rnaquantificationsets/search',
      SearchRnaQuantificationSetsRequest,  # noqa
      SearchRnaQuantificationSetsResponse),  # noqa
     ('/rnaquantifications/search',
      SearchRnaQuantificationsRequest,  # noqa
      SearchRnaQuantificationsResponse),  # noqa
     ('/expressionlevels/search',
      SearchExpressionLevelsRequest,  # noqa
      SearchExpressionLevelsResponse)]  # noqa
