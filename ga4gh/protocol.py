"""
Definitions of the GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import json
import inspect
import datetime
import itertools
from cStringIO import StringIO

import avro.io


def convertDatetime(t):
    """
    Converts the specified datetime object into its appropriate protocol
    value. This is the number of milliseconds from the epoch.
    """
    epoch = datetime.datetime.utcfromtimestamp(0)
    delta = t - epoch
    millis = delta.total_seconds() * 1000
    return int(millis)


class SearchResponseBuilder(object):
    """
    A class to allow sequential building of SearchResponse objects.
    This is a performance tweak which allows us to substantially
    reduce the number of live objects we require in the server when
    we are building responses, as we write the JSON representation
    of ProtocolElements directly to a buffer.
    """
    def __init__(self, responseClass, pageSize, maxResponseLength):
        """
        Allocates a new SearchResponseBuilder for the specified
        subclass of SearchResponse, with the specified
        user-requested pageSize and the system mandated
        maxResponseLength (in bytes). The maxResponseLength is an
        approximate limit on the overall length of the JSON
        response.
        """
        self._responseClass = responseClass
        self._pageSize = pageSize
        self._maxResponseLength = maxResponseLength
        self._valueListBuffer = StringIO()
        self._numElements = 0
        self._nextPageToken = None

    def getPageSize(self):
        """
        Returns the page size for this SearchResponseBuilder. This is the
        user-requested maximum size for the number of elements in the
        value list.
        """
        return self._pageSize

    def getMaxResponseLength(self):
        """
        Returns the approximate maximum response length. More precisely,
        this is the total length (in bytes) of the concatenated JSON
        representations of the values in the value list after which
        we consider the buffer to be full.
        """
        return self._maxResponseLength

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
        if self._numElements > 0:
            self._valueListBuffer.write(", ")
        self._numElements += 1
        self._valueListBuffer.write(protocolElement.toJsonString())

    def isFull(self):
        """
        Returns True if the response buffer is full, and False otherwise.
        The buffer is full if either (1) the number of items in the value
        list is >= pageSize or (2) the total length of the serialised
        elements in the page is >= maxResponseLength.
        """
        return (
            self._numElements >= self._pageSize or
            self._valueListBuffer.tell() >= self._maxResponseLength)

    def getJsonString(self):
        """
        Returns a string version of the SearchResponse that has
        been built by this SearchResponseBuilder. This is a fully
        formed JSON document, and consists of the pageToken and
        the value list.
        """
        pageListString = "[{}]".format(self._valueListBuffer.getvalue())
        return '{{"nextPageToken": {},"{}": {}}}'.format(
            json.dumps(self._nextPageToken),
            self._responseClass.getValueListName(), pageListString)


class ProtocolElementEncoder(json.JSONEncoder):
    """
    Class responsible for encoding ProtocolElements as JSON.
    """
    def default(self, obj):
        if isinstance(obj, ProtocolElement):
            ret = {a: getattr(obj, a) for a in obj.__slots__}
        else:
            ret = super(ProtocolElementEncoder, self).default(obj)
        return ret


class ProtocolElement(object):
    """
    Superclass of GA4GH protocol elements. These elements are in one-to-one
    correspondence with the Avro definitions, and provide the basic elements
    of the on-the-wire protocol.
    """
    def __str__(self):
        return "{0}({1})".format(self.__class__.__name__, self.toJsonString())

    def __eq__(self, other):
        """
        Returns True if all fields in this protocol element are equal to the
        fields in the specified protocol element.
        """
        if type(other) != type(self):
            return False

        fieldNames = itertools.imap(lambda f: f.name, self.schema.fields)
        return all(getattr(self, k) == getattr(other, k) for k in fieldNames)

    def __ne__(self, other):
        return not self == other

    def toJsonString(self):
        """
        Returns a JSON encoded string representation of this ProtocolElement.
        """
        return json.dumps(self, cls=ProtocolElementEncoder)

    def toJsonDict(self):
        """
        Returns a JSON dictionary representation of this ProtocolElement.
        """
        out = {}
        for field in self.schema.fields:
            val = getattr(self, field.name)
            if self.isEmbeddedType(field.name):
                if isinstance(val, list):
                    out[field.name] = list(el.toJsonDict() for el in val)
                elif val is None:
                    out[field.name] = None
                else:
                    out[field.name] = val.toJsonDict()
            elif isinstance(val, list):
                out[field.name] = list(val)
            else:
                out[field.name] = val
        return out

    @classmethod
    def validate(cls, jsonDict):
        """
        Validates the specified JSON dictionary to determine if it is an
        instance of this element's schema.
        """
        return avro.io.validate(cls.schema, jsonDict)

    @classmethod
    def fromJsonString(cls, jsonStr):
        """
        Returns a decoded ProtocolElement from the specified JSON string.
        """
        jsonDict = json.loads(jsonStr)
        return cls.fromJsonDict(jsonDict)

    @classmethod
    def fromJsonDict(cls, jsonDict):
        """
        Returns a decoded ProtocolElement from the specified JSON dictionary.
        """
        if jsonDict is None:
            raise ValueError("Required values not set in {0}".format(cls))

        instance = cls()
        for field in cls.schema.fields:
            instanceVal = field.default
            if field.name in jsonDict:
                val = jsonDict[field.name]
                if cls.isEmbeddedType(field.name):
                    instanceVal = cls._decodeEmbedded(field, val)
                else:
                    instanceVal = val
            setattr(instance, field.name, instanceVal)
        return instance

    @classmethod
    def _decodeEmbedded(cls, field, val):
        if val is None:
            return None

        embeddedType = cls.getEmbeddedType(field.name)
        if isinstance(field.type, avro.schema.ArraySchema):
            return list(embeddedType.fromJsonDict(elem) for elem in val)
        else:
            return embeddedType.fromJsonDict(val)


class SearchRequest(ProtocolElement):
    """
    The superclass of all SearchRequest classes in the protocol.
    """


class SearchResponse(ProtocolElement):
    """
    The superclass of all SearchResponse classes in the protocol.
    """
    @classmethod
    def getValueListName(cls):
        """
        Returns the name of the list used to store the values held
        in a page of results.
        """
        return cls._valueListName


def getProtocolClasses(superclass=ProtocolElement):
    """
    Returns all the protocol classes that are subclasses of the
    specified superclass. Only 'leaf' classes are returned,
    corresponding directly to the classes defined in the protocol.
    """
    # We keep a manual list of the superclasses that we define here
    # so we can filter them out when we're getting the protocol
    # classes.
    superclasses = set([
        ProtocolElement, SearchRequest, SearchResponse])
    thisModule = sys.modules[__name__]
    subclasses = []
    for name, class_ in inspect.getmembers(thisModule):
        if ((inspect.isclass(class_) and
                issubclass(class_, superclass) and
                class_ not in superclasses)):
            subclasses.append(class_)
    return subclasses


# We can now import the definitions of the protocol elements from the
# generated file.
from _protocol_definitions import *  # NOQA
