"""
Definitions of the GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import datetime

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


class ProtocolElementEncoder(json.JSONEncoder):
    """
    Class responsible for encoding ProtocolElements as JSON.
    """
    def default(self, obj):
        if isinstance(obj, ProtocolElement):
            ret = obj.__dict__
        else:
            ret = super(ProtocolElementEncoder, self).default(obj)
        return ret


class ProtocolElement(object):
    """
    Superclass of GA4GH protocol elements. These elements are in one-to-one
    correspondence with the Avro definitions, and provide the basic elements
    of the on-the-wire protocol.
    """
    # TODO We need to think the API for this through a little better.
    # Should the JSON converters return strings or dicts, or do we
    # need both? The equality operators also need some work and more
    # rigorous testing.

    def __str__(self):
        return "{0}({1})".format(self.__class__.__name__, self.toJSON())

    def __eq__(self, other):
        """
        Returns True if all fields in this protocol element are equal to the
        fields in the specified protocol element.

        TODO This feature is experimental and requires more testing; use
        with caution!!
        """
        # TODO This implementation is ugly! There must be a nicer way to do
        # this...
        ret = False
        if type(other) == type(self):
            ret = True
            for f in self.schema.fields:
                k = f.name
                a1 = getattr(self, k)
                a2 = getattr(other, k)
                if self.isEmbeddedType(k):
                    if isinstance(f.type, avro.schema.ArraySchema):
                        if len(a1) == len(a2):
                            ret = ret and all(x == y for x, y in zip(a1, a2))
                        else:
                            ret = False
                    else:
                        ret = ret and a1 == a2
                else:
                    ret = ret and a1 == a2
        return ret

    def __ne__(self, other):
        return not self == other

    def toJSON(self):
        """
        Returns a JSON encoded string representation of this ProtocolElement.

        TODO This method is deprecated in favour of the more explicit
        toJSONString.
        """
        return json.dumps(self, cls=ProtocolElementEncoder)

    def toJSONString(self):
        """
        Returns a JSON encoded string representation of this ProtocolElement.

        TODO is this a good API? Should we just let the client call the
        appropriate json method and use ProtocolElementEncoder?
        """
        return json.dumps(self, cls=ProtocolElementEncoder)

    def toJSONDict(self):
        """
        Returns a JSON dictionary representation of this ProtocolElement.
        """
        # TODO this is horrible! Need to streamline this and reuse some
        # of the code in _decode below
        return json.loads(json.dumps(self, cls=ProtocolElementEncoder))

    @classmethod
    def validate(cls, jsonDict):
        """
        Validates the specified json dictionary to determine if it is an
        instance of this element's schema.
        """
        return avro.io.validate(cls.schema, jsonDict)

    @classmethod
    def fromJSON(cls, jsonStr):
        """
        Returns a decoded ProtocolElement from the specified JSON string.

        TODO change the signature of this to take a json dictionary, and the
        fold in _decode in here.
        """
        d = json.loads(jsonStr)
        return cls._decode(d)

    @classmethod
    def _decode(cls, d):
        # TODO This code is needlessly obscure and needs refactoring.
        instance = cls()
        schema = cls.schema
        if d is None:
            raise ValueError("Required values not set in {0}".format(cls))
        for f in schema.fields:
            instance.__dict__[f.name] = f.default
            k = f.name
            if k in d:
                v = d[k]
                if cls.isEmbeddedType(k):
                    if isinstance(f.type, avro.schema.ArraySchema):
                        l = []
                        for element in v:
                            l.append(cls.getEmbeddedType(k)._decode(element))
                        v = l
                    elif isinstance(f.type, avro.schema.RecordSchema):
                        v = cls.getEmbeddedType(k)._decode(v)
                    elif isinstance(f.type, avro.schema.UnionSchema):
                        if v is not None:
                            v = cls.getEmbeddedType(k)._decode(v)
                    else:
                        raise Exception("Schema assumptions violated")
            instance.__dict__[k] = v
        return instance

# We can now import the definitions of the protocol elements from the
# generated file.
from _protocol_definitions import *
