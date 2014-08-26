"""
Definitions of the GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import datetime


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
    _embeddedTypes = {}

    def toJSON(self):
        """
        Returns a JSON encoded string representation of this ProtocolElement.
        """
        return json.dumps(self, cls=ProtocolElementEncoder)

    @classmethod
    def fromJSON(cls, json_str):
        """
        Returns a decoded ProtocolElement from the specified JSON string.
        """
        d = json.loads(json_str)
        return cls._decode(d)

    @classmethod
    def _decode(cls, d):
        instance = cls()
        for k, v in d.items():
            if k in cls._embeddedTypes:
                # TODO is this always a list?
                l = []
                for element in v:
                    l.append(cls._embeddedTypes[k]._decode(element))
                v = l
            instance.__dict__[k] = v
        return instance


class GAKeyValue(ProtocolElement):
    """
    A structure for encoding arbitrary Key-Value tuples, or tags, on other
    record types.

    TODO for convenience we pass the values in at the constructor, which is
    inconsistent with the other classes. We may want to remove this.
    """
    def __init__(self, key=None, value=None):
        self.key = key
        self.value = value


class GACall(ProtocolElement):
    """
    A GACall represents the determination of genotype with respect to a
    particular variant. It may include associated information such as quality
    and phasing. For example, a call might assign a probability of 0.32 to the
    occurrence of a SNP named rs1234 in a call set with the name NA12345
    """
    def __init__(self):
        self.callSetId = None
        self.callSetName = None
        self.genotype = []
        self.phaseset = None
        self.genotypeLikelihood = []
        self.info = []


class GAVariant(ProtocolElement):
    """
    A GAVariant represents a change in DNA sequence relative to some
    reference. For example, a variant could represent a SNP or an
    insertion. Variants belong to a GAVariantSet. This is equivalent to a
    row in VCF.
    """
    _embeddedTypes = {
        "calls": GACall,
        "info": GAKeyValue
    }

    def __init__(self):
        self.id = ""
        self.variantSetId = ""
        self.names = []
        self.created = None
        self.updated = None
        self.referenceName = ""
        self.start = None
        self.end = None
        self.referenceBases = ""
        self.alternateBases = []
        self.info = []
        self.calls = []


class GASearchVariantsRequest(ProtocolElement):
    """
    Search for variants.
    """
    def __init__(self):
        self.variantSetIds = []
        self.variantName = None
        self.callSetIds = []
        self.referenceName = None
        self.start = None
        self.end = None
        self.pageToken = None
        self.maxResults = 10  # Isn't this a bit small?


class GASearchVariantsResponse(ProtocolElement):
    _embeddedTypes = {"variants": GAVariant}

    def __init__(self):
        self.variants = []
        self.nextPageToken = None
