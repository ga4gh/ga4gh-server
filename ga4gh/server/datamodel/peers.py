"""
Peer datamodel for exchanging data about GA4GH services.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import re
import urlparse

import ga4gh.server.exceptions as exceptions

import ga4gh.schemas.protocol as protocol


def isUrl(urlString):
    """
    Attempts to return whether a given URL string is valid by checking
    for the presence of the URL scheme and netloc using the urlparse
    module, and then using a regex.

    From http://stackoverflow.com/questions/7160737/
    """
    parsed = urlparse.urlparse(urlString)
    urlparseValid = parsed.netloc != '' and parsed.scheme != ''
    regex = re.compile(
        r'^(?:http|ftp)s?://'  # http:// or https://
        r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)'
        r'+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'  # domain...
        r'localhost|'  # localhost...
        r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'  # ...or ip
        r'(?::\d+)?'  # optional port
        r'(?:/?|[/?]\S+)$', re.IGNORECASE)

    return regex.match(urlString) and urlparseValid


class Peer(object):
    """
    This class represents an abstract Peer object.
    It sets default values and getters, as well as the
    toProtocolElement function.
    """
    def __init__(self, url, attributes={}, record=None):
        self._url = ""
        self._attributes = {}
        self.setUrl(url) \
            .setAttributes(attributes)
        if record is not None:
            self.populateFromRow(record)

    def setUrl(self, url):
        """
        Attempt to safely set the URL by string.
        """
        if isUrl(url):
            self._url = url
        else:
            raise exceptions.BadUrlException(url)
        return self

    def getUrl(self):
        return self._url

    def setAttributes(self, attributes):
        """
        Sets the attributes message to the provided value.
        """
        self._attributes = attributes
        return self

    def setAttributesJson(self, attributesJson):
        """
        Sets the attributes dictionary from a JSON string.
        """
        try:
            self._attributes = json.loads(attributesJson)
        except:
            raise exceptions.InvalidJsonException(attributesJson)
        return self

    def serializeAttributes(self, msg):
        """
        Sets the attrbutes of a message during serialization.
        """
        attributes = self.getAttributes()
        for key in attributes:
            protocol.setAttribute(
                msg.attributes.attr[key].values, attributes[key])
        return msg

    def getAttributes(self):
        """
        Returns the attributes for the DatamodelObject.
        """
        return self._attributes

    def toProtocolElement(self):
        peer = protocol.Peer()
        peer.url = self._url
        self.serializeAttributes(peer)
        return peer

    def populateFromRow(self, peerRecord):
        """
        This method accepts a model record and sets class variables.
        """
        self.setUrl(peerRecord.url) \
            .setAttributesJson(peerRecord.attributes)
        return self
