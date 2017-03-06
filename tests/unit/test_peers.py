"""
Unit tests for the peer data model.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

from ga4gh.server.datamodel import peers
from ga4gh.server import exceptions


class TestPeers(unittest.TestCase):
    """
    Unit tests for the peer datamodel.
    """
    def testBadUrlException(self):
        peers.Peer("https://1kgenomes.ga4gh.org")
        with (self.assertRaises(exceptions.BadUrlException)):
            peers.Peer("ht9://1kgenomes.ga4gh.org")
        with (self.assertRaises(exceptions.BadUrlException)):
            peers.Peer("http://1kgen")
        with (self.assertRaises(exceptions.BadUrlException)):
            peers.Peer("http://.org")

    def testToProtocolElement(self):
        url = "http://1kgenomes.ga4gh.org"
        key = "testkey"
        value = "testvalue"
        attributes = {key: value}
        peer = peers.Peer(url, attributes)
        protocolElement = peer.toProtocolElement()
        self.assertEqual(url, protocolElement.url)
        self.assertEqual(
            attributes[key],
            protocolElement.attributes.attr[key].values[0].string_value)

    def testAttributesJson(self):
        url = "http://1kgenomes.ga4gh.org"
        with (self.assertRaises(exceptions.InvalidJsonException)):
            peer = peers.Peer(url)
            peer.setAttributesJson('a{"bad": "json"}')
        with (self.assertRaises(exceptions.InvalidJsonException)):
            peer = peers.Peer(url)
            peer.setAttributesJson('{"bad"; "json"}')
        with (self.assertRaises(exceptions.InvalidJsonException)):
            peer = peers.Peer(url)
            peer.setAttributesJson('{"bad": json}')
