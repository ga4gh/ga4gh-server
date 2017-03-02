"""
Provides methods to initialize the server's peer to peer connections.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import requests

import ga4gh.server.datamodel as datamodel
import ga4gh.server.exceptions as exceptions


def getInitialPeerList(filePath, logger=None):
    """
    Attempts to get a list of peers from a file specified in configuration.

    This file has one URL per line and can contain newlines and comments.

        # Main ga4gh node
        http://1kgenomes.ga4gh.org
        # Local intranet peer
        https://192.168.1.1

    The server will attempt to add URLs in this file to its registry at
    startup and will log a warning if the file isn't found.
    """
    ret = []
    with open(filePath) as textFile:
        ret = textFile.readlines()
    if len(ret) == 0:
        if logger:
            logger.warn("Couldn't load the initial "
                        "peer list. Try adding a "
                        "file named 'initial_peers.txt' "
                        "to {}".format(os.getcwd()))
    # Remove lines that start with a hash or are empty.
    return filter(lambda x: x != "" and not x.find("#") != -1, ret)


def insertInitialPeer(dataRepository, url, logger=None):
    """
    Takes the datarepository, a url, and an optional logger and attempts
    to add the peer into the repository.
    """
    insertPeer = dataRepository.insertPeer
    try:
        peer = datamodel.peers.Peer(url)
        insertPeer(peer)
    except exceptions.RepoManagerException as exc:
        if logger:
            logger.debug(
                "Peer already in registry {} {}".format(peer.getUrl(), exc))
    except exceptions.BadUrlException as exc:
        if logger:
            logger.debug("A URL in the initial "
                         "peer list {} was malformed. {}".format(url), exc)


def initialize(filePath, dataRepository, logger=None):
    for initialPeer in getInitialPeerList(filePath):
        try:
            headers = {'content-type': 'application/json'}
            data = {"peer": {"url": initialPeer}}
            url = "{}/announce".format(initialPeer.rstrip("/"))
            requests.post(url, headers=headers, json=data)
        except Exception as e:
            if logger:
                logger.info("Couldn't announce to initial peer {}".format(
                    (e, url)))
        insertInitialPeer(dataRepository, initialPeer, logger)
