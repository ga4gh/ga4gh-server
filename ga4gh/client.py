"""
Client classes for the GA4GH reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import requests
import posixpath

import ga4gh.protocol as protocol

# suppress warning about using https without cert verification
requests.packages.urllib3.disable_warnings()


class HTTPClient(object):

    workaroundGoogle = 'google'

    """
    Simple HTTP client for the GA4GH protocol.
    """
    def __init__(self, urlPrefix, debugLevel, workarounds, key):
        self._urlPrefix = urlPrefix
        self._debugLevel = debugLevel
        self._bytesRead = 0
        self._workarounds = workarounds
        self._key = key

    def runRequest(self, request, url, protocolClass, listAttr):
        """
        Runs the specified request at the specified url and instantiates
        an object of the specified class. We yield each object in listAttr.
        If pages of results are present, repeat this process until the
        pageToken is null.
        """
        notDone = True
        while notDone:
            jsonString = request.toJSONString()
            headers = {"Content-type": "application/json"}
            # make sure we correctly join with/out trailing slashes
            fullUrl = posixpath.join(self._urlPrefix, url)
            # TODO Can we get requests to output debugging information when
            # debugLevel > 0?
            authUrl = self._addAuth(fullUrl)
            response = requests.post(
                authUrl, jsonString, headers=headers, verify=False)
            response.raise_for_status()
            jsonString = response.text
            self._bytesRead += len(jsonString)
            if self._debugLevel > 1:
                # TODO use a logging output and integrate with HTTP client more
                # nicely.
                print("json response:")
                pp = json.dumps(
                    json.loads(jsonString), sort_keys=True, indent=4)
                print(pp)
            resp = protocolClass.fromJSONString(jsonString)
            # TODO handle HTTP errors from requests and display.
            for extract in getattr(resp, listAttr):
                yield extract
            request.pageToken = resp.nextPageToken
            notDone = resp.nextPageToken is not None

    def searchVariants(self, request):
        """
        Sends the specified GASearchVariantsRequest to the server and returns
        an iterator over the returned set of GAVariant objects. Result paging
        is handled transparently, so that several HTTP requests may be made
        while this method executes.
        """
        return self.runRequest(
            request, "variants/search", protocol.GASearchVariantsResponse,
            "variants")

    def searchVariantSets(self, request):
        """
        Returns an iterator over the VariantSets from the server.
        """
        return self.runRequest(
            request, "variantsets/search",
            protocol.GASearchVariantSetsResponse, "variantSets")

    def searchReferenceSets(self, request):
        """
        Returns an iterator over the ReferenceSets from the server.
        """
        return self.runRequest(
            request, "referencesets/search",
            protocol.GASearchReferenceSetsResponse, "referenceSets")

    def getBytesRead(self):
        """
        Returns the total number of (non HTTP) bytes read from the server
        by this client.
        """
        return self._bytesRead

    # TODO temporary auth solution
    def _addAuth(self, url):
        if self._usingWorkaroundsFor(self.workaroundGoogle):
            return url + "?key={0}".format(self._key)
        else:
            return url

    def _usingWorkaroundsFor(self, workaround):
        return workaround in self._workarounds
