"""
Client classes for the GA4GH reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import future
 
from future.standard_library import hooks
with hooks():
    import http.client

import ga4gh
import ga4gh.protocol as protocol

class HTTPClient(object):
    """
    Simple HTTP client for the GA4GH protocol.
    """
    def __init__(self, host, port, debugLevel):
        self.httpConnection = http.client.HTTPConnection(host, port)
        self.httpConnection.set_debuglevel(debugLevel)
        self.debugLevel = debugLevel

    def searchVariants(self, request):
        s = request.toJSON()
        self.httpConnection.request("POST", "variants/search", s)
        r = self.httpConnection.getresponse()
        s = r.read().decode()
        resp = ga4gh.protocol.GASearchVariantsResponse.fromJSON(s)
        return resp 

