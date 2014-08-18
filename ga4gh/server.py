"""
Server classes for the GA4GH reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import random
import datetime
import future
from future.standard_library import hooks
with hooks():
    import http.server

import ga4gh
import ga4gh.protocol as protocol


class VariantSimulator(object):
    """
    A class that simulates Variants that can be served by the GA4GH API.
    """
    def __init__(self, seed=0, numCalls=1):
        self._randomSeed = seed
        self._numCalls = numCalls
        self._variantSetId = "vs_sim"
        now = protocol.convertDatetime(datetime.datetime.now())
        self._created = now
        self._updated = now

    def generateVariant(self, referenceName, position, rng):
        """
        Generate a random variant for the specified position using the
        specified random number generator. This generator should be seeded
        with a value that is unique to this position so that the same variant
        will always be produced regardless of the order it is generated in.
        """
        v = protocol.GAVariant()
        # The id is the combination of the position, referenceName and variant
        # set id; this allows us to generate the variant from the position and
        # id.
        v.id = "{0}:{1}:{2}".format(self._variantSetId, referenceName,
                                    position)
        v.variantSetId = self._variantSetId
        v.referenceName = referenceName
        v.names = []  # What's a good model to generate these?
        v.created = self._created
        v.updated = self._updated
        v.start = position
        v.end = position + 1  # SNPs only for now
        bases = ["A", "C", "G", "T"]
        ref = rng.choice(bases)
        v.referenceBases = ref
        alt = rng.choice([b for b in bases if b != ref])
        v.alternateBases = [alt]
        v.calls = []
        for j in range(self._numCalls):
            c = protocol.GACall()
            # for now, the genotype is either [0,1], [1,1] or [1,0] with equal
            # probability; probably will want to do something more
            # sophisticated later.
            g = rng.choice([[0, 1], [1, 0], [1, 1]])
            c.genotype = g
            # TODO What is a reasonable model for generating these likelihoods?
            # Are these log-scaled? Spec does not say.
            c.genotypeLikelihood = [-100, -100, -100]
            v.calls.append(c)
        return v

    def searchVariants(self, request):
        """
        Serves the specified GASearchVariantsRequest and returns a
        GASearchVariantsResponse. If the number of variants to be returned is
        greater than request.maxResults then the nextPageToken is set to a
        non-null value. Subsequent request objects should provide this value in
        the pageToken attribute to obtain the next page of results.
        """
        response = protocol.GASearchVariantsResponse()
        rng = random.Random()
        v = []
        j = request.start
        if request.pageToken is not None:
            j = request.pageToken
        while j < request.end and len(v) != request.maxResults:
            rng.seed(self.randomSeed + j)
            if rng.random() < self.variantDensity:
                v.append(self.generateVariant(request.referenceName, j, rng))
            j += 1
        if j < request.end - 1:
            response.nextPageToken = j
        response.variants = v
        return response


class ProtocolHandler(object):
    """
    Class that handles the GA4GH protocol messages and responses.
    """
    def __init__(self, backend):
        self._backend = backend

    def searchVariants(self, jsonRequest):
        """
        Handles the specified JSON encoding of a GASearchVariantsRequest.
        and returns the corresponding JSON encoded GASearchVariantsResponse
        (in case of success) or GAException (in case of error).
        """
        request = protocol.GASearchVariantsRequest.fromJSON(jsonRequest)
        # TODO wrap this call in try: except and make a GAException object
        # out of the resulting Exception object. Two classes of exception
        # should be identified: those due to input errors and other expected
        # problems the backend must deal with, and other exceptions which
        # indicate a server error. The former type should all subclass a
        # an exception defined in the ga4gh package.
        resp = self._backend.searchVariants(request)
        s = resp.toJSON()
        return s


class HTTPRequestHandler(http.server.BaseHTTPRequestHandler):
    """
    Handler for the HTTP level of the GA4GH protocol.
    """

    def do_POST(self):
        """
        Handle a single POST request.
        """
        h = self.server.ga4ghProtocolHandler
        # TODO read the path and 404 if not correct
        length = int(self.headers['Content-Length'])
        # TODO is this safe encoding-wise? Do we need to specify an
        # explicit encoding?
        jsonRequest = self.rfile.read(length).decode()
        s = h.searchVariants(jsonRequest).encode()
        self.send_response(200)
        self.send_header("Content-type", "application/json")
        self.send_header("Content-Length", len(s))
        self.end_headers()
        self.wfile.write(s)


class HTTPServer(http.server.HTTPServer):
    """
    Basic HTTP server for the GA4GH protocol.
    """
    def __init__(self, serverAddress, backend):
        # Cannot use super() here because of Python 2 issues.
        http.server.HTTPServer.__init__(self, serverAddress,
                                        HTTPRequestHandler)
        self.ga4ghProtocolHandler = ProtocolHandler(backend)
