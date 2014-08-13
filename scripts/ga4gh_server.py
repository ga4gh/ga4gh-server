"""
Command line interface for the ga4gh reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import future
import argparse
 
from future.standard_library import hooks
with hooks():
    import http.server

import ga4gh
import ga4gh.protocol as protocol

#######################################
# TODO these classes need to be moved to within the namespace somewhere.

import random
import datetime

class VariantSimulator(object):
    """
    A class that simulates Variants that can be served by the GA4GH API.
    """
    def __init__(self, seed=0, numCalls=1):
        self.positionRng = random.Random()
        self.positionRng.seed(seed)
        self.randomSeed = seed
        self.numCalls = numCalls
        self.referenceName = "ref_sim"
        self.variantSetId = "vs_sim"
        now = protocol.convertDatetime(datetime.datetime.now()) 
        self.created = now
        self.updated = now
        

    def generateVariant(self, position):
        """
        Generate a random variant for the specified position. This is done
        by seeding the random number generator with the position and generating
        values for all other fields using this generator. The method is 
        therefore deterministic and reproducible.
        """
        rng = random.Random()
        seed = position + self.randomSeed
        rng.seed(position)
        v = protocol.GAVariant()
        # The id is the combination of the position, reference id and variant
        # set id; this allows us to generate the variant from the position and
        # id.
        v.id = "{0}:{1}:{2}".format(self.variantSetId, self.referenceName, 
                position)
        v.variantSetId = self.variantSetId 
        v.referenceName = self.referenceName
        v.created = self.created
        v.updated = self.updated
        v.start = position
        v.end = position + 1 # SNPs only for now
        bases = ["A", "C", "G", "T"]
        ref = rng.choice(bases) 
        v.referenceBases = ref
        alt = rng.choice([b for b in bases if b != ref])
        v.alternateBases = [alt]
        v.calls = []
        for j in range(self.numCalls):
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
        variants = []
        for j in range(request.start, request.end):
            # This is not quite right as depending on the range that we 
            # search, a given position will be or won't be a variant.
            # This must be fully deterministic.
            if self.positionRng.random() < self.variantDensity: 
                variants.append(self.generateVariant(j))
        response = protocol.GASearchVariantsResponse()
        response.variants = variants
        return response
       
class ProtocolHandler(object):
    """
    Class that handles the GA4GH protocol messages and responses.
    """
    def __init__(self, args):
        self.variantSimulator = VariantSimulator()
        self.variantSimulator.randomSeed = args.seed
        self.variantSimulator.variantDensity = args.variant_density

    def searchVariants(self, json_request):
        request = protocol.GASearchVariantsRequest.fromJSON(json_request)
        resp = self.variantSimulator.searchVariants(request) 
        s = resp.toJSON()
        return s
            

class GA4GHRequestHandler(http.server.BaseHTTPRequestHandler):
    
    def do_POST(self):
        # TODO read the path and 404 if not correct
        h = self.server.ga4gh_protocol_handler
        length = int(self.headers['Content-Length'])
        # TODO is this safe encoding-wise? Do we need to specify an
        # explicit encoding?
        json_request = self.rfile.read(length).decode()
        s = h.searchVariants(json_request).encode()
        self.send_response(200)
        self.send_header("Content-type", "application/json")
        self.send_header("Content-Length", len(s))
        self.end_headers()
        self.wfile.write(s)


class SimulateRunner(object):
    def __init__(self, args):
        server_address = ('', args.port)
        self.http_server = http.server.HTTPServer(server_address, 
                GA4GHRequestHandler)
        self.protocol_handler = ProtocolHandler(args)
        self.http_server.ga4gh_protocol_handler = self.protocol_handler

    def run(self):
        self.http_server.serve_forever()


def main():
    parser = argparse.ArgumentParser(description="GA4GH reference server")                                       
    # Add global options
    parser.add_argument("--port", "-P", default=8000, type=int,
            help="The port to listen on")
    parser.add_argument('--verbose', '-v', action='count')
    subparsers = parser.add_subparsers(title='subcommands',)                                             
                                                                                                                 
    # help  
    help_parser = subparsers.add_parser("help",
            description = "ga4gh_server help",
            help="show this help message and exit")
    
    sim_parser = subparsers.add_parser("simulate",
            description = "Runs the server and serves simulated data",
            help="Simulates data.")
    sim_parser.add_argument("--seed", "-s", default=0, type=int,
            help="The random seed for variants")
    sim_parser.add_argument("--variant-density", "-d", default=0.5, type=float,
            help="The probability a given position is a variant.")
    sim_parser.set_defaults(runner=SimulateRunner)

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

if __name__ == "__main__":
    main()
