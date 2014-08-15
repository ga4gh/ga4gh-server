"""
Command line interface for the ga4gh reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import future
 
import argparse

import ga4gh
import ga4gh.server

class SimulateRunner(object):
    def __init__(self, args):
        backend = ga4gh.server.VariantSimulator()
        backend.randomSeed = args.seed
        backend.variantDensity = args.variant_density
        server_address = ('', args.port)
        self.http_server = ga4gh.server.HTTPServer(server_address, backend)

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
