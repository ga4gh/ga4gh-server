"""
Command line interface for the ga4gh reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

import ga4gh
import ga4gh.server

import werkzeug.serving

class OldServerRunner(object):
    """
    Superclass of server runner; takes care of functionality common to
    all backends.
    """
    def __init__(self, args):
        hp = ('', args.port)
        backend = self.getBackend(args)
        self._httpServer = ga4gh.server.HTTPServer(hp, backend)

    def run(self):
        self._httpServer.serve_forever()


class ServerRunner(object):
    """
    Superclass of server runner; takes care of functionality common to
    all backends.
    """
    def __init__(self, args):
        backend = self.getBackend(args)
        self._port = args.port
        self._httpHandler = ga4gh.server.HTTPHandler(backend)

    def run(self):
        werkzeug.serving.run_simple(
            '', self._port, self._httpHandler.wsgiApplication,
            use_reloader=True)



class SimulateRunner(ServerRunner):
    """
    Runner class to start the a server using the simulated backend.
    """
    def getBackend(self, args):
        backend = ga4gh.server.VariantSimulator(
            args.seed, args.numCalls, args.variantDensity)
        return backend


class WormtableRunner(ServerRunner):
    """
    Runner class to run the server using the wormtable based backend.
    """
    def getBackend(self, args):
        backend = ga4gh.server.WormtableBackend(args.dataDir)
        return backend


def main():
    parser = argparse.ArgumentParser(description="GA4GH reference server")
    # Add global options
    parser.add_argument(
        "--port", "-P", default=8000, type=int,
        help="The port to listen on")
    parser.add_argument('--verbose', '-v', action='count', default=0)
    subparsers = parser.add_subparsers(title='subcommands',)

    # help
    helpParser = subparsers.add_parser(
        "help",
        description="ga4gh_server help",
        help="show this help message and exit")
    # Simplistic simulator
    simParser = subparsers.add_parser(
        "simulate",
        description="Serve a simplistic simulated model.",
        help="Serve simulated data")
    simParser.add_argument(
        "--seed", "-s", default=0, type=int,
        help="The random seed for variants")
    simParser.add_argument(
        "--variantDensity", "-d", default=0.5, type=float,
        help="The probability that a given position is a variant")
    simParser.add_argument(
        "--numCalls", "-c", default=1, type=int,
        help="The number of GACalls returned for each variant")
    simParser.set_defaults(runner=SimulateRunner)
    # Wormtable backend
    wtbParser = subparsers.add_parser(
        "wormtable",
        description="Serve the API using a wormtable based backend.",
        help="Serve data from tables.")
    wtbParser.add_argument(
        "dataDir",
        help="The directory containing the wormtables to be served.")
    wtbParser.set_defaults(runner=WormtableRunner)

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

if __name__ == "__main__":
    main()
