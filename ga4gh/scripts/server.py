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


class WormtableRunner(ServerRunner):
    """
    Runner class to run the server using the wormtable based backend.
    """
    def getBackend(self, args):
        backend = ga4gh.server.WormtableBackend(args.dataDir)
        return backend


class TabixRunner(ServerRunner):
    """
    Runner class to start the server using a tabix backend.
    """
    def getBackend(self, args):
        backend = ga4gh.server.TabixBackend(args.dataDir)
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
    # Wormtable backend
    wtbParser = subparsers.add_parser(
        "wormtable",
        description="Serve the API using a wormtable based backend.",
        help="Serve data from tables.")
    wtbParser.add_argument(
        "dataDir",
        help="The directory containing the wormtables to be served.")
    wtbParser.set_defaults(runner=WormtableRunner)
    # Tabix
    tabixParser = subparsers.add_parser(
        "tabix",
        description="Serve the API using a tabix based backend.",
        help="Serve data from Tabix indexed VCFs")
    tabixParser.add_argument(
        "dataDir",
        help="The directory containing VCFs")
    tabixParser.set_defaults(runner=TabixRunner)

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

if __name__ == "__main__":
    main()
