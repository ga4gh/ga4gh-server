"""
Command line interface for the ga4gh reference implementation.
"""
from __future__ import print_function
from __future__ import division

import future
import argparse
 
from future.standard_library import hooks
with hooks():
    import http.client

import ga4gh
import ga4gh.protocol

class VariantSearchRunner(object):
    def __init__(self, args):
        pass

    def run(self):
        svr = ga4gh.protocol.GASearchVariantsRequest([222], 0, 1000)
        s = svr.to_json()
        print("running client")
        print("sending", s)
        c = http.client.HTTPConnection('localhost', 8000)
        c.request("POST", "variants/search", s)
        r = c.getresponse()
        print("got response", r.status, r.reason)
        print(r.read())


def main():
    parser = argparse.ArgumentParser(description="GA4GH reference client")                                       
    subparsers = parser.add_subparsers(title='subcommands',)                                             
                                                                                                                 
    # help  
    help_parser = subparsers.add_parser("help",
            description = "ga4gh_client help",
            help="show this help message and exit")
    # variants/search 
    vs_parser = subparsers.add_parser("variants-search",
            description = "Search for variants",
            help="Search for variants.")
    vs_parser.set_defaults(runner=VariantSearchRunner)

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

if __name__ == "__main__":
    main()
