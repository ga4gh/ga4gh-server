"""
Command line interface for the ga4gh reference implementation.
"""
from __future__ import print_function
from __future__ import division

import future
import argparse
 
from future.standard_library import hooks
with hooks():
    import http.server

import ga4gh
import ga4gh.protocol

class GA4GHRequestHandler(http.server.BaseHTTPRequestHandler):
    
    def do_POST(self):
        print("POST received", self.path)
        length = int(self.headers['Content-Length'])
        data = self.rfile.read(length)
        print("got", data)
        resp = ga4gh.protocol.GASearchVariantsResponse([])
        s = resp.to_json()

        self.send_response(200)
        self.send_header("Content-type", "application/json")
        self.send_header("Content-Length", str(s))
        self.end_headers()
        self.wfile.write(s.encode())


class SimulateRunner(object):
    def __init__(self, args):
        pass

    def run(self):
        print("running server")
        server_address = ('', 8000)
        httpd = http.server.HTTPServer(server_address, GA4GHRequestHandler)
        httpd.serve_forever()


def main():
    parser = argparse.ArgumentParser(description="GA4GH reference server")                                       
    subparsers = parser.add_subparsers(title='subcommands',)                                             
                                                                                                                 
    # help  
    help_parser = subparsers.add_parser("help",
            description = "ga4gh_server help",
            help="show this help message and exit")
    
    sim_parser = subparsers.add_parser("simulate",
            description = "Runs the server and serves simulated data",
            help="Simulates data.")
    sim_parser.set_defaults(runner=SimulateRunner)

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

if __name__ == "__main__":
    main()
