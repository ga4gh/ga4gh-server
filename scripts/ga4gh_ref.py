"""
Command line interface for the ga4gh reference implementation.
"""
from __future__ import print_function
from __future__ import division

import argparse
import http.server
import http.client

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


class ServerRunner(object):
    def __init__(self, args):
        pass

    def run(self):
        print("running server")
        server_address = ('', 8000)
        httpd = http.server.HTTPServer(server_address, GA4GHRequestHandler)
        httpd.serve_forever()

class ClientRunner(object):
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
    parser = argparse.ArgumentParser(description="Reference client/server")                                       
    subparsers = parser.add_subparsers(title='subcommands',)                                             
                                                                                                                 
    # help  
    help_parser = subparsers.add_parser("help",
            description = "ga4gh_ref help",
            help="show this help message and exit")
    
    server_parser = subparsers.add_parser("server",
            description = "Runs the server process",
            help="Runs the server process")
    server_parser.set_defaults(runner=ServerRunner)

    client_parser = subparsers.add_parser("client",
            description = "Runs the client process",
            help="Runs the client process")
    client_parser.set_defaults(runner=ClientRunner)

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

if __name__ == "__main__":
    main()
