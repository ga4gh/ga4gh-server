"""
Server cli
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import requests

import ga4gh.server.cli as cli
import ga4gh.server.frontend as frontend

import ga4gh.common.cli as common_cli


def addServerOptions(parser):
    parser.add_argument(
        "--port", "-P", default=8000, type=int,
        help="The port to listen on")
    parser.add_argument(
        "--host", "-H", default="127.0.0.1",
        help="The server host string; use 0.0.0.0 to allow all connections.")
    parser.add_argument(
        "--config-file", "-f", type=str, default=None,
        help="The configuration file to use")
    cli.addVersionArgument(parser)
    cli.addDisableUrllibWarningsArgument(parser)


def getServerParser():
    parser = common_cli.createArgumentParser("GA4GH reference server")
    addServerOptions(parser)
    return parser


def server_main(args=None):
    parser = getServerParser()
    parsedArgs = parser.parse_args(args)
    if parsedArgs.disable_urllib_warnings:
        requests.packages.urllib3.disable_warnings()
    frontend.configure(
        parsedArgs.config_file, parsedArgs.port)
    sslContext = None
    if frontend.app.config.get('USE_TLS') or ("OIDC_PROVIDER" in frontend.app.config):
        sslContext = "adhoc"
    frontend.app.run(
        host=parsedArgs.host, port=parsedArgs.port,
        use_reloader=frontend.app.config.get('USE_RELOADER'),
        ssl_context=sslContext)
