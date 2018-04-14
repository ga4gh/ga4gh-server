"""
Server cli
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import requests
import multiprocessing
import gunicorn.app.base

import ga4gh.server.cli as cli
import ga4gh.server.frontend as frontend

import ga4gh.common.cli as common_cli


class StandaloneApplication(gunicorn.app.base.BaseApplication):
    def __init__(self, app, options=None):
        self.options = options or {}
        self.application = app
        super(StandaloneApplication, self).__init__()

    def load_config(self):
        config = dict(
            [(key, value) for key, value in self.options.iteritems()
             if key in self.cfg.settings and value is not None])
        for key, value in config.iteritems():
            self.cfg.set(key.lower(), value)

    def load(self):
        return self.application


def addServerOptions(parser):
    parser.add_argument(
        "--port", "-P", default=8000, type=int,
        help="The port to listen on")
    parser.add_argument(
        "--host", "-H", default="127.0.0.1",
        help="The server host string; use 0.0.0.0 to allow all connections.")
    parser.add_argument(
        "--config", "-c", default='DevelopmentConfig', type=str,
        help="The configuration to use")
    parser.add_argument(
        "--config-file", "-f", type=str, default=None,
        help="The configuration file to use")
    parser.add_argument(
        "--tls", "-t", action="store_true", default=False,
        help="Start in TLS (https) mode (for Flask debug)")
    parser.add_argument(
        "--dont-use-reloader", default=False, action="store_true",
        help="Don't use the Flask reloader (for Flask debug)")
    parser.add_argument(
        "--debug", "-d", action='store_true', default=False,
        help="Runs the server using the gunicorn WSGI server.")
    cli.addVersionArgument(parser)
    cli.addDisableUrllibWarningsArgument(parser)


def runGunicornServer(parsedArgs):
    options = {
        'bind': '%s:%s' % (parsedArgs.host, parsedArgs.port),
        'workers': number_of_workers(),
        'accesslog': '-',  # Puts the access log on stdout
        'errorlog': '-'    # Puts the error log on stdout
    }
    app = StandaloneApplication(frontend.app, options)
    app.run()
    return app


def getServerParser():
    """
    Used by sphinx.argparse.
    """
    parser = common_cli.createArgumentParser("GA4GH reference server")
    addServerOptions(parser)
    return parser


def number_of_workers():
    return (multiprocessing.cpu_count() * 2) + 1


def server_main(args=None):
    parser = getServerParser()
    parsedArgs = parser.parse_args(args)
    if parsedArgs.disable_urllib_warnings:
        requests.packages.urllib3.disable_warnings()
    frontend.configure(
        parsedArgs.config_file, parsedArgs.config, parsedArgs.port)
    sslContext = None
    if parsedArgs.tls or ("OIDC_PROVIDER" in frontend.app.config):
        sslContext = "adhoc"
    if parsedArgs.debug:
        frontend.app.run(host=parsedArgs.host,
                         port=parsedArgs.port,
                         use_reloader=not parsedArgs.dont_use_reloader,
                         ssl_context=sslContext)
    else:
        runGunicornServer(parsedArgs)
