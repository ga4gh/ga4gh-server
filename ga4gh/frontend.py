"""
The Flask frontend for the GA4GH API.

TODO Document properly.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import traceback
import datetime

import flask
import flask.ext.api as api
import flask.ext.cors as cors
import humanize

import ga4gh.frontend_exceptions as frontendExceptions
import ga4gh.protocol as protocol


app = flask.Flask(__name__)


class ServerStatus(object):

    def __init__(self):
        self.startupTime = datetime.datetime.now()

    def getStatusInfo(self):
        info = {
            'urls': self._getUrls(),
            'uptime': self._getUptime(),
            'version': protocol.version,
        }
        return info

    def _getUptime(self):
        uptime = {
            'natural': humanize.naturaltime(self.startupTime),
            'datetime': self.startupTime.strftime("%H:%M:%S %d %b %Y")
        }
        return uptime

    def _getUrls(self):
        urls = []
        rules = app.url_map.iter_rules()
        excluded_methods = ('OPTIONS', 'HEAD')
        excluded_rules = (
            '/', '/flask-api/static/<path:filename>',
            '/static/<path:filename>')
        for rule in rules:
            for method in rule.methods:
                if (method not in excluded_methods and
                        rule.rule not in excluded_rules):
                        urls.append((rule.rule, method))
        urls.sort()
        return urls


def configure(config="DefaultConfig", config_file=None):
    configStr = 'ga4gh.serverconfig:{0}'.format(config)
    app.config.from_object(configStr)
    if os.environ.get('GA4GH_CONFIGURATION') is not None:
        app.config.from_envvar('GA4GH_CONFIGURATION')
    if config_file is not None:
        app.config.from_pyfile(config_file)
    cors.CORS(app, allow_headers='Content-Type')
    app.serverStatus = ServerStatus()


def handleHttpPost(request, endpoint):
    """
    Handles the specified HTTP POST request, which maps to the specified
    protocol handler handpoint and protocol request class.
    """
    mimetype = "application/json"
    if request.mimetype != mimetype:
        raise api.exceptions.UnsupportedMediaType()
    responseStr = endpoint(request.get_data())
    return flask.Response(responseStr, status=200, mimetype=mimetype)


def handleHttpOptions():
    """
    Handles the specified HTTP OPTIONS request.
    """
    response = flask.Response("", mimetype="application/json")
    response.headers.add("Access-Control-Request-Methods", "GET,POST,OPTIONS")
    return response


@app.errorhandler(Exception)
def handleException(exception):
    exceptionClass = exception.__class__
    if exceptionClass in frontendExceptions.exceptionMap:
        newExceptionClass = frontendExceptions.exceptionMap[exceptionClass]
        exception = newExceptionClass()

    if not isinstance(exception, frontendExceptions.FrontendException):
        if app.config['DEBUG']:
            print(traceback.format_exc(exception))
        exception = frontendExceptions.ServerException()

    error = protocol.GAException()
    error.errorCode = exception.code
    error.message = exception.message
    response = flask.jsonify(error.toJsonDict())
    response.status_code = exception.httpStatus
    return response


@app.route('/')
def index():
    return flask.render_template(
        'index.html', info=app.serverStatus.getStatusInfo())


@app.route('/references/<id>', methods=['GET'])
def getReference(id):
    raise frontendExceptions.PathNotFoundException()


@app.route('/references/<id>/bases', methods=['GET'])
def getReferenceBases(id):
    raise frontendExceptions.PathNotFoundException()


@app.route('/referencesets/<id>', methods=['GET'])
def getReferenceSet(id):
    raise frontendExceptions.PathNotFoundException()


@app.route('/callsets/search', methods=['POST'])
def searchCallSets():
    raise frontendExceptions.PathNotFoundException()


@app.route('/readgroupsets/search', methods=['POST'])
def searchReadGroupSets():
    raise frontendExceptions.PathNotFoundException()


@app.route('/reads/search', methods=['POST'])
def searchReads():
    raise frontendExceptions.PathNotFoundException()


@app.route('/referencesets/search', methods=['POST'])
def searchReferenceSets():
    raise frontendExceptions.PathNotFoundException()


@app.route('/references/search', methods=['POST'])
def searchReferences():
    raise frontendExceptions.PathNotFoundException()


@app.route('/variantsets/search', methods=['POST', 'OPTIONS'])
def searchVariantSets():
    if flask.request.method == "POST":
        return handleHttpPost(flask.request, app.backend.searchVariantSets)
    else:
        return handleHttpOptions()


@app.route('/variants/search', methods=['POST', 'OPTIONS'])
def searchVariants():
    if flask.request.method == "POST":
        return handleHttpPost(flask.request, app.backend.searchVariants)
    else:
        return handleHttpOptions()
