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
import flask.ext.cors as cors
import humanize

import ga4gh.frontend_exceptions as frontendExceptions
import ga4gh.backend as backend
import ga4gh.protocol as protocol


app = flask.Flask(__name__)


class Version(object):
    """
    A major/minor/revision version tag
    """
    @classmethod
    def parseString(cls, versionString):
        versions = versionString.strip('vV').split('.')
        return Version(*versions)

    def __init__(self, major, minor, revision):
        self.version = (major, minor, revision)

    def __cmp__(self, other):
        return cmp(self.version, other.version)

    def __hash__(self):
        return hash(self.version)

    def __eq__(self, other):
        return self.version == other.version

    def __ne__(self, other):
        return not self.__eq__(other)


class ServerStatus(object):
    """
    Generates information about the status of the server for display
    """
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


def configure(config="DefaultConfig", configFile=None):
    configStr = 'ga4gh.serverconfig:{0}'.format(config)
    app.config.from_object(configStr)
    if os.environ.get('GA4GH_CONFIGURATION') is not None:
        app.config.from_envvar('GA4GH_CONFIGURATION')
    if configFile is not None:
        app.config.from_pyfile(configFile)
    cors.CORS(app, allow_headers='Content-Type')
    app.serverStatus = ServerStatus()
    # Allocate the backend
    # TODO is this a good way to determine what type of backend we should
    # instantiate? We should think carefully about this. The approach of
    # using the special strings __SIMULATED__ and __EMPTY__ seems OK for
    # now, but is certainly not ideal.
    dataSource = app.config["DATA_SOURCE"]
    if dataSource == "__SIMULATED__":
        randomSeed = app.config["SIMULATED_BACKEND_RANDOM_SEED"]
        numCalls = app.config["SIMULATED_BACKEND_NUM_CALLS"]
        variantDensity = app.config["SIMULATED_BACKEND_VARIANT_DENSITY"]
        numVariantSets = app.config["SIMULATED_BACKEND_NUM_VARIANT_SETS"]
        theBackend = backend.SimulatedBackend(
            randomSeed, numCalls, variantDensity, numVariantSets)
    elif dataSource == "__EMPTY__":
        theBackend = backend.EmptyBackend()
    else:
        theBackend = backend.FileSystemBackend(dataSource)
    theBackend.setRequestValidation(app.config["REQUEST_VALIDATION"])
    theBackend.setResponseValidation(app.config["RESPONSE_VALIDATION"])
    app.backend = theBackend


def handleHttpPost(request, endpoint):
    """
    Handles the specified HTTP POST request, which maps to the specified
    protocol handler handpoint and protocol request class.
    """
    mimetype = "application/json"
    if request.mimetype != mimetype:
        raise frontendExceptions.UnsupportedMediaTypeException()
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
    # if the caught exception implements toErrorMessage, extract the
    # message to be returned to the client at this point (because
    # the exception object may get overwritten later in the method)
    errorMessage = None
    if hasattr(exception, 'toErrorMessage'):
        errorMessage = exception.toErrorMessage()

    # if the type of exception is one we have registered in the exceptionMap,
    # convert the exception to one that we want to return to the client
    exceptionClass = exception.__class__
    if exceptionClass in frontendExceptions.exceptionMap:
        newExceptionClass = frontendExceptions.exceptionMap[exceptionClass]
        exception = newExceptionClass()

    # if at this point the exception is still unrecognized
    # send the client a 500 error
    if not isinstance(exception, frontendExceptions.FrontendException):
        if app.config['DEBUG']:
            print(traceback.format_exc(exception))
        exception = frontendExceptions.ServerException()

    # serialize the exception
    error = exception.toGaException()
    if errorMessage is not None:
        error.message = errorMessage
    response = flask.jsonify(error.toJsonDict())
    response.status_code = exception.httpStatus
    return response


def handleFlaskPostRequest(version, flaskRequest, endpoint):
    if Version.parseString(version) != Version.parseString(protocol.version):
        raise frontendExceptions.VersionNotSupportedException()

    if flaskRequest.method == "POST":
        return handleHttpPost(flaskRequest, endpoint)
    elif flaskRequest.method == "OPTIONS":
        return handleHttpOptions()
    else:
        raise frontendExceptions.MethodNotAllowedException()


@app.route('/')
def index():
    return flask.render_template(
        'index.html', info=app.serverStatus.getStatusInfo())


@app.route('/<version>/references/<id>', methods=['GET'])
def getReference(version, id):
    raise frontendExceptions.PathNotFoundException()


@app.route('/<version>/references/<id>/bases', methods=['GET'])
def getReferenceBases(version, id):
    raise frontendExceptions.PathNotFoundException()


@app.route('/<version>/referencesets/<id>', methods=['GET'])
def getReferenceSet(version, id):
    raise frontendExceptions.PathNotFoundException()


@app.route('/<version>/callsets/search', methods=['POST'])
def searchCallSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchCallSets)


@app.route('/<version>/readgroupsets/search', methods=['POST'])
def searchReadGroupSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReadGroupSets)


@app.route('/<version>/reads/search', methods=['POST'])
def searchReads(version):
    raise frontendExceptions.PathNotFoundException()


@app.route('/<version>/referencesets/search', methods=['POST'])
def searchReferenceSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReferenceSets)


@app.route('/<version>/references/search', methods=['POST'])
def searchReferences(version):
    raise frontendExceptions.PathNotFoundException()


@app.route('/<version>/variantsets/search', methods=['POST', 'OPTIONS'])
def searchVariantSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchVariantSets)


@app.route('/<version>/variants/search', methods=['POST', 'OPTIONS'])
def searchVariants(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchVariants)


# The below methods ensure that JSON is returned for various errors
# instead of the default, html


@app.errorhandler(404)
def pathNotFoundHandler(errorString):
    return handleHttpError(frontendExceptions.PathNotFoundException())


@app.errorhandler(405)
def methodNotAllowedHandler(errorString):
    return handleHttpError(frontendExceptions.MethodNotAllowedException())


def handleHttpError(exception):
    error = exception.toGaException()
    response = flask.jsonify(error.toJsonDict())
    response.status_code = exception.httpStatus
    return response
