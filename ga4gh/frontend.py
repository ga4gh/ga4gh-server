"""
The Flask frontend for the GA4GH API.

TODO Document properly.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import datetime

import flask
import flask.ext.cors as cors
import humanize
import werkzeug

import ga4gh
import ga4gh.backend as backend
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions

MIMETYPE = "application/json"
SEARCH_ENDPOINT_METHODS = ['POST', 'OPTIONS']


app = flask.Flask(__name__)


class NoConverter(werkzeug.routing.BaseConverter):
    """
    A converter that allows the routing matching algorithm to not
    match on certain literal terms

    This is needed because if there are e.g. two routes:

    /<version>/callsets/search
    /<version>/callsets/<id>

    A request for /someVersion/callsets/search will get routed to
    the second, which is not what we want.
    """
    def __init__(self, map, *items):
        werkzeug.routing.BaseConverter.__init__(self, map)
        self.items = items

    def to_python(self, value):
        if value in self.items:
            raise werkzeug.routing.ValidationError()
        return value


app.url_map.converters['no'] = NoConverter


class Version(object):
    """
    A major/minor/revision version tag
    """
    currentString = "current"

    @classmethod
    def isCurrentVersion(cls, versionString):
        if versionString == cls.currentString:
            return True
        return (Version.parseString(versionString) ==
                Version.parseString(protocol.version))

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

    @classmethod
    def getVersionForUrl(cls, versionString):
        """
        Returns the specfied version string in a form suitable for using
        within a URL. This involved prefixing with 'v'.
        """
        ret = versionString
        if not ret.startswith("v"):
            ret = "v{}".format(versionString)
        return ret


class ServerStatus(object):
    """
    Generates information about the status of the server for display
    """
    def __init__(self):
        self.startupTime = datetime.datetime.now()

    def getConfiguration(self):
        """
        Returns a list of configuration (key, value) tuples
        that are useful for users to view on an information page.
        Note that we should be careful here not to leak sensitive
        information. For example, keys and paths of data files should
        not be returned.
        """
        # TODO what other config keys are appropriate to export here?
        keys = [
            'DEBUG', 'REQUEST_VALIDATION', 'RESPONSE_VALIDATION',
            'DEFAULT_PAGE_SIZE', 'MAX_RESPONSE_LENGTH',
        ]
        return [(k, app.config[k]) for k in keys]

    def getPreciseUptime(self):
        """
        Returns the server precisely.
        """
        return self.startupTime.strftime("%H:%M:%S %d %b %Y")

    def getNaturalUptime(self):
        """
        Returns the uptime in a human-readable format.
        """
        return humanize.naturaltime(self.startupTime)

    def getProtocolVersion(self):
        """
        Returns the GA4GH protocol version we support.
        """
        return protocol.version

    def getServerVersion(self):
        """
        Returns the software version of this server.
        """
        return ga4gh.__version__

    def getUrls(self):
        """
        Returns the list of (httpMethod, URL) tuples that this server
        supports.
        """
        urls = []
        rules = app.url_map.iter_rules()
        excluded_methods = ('OPTIONS', 'HEAD')
        excluded_rules = (
            '/', '/flask-api/static/<path:filename>',
            '/static/<path:filename>')
        version = Version.getVersionForUrl(protocol.version)
        for rule in rules:
            for method in rule.methods:
                if (method not in excluded_methods and
                        rule.rule not in excluded_rules):
                    versionedRule = rule.rule.replace(
                        "<version>", version)
                    urls.append((versionedRule, method))
        urls.sort()
        return urls

    def getDatasetIds(self):
        """
        Returns the list of datasetIds for this backend
        """
        return app.backend.getDatasetIds()

    def getVariantSets(self, datasetId):
        """
        Returns the list of variant sets for the dataset
        """
        return app.backend.getDataset(datasetId).getVariantSets()

    def getReadGroupSets(self, datasetId):
        """
        Returns the list of ReadGroupSets for the dataset
        """
        return app.backend.getDataset(datasetId).getReadGroupSets()

    def getReferenceSets(self):
        """
        Returns the list of ReferenceSets for this server.
        """
        return app.backend.getReferenceSets()


def configure(configFile=None, baseConfig="ProductionConfig", extraConfig={}):
    """
    TODO Document this critical function! What does it do? What does
    it assume?
    """
    configStr = 'ga4gh.serverconfig:{0}'.format(baseConfig)
    app.config.from_object(configStr)
    if os.environ.get('GA4GH_CONFIGURATION') is not None:
        app.config.from_envvar('GA4GH_CONFIGURATION')
    if configFile is not None:
        app.config.from_pyfile(configFile)
    app.config.update(extraConfig.items())
    # Setup CORS
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
    theBackend.setDefaultPageSize(app.config["DEFAULT_PAGE_SIZE"])
    theBackend.setMaxResponseLength(app.config["MAX_RESPONSE_LENGTH"])
    app.backend = theBackend


def getFlaskResponse(responseString, httpStatus=200):
    """
    Returns a Flask response object for the specified data and HTTP status.
    """
    return flask.Response(responseString, status=httpStatus, mimetype=MIMETYPE)


def handleHttpPost(request, endpoint):
    """
    Handles the specified HTTP POST request, which maps to the specified
    protocol handler endpoint and protocol request class.
    """
    if request.mimetype != MIMETYPE:
        raise exceptions.UnsupportedMediaTypeException()
    responseStr = endpoint(request.get_data())
    return getFlaskResponse(responseStr)


def handleList(id_, endpoint, request):
    """
    Handles the specified HTTP GET request, mapping to a list request
    """
    responseStr = endpoint(id_, request.args)
    return getFlaskResponse(responseStr)


def handleHttpGet(id_, endpoint):
    """
    Handles the specified HTTP GET request, which maps to the specified
    protocol handler endpoint and protocol request class
    """
    responseStr = endpoint(id_)
    return getFlaskResponse(responseStr)


def handleHttpOptions():
    """
    Handles the specified HTTP OPTIONS request.
    """
    response = getFlaskResponse("")
    response.headers.add("Access-Control-Request-Methods", "GET,POST,OPTIONS")
    return response


@app.errorhandler(Exception)
def handleException(exception):
    """
    Handles an exception that occurs somewhere in the process of handling
    a request.
    """
    if app.config['DEBUG']:
        app.log_exception(exception)
    serverException = exception
    if not isinstance(exception, exceptions.BaseServerException):
        serverException = exceptions.getServerError(exception)
    responseStr = serverException.toProtocolElement().toJsonString()
    return getFlaskResponse(responseStr, serverException.httpStatus)


def assertCorrectVersion(version):
    if not Version.isCurrentVersion(version):
        raise exceptions.VersionNotSupportedException()


def handleFlaskGetRequest(version, id_, flaskRequest, endpoint):
    """
    Handles the specified flask request for one of the GET URLs
    at the specified version.  Invokes the specified endpoint to
    generate a response.
    """
    assertCorrectVersion(version)
    if flaskRequest.method == "GET":
        return handleHttpGet(id_, endpoint)
    else:
        raise exceptions.MethodNotAllowedException()


def handleFlaskListRequest(version, id_, flaskRequest, endpoint):
    """
    Handles the specified flask list request for one of the GET URLs
    at the specified version.  Invokes the specified endpoint to
    generate a response.
    """
    assertCorrectVersion(version)
    if flaskRequest.method == "GET":
        return handleList(id_, endpoint, flaskRequest)
    else:
        raise exceptions.MethodNotAllowedException()


def handleFlaskPostRequest(version, flaskRequest, endpoint):
    """
    Handles the specified flask request for one of the POST URLS
    at the specified version. Invokes the specified endpoint to
    generate a response.
    """
    assertCorrectVersion(version)
    if flaskRequest.method == "POST":
        return handleHttpPost(flaskRequest, endpoint)
    elif flaskRequest.method == "OPTIONS":
        return handleHttpOptions()
    else:
        raise exceptions.MethodNotAllowedException()


@app.route('/')
def index():
    return flask.render_template('index.html', info=app.serverStatus)


@app.route('/<version>')
def indexRedirect(version):
    try:
        isCurrentVersion = Version.isCurrentVersion(version)
    except TypeError:  # malformed "version string"
        raise exceptions.PathNotFoundException()
    if isCurrentVersion:
        return index()
    else:
        raise exceptions.PathNotFoundException()


@app.route('/<version>/references/<id>')
def getReference(version, id):
    return handleFlaskGetRequest(
        version, id, flask.request, app.backend.getReference)


@app.route('/<version>/referencesets/<id>')
def getReferenceSet(version, id):
    return handleFlaskGetRequest(
        version, id, flask.request, app.backend.getReferenceSet)


@app.route('/<version>/references/<id>/bases')
def listReferenceBases(version, id):
    return handleFlaskListRequest(
        version, id, flask.request, app.backend.listReferenceBases)


@app.route('/<version>/callsets/search', methods=SEARCH_ENDPOINT_METHODS)
def searchCallSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchCallSets)


@app.route('/<version>/readgroupsets/search', methods=SEARCH_ENDPOINT_METHODS)
def searchReadGroupSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReadGroupSets)


@app.route('/<version>/reads/search', methods=SEARCH_ENDPOINT_METHODS)
def searchReads(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReads)


@app.route('/<version>/referencesets/search', methods=SEARCH_ENDPOINT_METHODS)
def searchReferenceSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReferenceSets)


@app.route('/<version>/references/search', methods=SEARCH_ENDPOINT_METHODS)
def searchReferences(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReferences)


@app.route('/<version>/variantsets/search', methods=SEARCH_ENDPOINT_METHODS)
def searchVariantSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchVariantSets)


@app.route('/<version>/variants/search', methods=SEARCH_ENDPOINT_METHODS)
def searchVariants(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchVariants)


@app.route('/<version>/datasets/search', methods=SEARCH_ENDPOINT_METHODS)
def searchDatasets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchDatasets)


# The below paths have not yet been implemented


@app.route('/<version>/callsets/<no(search):id>')
def getCallset(version, id):
    raise exceptions.NotImplementedException()


@app.route('/<version>/alleles/<no(search):id>')
def getAllele(version, id):
    raise exceptions.NotImplementedException()


@app.route('/<version>/variants/<no(search):id>')
def getVariant(version, id):
    raise exceptions.NotImplementedException()


@app.route('/<version>/variantsets/<vsid>/sequences/<sid>')
def getVariantSetSequence(version, vsid, sid):
    raise exceptions.NotImplementedException()


@app.route('/<version>/variantsets/<no(search):id>')
def getVariantSet(version, id):
    raise exceptions.NotImplementedException()


@app.route('/<version>/feature/<id>')
def getFeature(version, id):
    raise exceptions.NotImplementedException()


@app.route('/<version>/sequences/<id>/bases')
def getSequenceBases(version, id):
    raise exceptions.NotImplementedException()


@app.route('/<version>/mode/<mode>')
def getMode(version, mode):
    raise exceptions.NotImplementedException()


@app.route('/<version>/datasets/<no(search):id>')
def getDataset(version, id):
    raise exceptions.NotImplementedException()


@app.route('/<version>/readgroupsets/<no(search):id>')
def getReadGroupSet(version, id):
    raise exceptions.NotImplementedException()


@app.route('/<version>/readgroups/<id>')
def getReadGroup(version, id):
    raise exceptions.NotImplementedException()


@app.route(
    '/<version>/genotypephenotype/search',
    methods=SEARCH_ENDPOINT_METHODS)
def searchGenotypePephenotype(version):
    raise exceptions.NotImplementedException()


@app.route('/<version>/individuals/search', methods=SEARCH_ENDPOINT_METHODS)
def searchIndividuals(version):
    raise exceptions.NotImplementedException()


@app.route('/<version>/samples/search', methods=SEARCH_ENDPOINT_METHODS)
def searchSamples(version):
    raise exceptions.NotImplementedException()


@app.route('/<version>/experiments/search', methods=SEARCH_ENDPOINT_METHODS)
def searchExperiments(version):
    raise exceptions.NotImplementedException()


@app.route(
    '/<version>/individualgroups/search',
    methods=SEARCH_ENDPOINT_METHODS)
def searchIndividualGroups(version):
    raise exceptions.NotImplementedException()


@app.route('/<version>/analyses/search', methods=SEARCH_ENDPOINT_METHODS)
def searchAnalyses(version):
    raise exceptions.NotImplementedException()


@app.route('/<version>/sequences/search', methods=SEARCH_ENDPOINT_METHODS)
def searchSequences(version):
    raise exceptions.NotImplementedException()


@app.route('/<version>/joins/search', methods=SEARCH_ENDPOINT_METHODS)
def searchJoins(version):
    raise exceptions.NotImplementedException()


@app.route('/<version>/subgraph/segments', methods=SEARCH_ENDPOINT_METHODS)
def subgraphSegments(version):
    raise exceptions.NotImplementedException()


@app.route('/<version>/subgraph/joins', methods=SEARCH_ENDPOINT_METHODS)
def subgraphJoins(version):
    raise exceptions.NotImplementedException()


@app.route('/<version>/features/search', methods=SEARCH_ENDPOINT_METHODS)
def searchFeatures(version):
    raise exceptions.NotImplementedException()


@app.route(
    '/<version>/variantsets/<id>/sequences/search',
    methods=SEARCH_ENDPOINT_METHODS)
def searchVariantSetSequences(version, id):
    raise exceptions.NotImplementedException()


@app.route('/<version>/alleles/search', methods=SEARCH_ENDPOINT_METHODS)
def searchAlleles(version):
    raise exceptions.NotImplementedException()


# The below methods ensure that JSON is returned for various errors
# instead of the default, html


@app.errorhandler(404)
def pathNotFoundHandler(errorString):
    return handleException(exceptions.PathNotFoundException())


@app.errorhandler(405)
def methodNotAllowedHandler(errorString):
    return handleException(exceptions.MethodNotAllowedException())
