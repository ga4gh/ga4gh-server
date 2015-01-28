"""
The Flask frontend for the GA4GH API.

TODO Document properly.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import flask
import flask.ext.api as api
import flask.ext.cors as cors

app = api.FlaskAPI(__name__)


def configure(config="DefaultConfig", config_file=None):
    configStr = 'ga4gh.serverconfig:{0}'.format(config)
    app.config.from_object(configStr)
    if os.environ.get('GA4GH_CONFIGURATION') is not None:
        app.config.from_envvar('GA4GH_CONFIGURATION')
    if config_file is not None:
        app.config.from_pyfile(args.config_file)
    cors.CORS(app, allow_headers='Content-Type')


def handleHTTPPost(request, endpoint):
    """
    Handles the specified HTTP POST request, which maps to the specified
    protocol handler handpoint and protocol request class.
    """
    mimetype = "application/json"
    if request.mimetype != mimetype:
        raise api.exceptions.UnsupportedMediaType()
    responseStr = endpoint(request.get_data())
    return flask.Response(responseStr, status=200, mimetype=mimetype)


def handleHTTPOptions():
    """
    Handles the specified HTTP OPTIONS request.
    """
    response = flask.Response("", mimetype="application/json")
    response.headers.add("Access-Control-Request-Methods", "GET,POST,OPTIONS")
    return response


@app.route('/')
def index():
    flask.abort(404)


@app.route('/references/<id>', methods=['GET'])
def getReference(id):
    flask.abort(404)


@app.route('/references/<id>/bases', methods=['GET'])
def getReferenceBases(id):
    flask.abort(404)


@app.route('/referencesets/<id>', methods=['GET'])
def getReferenceSet(id):
    flask.abort(404)


@app.route('/callsets/search', methods=['POST'])
def searchCallSets():
    flask.abort(404)


@app.route('/readgroupsets/search', methods=['POST'])
def searchReadGroupSets():
    flask.abort(404)


@app.route('/reads/search', methods=['POST'])
def searchReads():
    flask.abort(404)


@app.route('/referencesets/search', methods=['POST'])
def searchReferenceSets():
    flask.abort(404)


@app.route('/references/search', methods=['POST'])
def searchReferences():
    flask.abort(404)


@app.route('/variantsets/search', methods=['POST', 'OPTIONS'])
def searchVariantSets():
    if flask.request.method == "POST":
        return handleHTTPPost(flask.request, app.backend.searchVariantSets)
    else:
        return handleHTTPOptions()


@app.route('/variants/search', methods=['POST', 'OPTIONS'])
def searchVariants():
    if flask.request.method == "POST":
        return handleHTTPPost(flask.request, app.backend.searchVariants)
    else:
        return handleHTTPOptions()
