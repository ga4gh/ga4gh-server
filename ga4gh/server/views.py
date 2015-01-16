"""
The Flask views for the GA4GH HTTP API.

TODO Document properly.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from . import app
from flask import abort, request, Response
from flask.ext.api import exceptions


def handleHTTPPost(request, endpoint):
    """
    Handles the specified HTTP POST request, which maps to the specified
    protocol handler handpoint and protocol request class.
    """
    mimetype = "application/json"
    if request.mimetype != mimetype:
        raise exceptions.UnsupportedMediaType()
    responseStr = endpoint(request.get_data())
    return Response(responseStr, status=200, mimetype=mimetype)


def handleHTTPOptions():
    """
    Handles the specified HTTP OPTIONS request.
    """
    response = Response("", mimetype="application/json")
    response.headers.add("Access-Control-Request-Methods", "GET,POST,OPTIONS")
    return response


@app.route('/')
def index():
    abort(404)


@app.route('/references/<id>', methods=['GET'])
def getReference(id):
    abort(404)


@app.route('/references/<id>/bases', methods=['GET'])
def getReferenceBases(id):
    abort(404)


@app.route('/referencesets/<id>', methods=['GET'])
def getReferenceSet(id):
    abort(404)


@app.route('/callsets/search', methods=['POST'])
def searchCallSets():
    abort(404)


@app.route('/readgroupsets/search', methods=['POST'])
def searchReadGroupSets():
    abort(404)


@app.route('/reads/search', methods=['POST'])
def searchReads():
    abort(404)


@app.route('/referencesets/search', methods=['POST'])
def searchReferenceSets():
    abort(404)


@app.route('/references/search', methods=['POST'])
def searchReferences():
    abort(404)


@app.route('/variantsets/search', methods=['POST', 'OPTIONS'])
def searchVariantSets():
    if request.method == "POST":
        return handleHTTPPost(request, app.backend.searchVariantSets)
    else:
        return handleHTTPOptions()


@app.route('/variants/search', methods=['POST', 'OPTIONS'])
def searchVariants():
    if request.method == "POST":
        return handleHTTPPost(request, app.backend.searchVariants)
    else:
        return handleHTTPOptions()
