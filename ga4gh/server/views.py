from . import app
from flask import abort, request, render_template, Response
from flask.ext.api import exceptions
import ga4gh.protocol as protocol


def handleHTTPPost(request, endpoint, protocolClass):
    """
    Handles the specified HTTP POST request, which maps to the specified
    protocol handler handpoint and protocol request class.
    """
    if request.mimetype != "application/json":
        raise exceptions.UnsupportedMediaType()
    content_length_ = request.content_length
    data = request.get_data()
    # TODO this should be a more specific Exception for JSON
    # parse errors; malformed JSON input is a HTTP error, whereas
    # anything after this gives a HTTP success, with a GAException
    # response.
    try:
        protocolRequest = protocolClass.fromJSON(data)
    except ValueError:
        raise exceptions.ParseError()
    protocolResponse = endpoint(protocolRequest)
    s = protocolResponse.toJSON()
    response = Response(s, status=200, mimetype="application/json")
    return response


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
        return handleHTTPPost(
            request,
            app.config["VariantBackend"].searchVariantSets,
            protocol.GASearchVariantSetsRequest)
    else:
        return handleHTTPOptions()


@app.route('/variants/search', methods=['POST', 'OPTIONS'])
def searchVariants():
    if request.method == "POST":
        return handleHTTPPost(
            request,
            app.config["VariantBackend"].searchVariants,
            protocol.GASearchVariantsRequest)
    else:
        return handleHTTPOptions()
