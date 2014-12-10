from . import app
from flask import request, render_template, abort
import ga4gh.protocol as protocol
import werkzeug.wrappers as wzw
import werkzeug.exceptions as wze


def handleHTTPPost(request, endpoint, protocolClass):
    """
    Handles the specified HTTP POST request, which maps to the specified
    protocol handler handpoint and protocol request class.
    """
    if request.mimetype != "application/json":
        # TODO is this the correct HTTP response?
        raise wze.UnsupportedMediaType()
    # Make sure we don't get tricked into reading in large volumes
    # of data, exhausting server memory
    if request.content_length >= 8192:  # FIXME!
        raise wze.RequestEntityTooLarge()
    data = request.get_data()
    # TODO this should be a more specific Exception for JSON
    # parse errors; malformed JSON input is a HTTP error, whereas
    # anything after this gives a HTTP success, with a GAException
    # response.
    try:
        protocolRequest = protocolClass.fromJSON(data)
    except ValueError:
        raise wze.BadRequest()
    protocolResponse = endpoint(protocolRequest)
    s = protocolResponse.toJSON()
    response = wzw.Response(s, mimetype="application/json")
    # TODO is this correct CORS support?
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


def handleHTTPOptions():
    """
    Handles the specified HTTP OPTIONS request returing a werkzeug
    response.
    """
    response = wzw.Response("", mimetype="application/json")
    # TODO is this correct CORS support?
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add(
        "Access-Control-Request-Methods", "GET,POST,OPTIONS")
    response.headers.add("Access-Control-Allow-Headers", "Content-Type")
    return response


@app.route('/')
def index():
    return "ATGC"


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
