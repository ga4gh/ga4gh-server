from . import app
from flask import abort, request, render_template, Response
import ga4gh.protocol as protocol


def handleHTTPPost(request, endpoint, protocolClass):
    """
    Handles the specified HTTP POST request, which maps to the specified
    protocol handler handpoint and protocol request class.
    """
    if request.mimetype != "application/json":
        abort(415) # Unsupported Media Type.
    # TODO: update app.config with user configurable MAX_CONTENT_LENGTH.
    # Setting to 2MB for now.
    content_length_ = request.content_length
    if content_length_ is not None and content_length_ > 2 * 1024 * 1024:
	abort(413) # Request Entity Too Large.
    data = request.get_data()
    try:
        protocolRequest = protocolClass.fromJSON(data)
    except ValueError:
        abort(400)
    protocolResponse = endpoint(protocolRequest)
    s = protocolResponse.toJSON()
    return(s, status=200, mimetype="application/json")


def handleHTTPOptions():
    """
    Handles the specified HTTP OPTIONS request.
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
