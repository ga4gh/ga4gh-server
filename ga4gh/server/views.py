from . import app
from flask import Response, request, render_template, abort
import ga4gh.protocol as protocol
from ga4gh.backends import WormtableBackend, TabixBackend


# Helper function
# TODO: Does this need http://flask-cors.readthedocs.org/en/latest/
def handleHTTPPost(request, endpoint, protocolClass):
    data = request.get_data()
    try:
        protocolRequest = protocolClass.fromJSON(data)
    except ValueError:
        return "Bad Request", 404
    protocolResponse = endpoint(protocolRequest)
    ret = protocolResponse.toJSON()
    resp = Response(response=ret,
                    status=200,
                    mimetype="application/json")
    return resp


@app.route('/')
def index():
    return render_template("index.haml")


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


@app.route('/variantsets/search', methods=['POST'])
def searchVariantSets():
    return handleHTTPPost(request,
                          app.config["VariantBackend"].searchVariantSets,
                          protocol.GASearchVariantSetsRequest)


@app.route('/variants/search', methods=['POST'])
def searchVariants():
    return handleHTTPPost(request,
                          app.config["VariantBackend"].searchVariants,
                          protocol.GASearchVariantsRequest)
