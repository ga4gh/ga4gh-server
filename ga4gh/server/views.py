from . import app
from flask import Response, request, render_template
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
    return "TODO"


@app.route('/references/<id>/bases', methods=['GET'])
def getReferenceBases(id):
    return "TODO"


@app.route('/referencesets/<id>', methods=['GET'])
def getReferenceSet(id):
    return "TODO"


@app.route('/callsets/search', methods=['POST'])
def searchCallSets():
    return "TODO"


@app.route('/readgroupsets/search', methods=['POST'])
def searchReadGroupSets():
    return "TODO"


@app.route('/reads/search', methods=['POST'])
def searchReads():
    return "TODO"


@app.route('/referencesets/search', methods=['POST'])
def searchReferenceSets():
    return "TODO"


@app.route('/references/search', methods=['POST'])
def searchReferences():
    return "TODO"


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
