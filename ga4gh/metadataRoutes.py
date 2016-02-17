import flask
import json

import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol
import ga4gh.sqlbackend as sqlbackend

metadataRoutes = flask.Blueprint('metadataRoutes', __name__)

# metadata objects

MIMETYPE = "application/json"

metadataRecordTypes = (
    "biosample"
)
a = protocol.BioSample()


def validateMetadataCreateRequest(metadataType, metadataCreateRequest):
    """
    Performs validation of a create request given a metadata
    type. Raises bad request exceptions as needed.
    :param metadataType: string representing class
    :param metadataCreateRequest: new object to be created
    :return: boolean
    """
    # TODO protocol based validation
    if 'id' in metadataCreateRequest:
        # requests cant have an ID
        raise exceptions.BadRequestException
    if 'metadataType' in metadataCreateRequest and \
            metadataType != metadataCreateRequest['metadataType']:
        raise exceptions.BadRequestException
    metadataCreateRequest['metadataType'] = metadataType
    if metadataCreateRequest['metadataType'] not in metadataRecordTypes:
        raise exceptions.BadRequestException
    return True


def getMetadataSQL(metadataId):
    mb = sqlbackend.SQLBackend(protocol.BioSample, 'metadata.db')
    record = mb.getById(metadataId)
    return record


def insertMetadataSQL(metadataObject):
    mb = sqlbackend.SQLBackend(protocol.BioSample, 'metadata.db')
    return mb.insert(metadataObject)


def updateMetadataSQL(metadataObject):
    mb = sqlbackend.SQLBackend(protocol.BioSample, 'metadata.db')
    return mb.update(metadataObject)


def getMetadata(metadataType, metadataId):
    # get file by ID
    obj = getMetadataSQL(metadataId)
    # return string
    return flask.Response(json.dumps(obj.toJsonDict()), status=200, mimetype=MIMETYPE)


def createMetadata(metadataType, metadataCreateRequest):
    response = {}
    # validate that metadataType and metadataObject.type don't disagree
    validated = validateMetadataCreateRequest(
        metadataType, metadataCreateRequest)
    # add to internal representation
    if validated:
        # check to see if it exists first
        writeResponse = insertMetadataSQL(metadataCreateRequest)
        response[metadataType] = writeResponse.toJsonDict()
        response['message'] = "success"
        # add to sql lite
    else:
        response['message'] = "fail"
        # an exception that should never happen
    # return write response
    return flask.Response(json.dumps(response), status=201, mimetype=MIMETYPE)


def updateMetadata(metadataType, metadataId, metadataUpdateRequest):
    response = {}
    # TODO validate update request (dont collide ids from req and path)
    metadataUpdateRequest['id'] = metadataId
    # TODO pivot on metadataType

    updateResponse = updateMetadataSQL(metadataUpdateRequest)
    if updateResponse:
        response[metadataType] = updateResponse.toJsonDict()
        response['message'] = "success"
    else:
        response['message'] = "fail"
    # TODO whats the updated http code?
    return flask.Response(json.dumps(response), status=201, mimetype=MIMETYPE)


def deleteMetadata(metadataType, metadataId):
    response = {}
    # TODO
    return flask.Response(json.dumps(response), status=999)


def handleCreateMetadataRequest(metadataType, request, handler):
    requestData = json.loads(request.data)
    return handler(metadataType, requestData)


def handleGetMetadataRequest(metadataType, metadataId, handler):
    return handler(metadataType, metadataId)


def handleUpdateMetadataRequest(metadataType, metadataId, request, handler):
    requestData = json.loads(request.data)
    return handler(metadataType, metadataId, requestData)


def handleDeleteMetadataRequest(metadataType, metadataId, handler):
    return handler(metadataType, metadataId)


@metadataRoutes.route('/metadata/<metadataType>',
                      methods=['POST'])
def createMetadataRoute(metadataType):
    return handleCreateMetadataRequest(
        metadataType,
        flask.request,
        createMetadata)


@metadataRoutes.route('/metadata/<metadataType>/<metadataId>',
                      methods=['GET'])
def getMetadataRoute(metadataType, metadataId):
    return handleGetMetadataRequest(
        metadataType,
        metadataId,
        getMetadata)


@metadataRoutes.route('/metadata/<metadataType>/<metadataId>',
                      methods=['POST'])
def updateMetadataRoute(metadataType, metadataId):
    return handleUpdateMetadataRequest(
        metadataType,
        metadataId,
        flask.request,
        updateMetadata)


@metadataRoutes.route('/metadata/<metadataType>/<metadataId>',
                      methods=['DELETE'])
def deleteMetadataRoute(metadataType, metadataId):
    return handleUpdateMetadataRequest(
        metadataType,
        metadataId,
        deleteMetadata)
