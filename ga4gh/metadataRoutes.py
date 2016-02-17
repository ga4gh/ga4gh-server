import flask
import json
import sqlite3
import uuid

import ga4gh.exceptions as exceptions

metadataRoutes = flask.Blueprint('metadataRoutes', __name__)

# metadata objects

MIMETYPE = "application/json"

metadataRecordTypes = (
    "biosample"
)

class MetadataBackend():
    # TODO generalize for the rest of the API?
    columns = ("id", "name", "description")
    tableName = "biosample"

    def __init__(self, dbFile):
        self._dbconn = sqlite3.connect(dbFile)
        dbconn = self._dbconn
        dbcur = dbconn.cursor()
        # TODO generalize table creation based on
        dbcur.execute("create table if not exists " + self.tableName + " (id text, name text, description text)")
        dbcur.close()

    def hydrate(self, row):
        return dict(zip(self.columns, row))

    def getById(self, id):
        sql = "SELECT * FROM biosample where id=='" + id + "';"
        query = self._dbconn.execute(sql)
        results = query.fetchall()
        if results:
            return self.hydrate(results[0])
        else:
            return None

    def delete(self, obj):
        # TODO
        pass

    def insert(self, obj):
        # generate ID
        obj['id'] = str(uuid.uuid4())
        if 'description' not in obj:
            obj['description'] = ""
        if 'name' not in obj:
            obj['name'] = ""
        # TODO generalize this line and clean inputs
        sql = "INSERT INTO " + self.tableName + " VALUES('" + obj['id'] + "', '" + obj['name'] + "', '" + obj['description'] + "')"
        query = self._dbconn.execute(sql)
        self._dbconn.commit()
        return self.getById(obj['id']) # TODO add message using cursor dereferencing

    def update(self, obj):
        # get the record we want to update
        old = self.getById(obj['id'])
        if not old:
            # there should be an existing record to update
            raise exceptions.BadRequestException
        # overwrite dict values with new
        for k in old.keys():
            if k not in obj:
                obj[k] = old[k]
        # TODO generalize this line and clean inputs
        sql = "UPDATE " + self.tableName + " SET name = '" + obj['name'] + "', description = '" + obj['description'] + "' WHERE id=='" + obj['id'] + "'"
        print(sql)
        query = self._dbconn.execute(sql)
        self._dbconn.commit()
        # update SQL
        return self.getById(obj['id'])


def validateMetadataCreateRequest(metadataType, metadataCreateRequest):
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
    mb = MetadataBackend('metadata.db')
    record = mb.getById(metadataId)
    return record

def insertMetadataSQL(metadataObject):
    mb = MetadataBackend('metadata.db')
    return mb.insert(metadataObject)

def updateMetadataSQL(metadataObject):
    mb = MetadataBackend('metadata.db')
    return mb.update(metadataObject)


def getMetadata(metadataType, metadataId):
    # get file by ID
    obj = getMetadataSQL(metadataId)
    # parse into JSON
    print(obj)
    # return string
    return flask.Response(json.dumps(obj), status=200, mimetype=MIMETYPE)

def createMetadata(metadataType, metadataCreateRequest):
    response = {}
    # validate that metadataType and metadataObject.type don't disagree
    validated = validateMetadataCreateRequest(metadataType, metadataCreateRequest)
    
    # add to internal representation
    if validated:
        # check to see if it exists first
        writeResponse = insertMetadataSQL(metadataCreateRequest)
        response[metadataType] = writeResponse
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
        response[metadataType] = updateResponse
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

@metadataRoutes.route('/metadata/<metadataType>', methods=['POST'])
def createMetadataRoute(metadataType):
    return handleCreateMetadataRequest(
        metadataType,
        flask.request,
        createMetadata)
        
@metadataRoutes.route('/metadata/<metadataType>/<metadataId>', methods=['GET'])
def getMetadataRoute(metadataType, metadataId):
    return handleGetMetadataRequest(
        metadataType,
        metadataId,
        getMetadata)

@metadataRoutes.route('/metadata/<metadataType>/<metadataId>', methods=['POST'])
def updateMetadataRoute(metadataType, metadataId):
    return handleUpdateMetadataRequest(
        metadataType,
        metadataId,
        flask.request,
        updateMetadata)

@metadataRoutes.route('/metadata/<metadataType>/<metadataId>', methods=['POST'])
def deleteMetadataRoute(metadataType, metadataId):
    return handleUpdateMetadataRequest(
        metadataType,
        metadataId,
        deleteMetadata)