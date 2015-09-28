"""
Client classes for the GA4GH reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import requests
import posixpath
import logging

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


class HttpClient(object):
    """
    GA4GH Http Client
    """

    def __init__(self, urlPrefix, debugLevel=0, key=None):
        self._urlPrefix = urlPrefix
        self._debugLevel = debugLevel
        self._bytesRead = 0
        self._key = key
        self._pageSize = None

        # logging config
        # TODO we need to revisit this logging setup so that we can
        # disentangle our logs from urllib3's.
        logging.basicConfig()
        self._logger = logging.getLogger(__name__)
        if self._debugLevel == 0:
            logLevel = logging.WARN
        elif self._debugLevel == 1:
            logLevel = logging.INFO
        else:
            logLevel = logging.DEBUG
        self._logger.setLevel(logLevel)

        requestsLog = logging.getLogger("requests.packages.urllib3")
        requestsLog.setLevel(logLevel)
        if self._debugLevel == 0:
            # suppress warning about using https without cert verification
            requests.packages.urllib3.disable_warnings()
        requestsLog.propagate = True

    def getPageSize(self):
        """
        Returns the suggested maximum size of pages of results returned by
        the server.
        """
        return self._pageSize

    def getBytesRead(self):
        """
        Returns the total number of (non HTTP) bytes read from the server
        by this client.
        """
        return self._bytesRead

    def setPageSize(self, pageSize):
        """
        Sets the requested maximum size of pages of results returned by the
        server to the specified value.
        """
        self._pageSize = pageSize

    def _getAuth(self):
        return {'key': self._key}

    # Ordinarily logger's implementation will take care of if log messages
    # should be emitted based on the log level.  The _shouldLog* methods
    # are only used if there are additional performance hits involved in
    # creating the log message that we want to avoid otherwise.
    def _shouldLogDebug(self):
        return self._debugLevel > 1

    def _shouldLogInfo(self):
        return self._debugLevel > 0

    def _debugResponse(self, jsonString):
        if self._shouldLogDebug():
            self._logger.debug("json response:")
            prettyString = self._prettyJsonString(jsonString)
            self._logger.debug(prettyString)

    def _debugRequest(self, jsonString):
        if self._shouldLogDebug():
            self._logger.debug("json request:")
            prettyString = self._prettyJsonString(jsonString)
            self._logger.debug(prettyString)

    def _prettyJsonString(self, jsonString):
        # note: expensive method
        return json.dumps(json.loads(jsonString), sort_keys=True, indent=4)

    def _checkStatus(self, response):
        if response.status_code != requests.codes.ok:
            self._logger.error("%s %s", response.status_code, response.text)
            raise exceptions.RequestNonSuccessException(
                "Url {0} had status_code {1}".format(
                    response.url, response.status_code))

    def _updateBytesRead(self, jsonString):
        self._bytesRead += len(jsonString)

    def _deserializeResponse(self, response, protocolResponseClass):
        jsonResponseString = response.text
        self._updateBytesRead(jsonResponseString)
        self._debugResponse(jsonResponseString)
        if jsonResponseString == '':
            raise exceptions.EmptyResponseException()
        responseObject = protocolResponseClass.fromJsonString(
            jsonResponseString)
        return responseObject

    def _doRequest(self, httpMethod, url, protocolResponseClass,
                   httpParams={}, httpData=None):
        """
        Performs a request to the server and returns the response
        """
        headers = {}
        params = self._getAuth()
        params.update(httpParams)
        self._logger.info("{0} {1}".format(httpMethod, url))
        if httpData is not None:
            headers.update({"Content-type": "application/json"})
            self._debugRequest(httpData)
        response = requests.request(
            httpMethod, url, params=params, data=httpData, headers=headers,
            verify=False)
        self._checkStatus(response)
        return self._deserializeResponse(response, protocolResponseClass)

    def runSearchRequest(self, protocolRequest, objectName,
                         protocolResponseClass):
        """
        Runs the specified request at the specified objectName and instantiates
        an object of the specified class. We yield each object in listAttr.
        If pages of results are present, repeat this process until the
        pageToken is null.
        """
        fullUrl = posixpath.join(self._urlPrefix, objectName + '/search')
        notDone = True
        while notDone:
            data = protocolRequest.toJsonString()
            responseObject = self._doRequest(
                'POST', fullUrl, protocolResponseClass, httpData=data)
            valueList = getattr(
                responseObject, protocolResponseClass.getValueListName())
            self._logger.info("Response pageSize={}".format(len(valueList)))
            for extract in valueList:
                yield extract
            notDone = responseObject.nextPageToken is not None
            protocolRequest.pageToken = responseObject.nextPageToken

    def listReferenceBases(self, id_, start=None, end=None):
        """
        Returns an iterator over the bases from the server in the form
        of consecutive strings. This command does not conform to the
        patterns of the other search and get requests, and is implemented
        differently.
        """
        url = "references/{id}/bases"
        fullUrl = posixpath.join(self._urlPrefix, url).format(id=id_)
        request = protocol.ListReferenceBasesRequest()
        request.start = start
        request.end = end
        notDone = True
        while notDone:
            response = self._doRequest(
                'GET', fullUrl, protocol.ListReferenceBasesResponse,
                request.toJsonDict())
            self._logger.info("Response pageSize={}".format(
                len(response.sequence)))
            yield response.sequence
            notDone = response.nextPageToken is not None
            request.pageToken = response.nextPageToken

    def runGetRequest(self, objectName, protocolResponseClass, id_):
        """
        Requests an object from the server and returns the object of
        type protocolResponseClass that has id id_.
        Used for requests where a single object is the expected response.
        """
        url = "{objectName}/{id}"
        fullUrl = posixpath.join(
            self._urlPrefix, url).format(id=id_, objectName=objectName)
        return self._doRequest('GET', fullUrl, protocolResponseClass)

    def getDataset(self, id_):
        """
        Returns a dataset from the server
        """
        return self.runGetRequest("datasets", protocol.Dataset, id_)

    def getReferenceSet(self, id_):
        """
        Returns a referenceSet from the server
        """
        return self.runGetRequest("referencesets", protocol.ReferenceSet, id_)

    def getReference(self, id_):
        """
        Returns a reference from the server
        """
        return self.runGetRequest("references", protocol.Reference, id_)

    def getReadGroupSet(self, id_):
        """
        Returns a read group set from the server
        """
        return self.runGetRequest(
            "readgroupsets", protocol.ReadGroupSet, id_)

    def getReadGroup(self, id_):
        """
        Returns a read group from the server
        """
        return self.runGetRequest("readgroups", protocol.ReadGroup, id_)

    def getCallset(self, id_):
        """
        Returns a callset from the server
        """
        return self.runGetRequest("callsets", protocol.CallSet, id_)

    def getVariant(self, id_):
        """
        Returns a variant from the server
        """
        return self.runGetRequest("variants", protocol.Variant, id_)

    def searchVariants(
            self, variantSetId, start=None, end=None, referenceName=None,
            callSetIds=None):
        """
        Returns an iterator over the Variants from the server
        """
        request = protocol.SearchVariantsRequest()
        request.referenceName = referenceName
        request.start = start
        request.end = end
        request.variantSetId = variantSetId
        request.callSetIds = callSetIds
        request.pageSize = self._pageSize
        return self.runSearchRequest(
            request, "variants", protocol.SearchVariantsResponse)

    def getVariantSet(self, id_):
        """
        Returns a variantSet from the server
        """
        return self.runGetRequest("variantsets", protocol.VariantSet, id_)

    def searchVariantSets(self, datasetId):
        """
        Returns an iterator over the VariantSets on the server. If datasetId
        is specified, return only the VariantSets in this dataset.
        """
        request = protocol.SearchVariantSetsRequest()
        request.datasetId = datasetId
        request.pageSize = self._pageSize
        return self.runSearchRequest(
            request, "variantsets", protocol.SearchVariantSetsResponse)

    def searchReferenceSets(
            self, accession=None, md5checksum=None, assemblyId=None):
        """
        Returns an iterator over the ReferenceSets from the server.
        """
        request = protocol.SearchReferenceSetsRequest()
        request.accession = accession
        request.md5checksum = md5checksum
        request.assemblyId = assemblyId
        request.pageSize = self._pageSize
        return self.runSearchRequest(
            request, "referencesets", protocol.SearchReferenceSetsResponse)

    def searchReferences(
            self, referenceSetId, accession=None, md5checksum=None):
        """
        Returns an iterator over the References from the server
        """
        request = protocol.SearchReferencesRequest()
        request.referenceSetId = referenceSetId
        request.accession = accession
        request.md5checksum = md5checksum
        request.pageSize = self._pageSize
        return self.runSearchRequest(
            request, "references", protocol.SearchReferencesResponse)

    def searchCallSets(self, variantSetId, name=None):
        """
        Returns an iterator over the CallSets from the server
        """
        request = protocol.SearchCallSetsRequest()
        request.variantSetId = variantSetId
        request.name = name
        request.pageSize = self._pageSize
        return self.runSearchRequest(
            request, "callsets", protocol.SearchCallSetsResponse)

    def searchReadGroupSets(self, datasetId, name=None):
        """
        Returns an iterator over the ReadGroupSets from the server
        """
        request = protocol.SearchReadGroupSetsRequest()
        request.datasetId = datasetId
        request.name = name
        request.pageSize = self._pageSize
        return self.runSearchRequest(
            request, "readgroupsets", protocol.SearchReadGroupSetsResponse)

    def searchReads(
            self, readGroupIds, referenceId=None, start=None, end=None):
        """
        Returns an iterator over the Reads from the server
        """
        request = protocol.SearchReadsRequest()
        request.readGroupIds = readGroupIds
        request.referenceId = referenceId
        request.start = start
        request.end = end
        request.pageSize = self._pageSize
        return self.runSearchRequest(
            request, "reads", protocol.SearchReadsResponse)

    def searchDatasets(self):
        """
        Returns an iterator over the Datasets from the server
        """
        request = protocol.SearchDatasetsRequest()
        request.pageSize = self._pageSize
        return self.runSearchRequest(
            request, "datasets", protocol.SearchDatasetsResponse)
