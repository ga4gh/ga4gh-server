"""
Client classes for the GA4GH reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import requests
import posixpath
import logging

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


class AbstractClient(object):
    """
    The abstract superclass of GA4GH Client objects.
    """

    def __init__(self, logLevel=0):
        self._pageSize = None
        self._logLevel = logLevel
        self._protocolBytesReceived = 0
        logging.basicConfig()
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(logLevel)

    def _deserializeResponse(self, jsonResponseString, protocolResponseClass):
        self._protocolBytesReceived += len(jsonResponseString)
        self._logger.debug("response:{}".format(jsonResponseString))
        if jsonResponseString == '':
            raise exceptions.EmptyResponseException()
        responseObject = protocolResponseClass.fromJsonString(
            jsonResponseString)
        return responseObject

    def _runSearchPageRequest(
            self, protocolRequest, objectName, protocolResponseClass):
        """
        Runs a complete transaction with the server to obtain a single
        page of search results.
        """
        raise NotImplemented()

    def _runSearchRequest(
            self, protocolRequest, objectName, protocolResponseClass):
        """
        Runs the specified request at the specified objectName and instantiates
        an object of the specified class. We yield each object in listAttr.
        If pages of results are present, repeat this process until the
        pageToken is null.
        """
        notDone = True
        while notDone:
            responseObject = self._runSearchPageRequest(
                protocolRequest, objectName, protocolResponseClass)
            valueList = getattr(
                responseObject, protocolResponseClass.getValueListName())
            for extract in valueList:
                yield extract
            notDone = responseObject.nextPageToken is not None
            protocolRequest.pageToken = responseObject.nextPageToken

    def _runListReferenceBasesPageRequest(self, id_, protocolRequest):
        """
        Runs a complete transaction with the server to get a single
        page of results for the specified ListReferenceBasesRequest.
        """
        raise NotImplemented()

    def listReferenceBases(self, id_, start=0, end=None):
        """
        Returns an iterator over the bases from the server in the form
        of consecutive strings. This command does not conform to the
        patterns of the other search and get requests, and is implemented
        differently.
        """
        request = protocol.ListReferenceBasesRequest()
        request.start = start
        request.end = end
        notDone = True
        # TODO We should probably use a StringIO here to make string buffering
        # a bit more efficient.
        basesList = []
        while notDone:
            response = self._runListReferenceBasesPageRequest(id_, request)
            basesList.append(response.sequence)
            notDone = response.nextPageToken is not None
            request.pageToken = response.nextPageToken
        return "".join(basesList)

    def _runGetRequest(self, objectName, protocolResponseClass, id_):
        """
        Requests an object from the server and returns the object of
        type protocolResponseClass that has id id_.
        Used for requests where a single object is the expected response.
        """
        raise NotImplemented()

    def getPageSize(self):
        """
        Returns the suggested maximum size of pages of results returned by
        the server.
        """
        return self._pageSize

    def setPageSize(self, pageSize):
        """
        Sets the requested maximum size of pages of results returned by the
        server to the specified value.
        """
        self._pageSize = pageSize

    def getProtocolBytesReceived(self):
        """
        Returns the total number of protocol bytes received from the server
        by this client.

        :return: The number of bytes consumed by protocol traffic read from
            the server during the lifetime of this client.
        :rtype: int
        """
        return self._protocolBytesReceived

    def getDataset(self, datasetId):
        """
        Returns the Dataset with the specified ID from the server.

        :param str datasetId: The ID of the Dataset of interest.
        :return: The Dataset of interest.
        :rtype: :class:`ga4gh.protocol.Dataset`
        """
        return self._runGetRequest("datasets", protocol.Dataset, datasetId)

    def getReferenceSet(self, referenceSetId):
        """
        Returns the ReferenceSet with the specified ID from the server.

        :param str referenceSetId: The ID of the ReferenceSet of interest.
        :return: The ReferenceSet of interest.
        :rtype: :class:`ga4gh.protocol.ReferenceSet`
        """
        return self._runGetRequest(
            "referencesets", protocol.ReferenceSet, referenceSetId)

    def getReference(self, referenceId):
        """
        Returns the Reference with the specified ID from the server.

        :param str referenceId: The ID of the Reference of interest.
        :return: The Reference of interest.
        :rtype: :class:`ga4gh.protocol.Reference`
        """
        return self._runGetRequest(
            "references", protocol.Reference, referenceId)

    def getReadGroupSet(self, readGroupSetId):
        """
        Returns the ReadGroupSet with the specified ID from the server.

        :param str readGroupSetId: The ID of the ReadGroupSet of interest.
        :return: The ReadGroupSet of interest.
        :rtype: :class:`ga4gh.protocol.ReadGroupSet`
        """
        return self._runGetRequest(
            "readgroupsets", protocol.ReadGroupSet, readGroupSetId)

    def getReadGroup(self, readGroupId):
        """
        Returns the ReadGroup with the specified ID from the server.

        :param str readGroupId: The ID of the ReadGroup of interest.
        :return: The ReadGroup of interest.
        :rtype: :class:`ga4gh.protocol.ReadGroup`
        """
        return self._runGetRequest(
            "readgroups", protocol.ReadGroup, readGroupId)

    def getCallSet(self, callSetId):
        """
        Returns the CallSet with the specified ID from the server.

        :param str callSetId: The ID of the CallSet of interest.
        :return: The CallSet of interest.
        :rtype: :class:`ga4gh.protocol.CallSet`
        """
        return self._runGetRequest("callsets", protocol.CallSet, callSetId)

    def getVariant(self, variantId):
        """
        Returns the Variant with the specified ID from the server.

        :param str variantId: The ID of the Variant of interest.
        :return: The Variant of interest.
        :rtype: :class:`ga4gh.protocol.Variant`
        """
        return self._runGetRequest("variants", protocol.Variant, variantId)

    def getVariantSet(self, variantSetId):
        """
        Returns the VariantSet with the specified ID from the server.

        :param str variantSetId: The ID of the VariantSet of interest.
        :return: The VariantSet of interest.
        :rtype: :class:`ga4gh.protocol.VariantSet`
        """
        return self._runGetRequest(
            "variantsets", protocol.VariantSet, variantSetId)

    def searchVariants(
            self, variantSetId, start=None, end=None, referenceName=None,
            callSetIds=None):
        """
        Returns an iterator over the Variants fulfilling the specified
        conditions from the specified VariantSet.

        :param str variantSetId: The ID of the
            :class:`ga4gh.protocol.VariantSet` of interest.
        :param int start: Required. The beginning of the window (0-based,
            inclusive) for which overlapping variants should be returned.
            Genomic positions are non-negative integers less than reference
            length. Requests spanning the join of circular genomes are
            represented as two requests one on each side of the join
            (position 0).
        :param int end: Required. The end of the window (0-based, exclusive)
            for which overlapping variants should be returned.
        :param str referenceName: The name of the
            :class:`ga4gh.protocol.Reference` we wish to return variants from.
        :param list callSetIds: Only return variant calls which belong to call
            sets with these IDs. If an empty array, returns variants without
            any call objects. If null, returns all variant calls.

        :return: An iterator over the :class:`ga4gh.protocol.Variant` objects
            defined by the query parameters.
        :rtype: iter
        """
        request = protocol.SearchVariantsRequest()
        request.referenceName = referenceName
        request.start = start
        request.end = end
        request.variantSetId = variantSetId
        request.callSetIds = callSetIds
        request.pageSize = self._pageSize
        return self._runSearchRequest(
            request, "variants", protocol.SearchVariantsResponse)

    def searchDatasets(self):
        """
        Returns an iterator over the Datasets on the server.

        :return: An iterator over the :class:`ga4gh.protocol.Dataset`
            objects on the server.
        """
        request = protocol.SearchDatasetsRequest()
        request.pageSize = self._pageSize
        return self._runSearchRequest(
            request, "datasets", protocol.SearchDatasetsResponse)

    def searchVariantSets(self, datasetId):
        """
        Returns an iterator over the VariantSets fulfilling the specified
        conditions from the specified Dataset.

        :param str datasetId: The ID of the :class:`ga4gh.protocol.Dataset`
            of interest.
        :return: An iterator over the :class:`ga4gh.protocol.VariantSet`
            objects defined by the query parameters.
        """
        request = protocol.SearchVariantSetsRequest()
        request.datasetId = datasetId
        request.pageSize = self._pageSize
        return self._runSearchRequest(
            request, "variantsets", protocol.SearchVariantSetsResponse)

    def searchReferenceSets(
            self, accession=None, md5checksum=None, assemblyId=None):
        """
        Returns an iterator over the ReferenceSets fulfilling the specified
        conditions.

        :param str accession: If not null, return the reference sets for which
            the `accession` matches this string (case-sensitive, exact match).
        :param str md5checksum: If not null, return the reference sets for
            which the `md5checksum` matches this string (case-sensitive, exact
            match). See :class:`ga4gh.protocol.ReferenceSet::md5checksum` for
            details.
        :param str assemblyId: If not null, return the reference sets for which
            the `assemblyId` matches this string (case-sensitive, exact match).
        :return: An iterator over the :class:`ga4gh.protocol.ReferenceSet`
            objects defined by the query parameters.
        """
        request = protocol.SearchReferenceSetsRequest()
        request.accession = accession
        request.md5checksum = md5checksum
        request.assemblyId = assemblyId
        request.pageSize = self._pageSize
        return self._runSearchRequest(
            request, "referencesets", protocol.SearchReferenceSetsResponse)

    def searchReferences(
            self, referenceSetId, accession=None, md5checksum=None):
        """
        Returns an iterator over the References fulfilling the specified
        conditions from the specified Dataset.

        :param str referenceSetId: The ReferenceSet to search.
        :param str accession: If not None, return the references for which the
            `accession` matches this string (case-sensitive, exact match).
        :param str md5checksum: If not None, return the references for which
            the `md5checksum` matches this string (case-sensitive, exact
            match).
        :return: An iterator over the :class:`ga4gh.protocol.Reference`
            objects defined by the query parameters.
        """
        request = protocol.SearchReferencesRequest()
        request.referenceSetId = referenceSetId
        request.accession = accession
        request.md5checksum = md5checksum
        request.pageSize = self._pageSize
        return self._runSearchRequest(
            request, "references", protocol.SearchReferencesResponse)

    def searchCallSets(self, variantSetId, name=None):
        """
        Returns an iterator over the CallSets fulfilling the specified
        conditions from the specified VariantSet.

        :param str name: Only CallSets matching the specified name will
            be returned.
        :return: An iterator over the :class:`ga4gh.protocol.CallSet`
            objects defined by the query parameters.
        """
        request = protocol.SearchCallSetsRequest()
        request.variantSetId = variantSetId
        request.name = name
        request.pageSize = self._pageSize
        return self._runSearchRequest(
            request, "callsets", protocol.SearchCallSetsResponse)

    def searchReadGroupSets(self, datasetId, name=None):
        """
        Returns an iterator over the ReadGroupSets fulfilling the specified
        conditions from the specified Dataset.

        :param str name: Only ReadGroupSets matching the specified name
            will be returned.
        :return: An iterator over the :class:`ga4gh.protocol.ReadGroupSet`
            objects defined by the query parameters.
        :rtype: iter
        """
        request = protocol.SearchReadGroupSetsRequest()
        request.datasetId = datasetId
        request.name = name
        request.pageSize = self._pageSize
        return self._runSearchRequest(
            request, "readgroupsets", protocol.SearchReadGroupSetsResponse)

    def searchReads(
            self, readGroupIds, referenceId=None, start=None, end=None):
        """
        Returns an iterator over the Reads fulfilling the specified
        conditions from the specified ReadGroupIds.

        :param str readGroupIds: The IDs of the
            :class:`ga4gh.protocol.ReadGroup` of interest.
        :param str referenceId: The name of the
            :class:`ga4gh.protocol.Reference` we wish to return reads
            mapped to.
        :param int start: The start position (0-based) of this query. If a
            reference is specified, this defaults to 0. Genomic positions are
            non-negative integers less than reference length. Requests spanning
            the join of circular genomes are represented as two requests one on
            each side of the join (position 0).
        :param int end: The end position (0-based, exclusive) of this query.
            If a reference is specified, this defaults to the reference's
            length.
        :return: An iterator over the
            :class:`ga4gh.protocol.ReadAlignment` objects defined by
            the query parameters.
        :rtype: iter
        """
        request = protocol.SearchReadsRequest()
        request.readGroupIds = readGroupIds
        request.referenceId = referenceId
        request.start = start
        request.end = end
        request.pageSize = self._pageSize
        return self._runSearchRequest(
            request, "reads", protocol.SearchReadsResponse)


class HttpClient(AbstractClient):
    """
    The GA4GH HTTP client. This class provides methods corresponding to the
    GA4GH search and object GET methods.

    .. todo:: Add a better description of the role of this class and include
        links to the high-level API documention.

    :param str urlPrefix: The base URL of the GA4GH server we wish to
        communicate with. This should include the 'http' or 'https' prefix.
    :param int logLevel: The amount of debugging information to log using
        the :mod:`logging` module. This is :data:`logging.WARNING` by default.
    :param str authenticationKey: The authentication key provided by the
        server after logging in.
    """

    def __init__(
            self, urlPrefix, logLevel=logging.WARNING, authenticationKey=None):
        super(HttpClient, self).__init__(logLevel)
        self._urlPrefix = urlPrefix
        self._authenticationKey = authenticationKey
        self._session = requests.Session()
        self._setupHttpSession()
        requestsLog = logging.getLogger("requests.packages.urllib3")
        requestsLog.setLevel(logLevel)
        requestsLog.propagate = True

    def _setupHttpSession(self):
        """
        Sets up the common HTTP session parameters used by requests.
        """
        headers = {"Content-type": "application/json"}
        self._session.headers.update(headers)
        # TODO is this unsafe????
        self._session.verify = False

    def _checkResponseStatus(self, response):
        """
        Checks the speficied HTTP response from the requests package and
        raises an exception if a non-200 HTTP code was returned by the
        server.
        """
        if response.status_code != requests.codes.ok:
            self._logger.error("%s %s", response.status_code, response.text)
            raise exceptions.RequestNonSuccessException(
                "Url {0} had status_code {1}".format(
                    response.url, response.status_code))

    def _getHttpParameters(self):
        """
        Returns the basic HTTP parameters we need all requests.
        """
        return {'key': self._authenticationKey}

    def _runSearchPageRequest(
            self, protocolRequest, objectName, protocolResponseClass):
        url = posixpath.join(self._urlPrefix, objectName + '/search')
        data = protocolRequest.toJsonString()
        self._logger.debug("request:{}".format(data))
        response = self._session.post(
            url, params=self._getHttpParameters(), data=data)
        self._checkResponseStatus(response)
        return self._deserializeResponse(response.text, protocolResponseClass)

    def _runGetRequest(self, objectName, protocolResponseClass, id_):
        urlSuffix = "{objectName}/{id}".format(objectName=objectName, id=id_)
        url = posixpath.join(self._urlPrefix, urlSuffix)
        response = self._session.get(url, params=self._getHttpParameters())
        self._checkResponseStatus(response)
        return self._deserializeResponse(response.text, protocolResponseClass)

    def _runListReferenceBasesPageRequest(self, id_, request):
        urlSuffix = "references/{id}/bases".format(id=id_)
        url = posixpath.join(self._urlPrefix, urlSuffix)
        params = self._getHttpParameters()
        params.update(request.toJsonDict())
        response = self._session.get(url, params=params)
        self._checkResponseStatus(response)
        return self._deserializeResponse(
            response.text, protocol.ListReferenceBasesResponse)


class LocalClient(AbstractClient):

    def __init__(self, backend):
        super(LocalClient, self).__init__()
        self._backend = backend
        self._getMethodMap = {
            "callsets": self._backend.runGetCallSet,
            "datasets": self._backend.runGetDataset,
            "referencesets": self._backend.runGetReferenceSet,
            "references": self._backend.runGetReference,
            "variantsets": self._backend.runGetVariantSet,
            "variants": self._backend.runGetVariant,
            "readgroupsets": self._backend.runGetReadGroupSet,
            "readgroups": self._backend.runGetReadGroup,
        }
        self._searchMethodMap = {
            "callsets": self._backend.runSearchCallSets,
            "datasets": self._backend.runSearchDatasets,
            "referencesets": self._backend.runSearchReferenceSets,
            "references": self._backend.runSearchReferences,
            "variantsets": self._backend.runSearchVariantSets,
            "variants": self._backend.runSearchVariants,
            "readgroupsets": self._backend.runSearchReadGroupSets,
            "reads": self._backend.runSearchReads,
        }

    def _runGetRequest(self, objectName, protocolResponseClass, id_):
        getMethod = self._getMethodMap[objectName]
        responseJson = getMethod(id_)
        return self._deserializeResponse(responseJson, protocolResponseClass)

    def _runSearchPageRequest(
            self, protocolRequest, objectName, protocolResponseClass):
        searchMethod = self._searchMethodMap[objectName]
        responseJson = searchMethod(protocolRequest.toJsonString())
        return self._deserializeResponse(responseJson, protocolResponseClass)

    def _runListReferenceBasesPageRequest(self, id_, request):
        requestArgs = request.toJsonDict()
        # We need to remove end from this dict if it's not specified because
        # of the way we're interacting with Flask and HTTP GET params.
        # TODO: This is a really nasty way of doing things; we really
        # should just have a request object and pass that around instead of an
        # arguments dictionary.
        if request.end is None:
            del requestArgs["end"]
        if request.pageToken is None:
            del requestArgs["pageToken"]
        responseJson = self._backend.runListReferenceBases(id_, requestArgs)
        return self._deserializeResponse(
            responseJson, protocol.ListReferenceBasesResponse)
