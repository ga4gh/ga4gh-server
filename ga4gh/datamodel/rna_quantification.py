"""
Module responsible for translating feature expression data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.datamodel as datamodel
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
import ga4gh.sqliteBackend as sqliteBackend
import ga4gh.pb as pb


"""
    The RNA Quantifications associated with a GA4GH dataset reside in a sqlite
    database which is contained in the rnaQuant subdirectory of the dataset
    directory.

    The sqlite .db file has 2 tables:
    RNAQUANTIFICATION : contains rnaQuantification data
    EXPRESSION : contains feature level expression data

    Desired GA4GH objects will be generated on the fly by the dictionaries
    returned by database queries and sent to the backend.
"""


class AbstractExpressionLevel(datamodel.DatamodelObject):
    """
    An abstract base class of a expression level
    """
    compoundIdClass = datamodel.ExpressionLevelCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractExpressionLevel, self).__init__(
            parentContainer, localId)
        self._expression = 0.0
        self._quantificationGroupId = ""
        self._isNormalized = ""
        self._rawReadCount = 0.0
        self._score = 0.0
        self._units = ""
        self._name = localId
        self._confInterval = []

    def toProtocolElement(self):
        protocolElement = protocol.ExpressionLevel()
        protocolElement.id = self.getId()
        protocolElement.name = self._name
        protocolElement.quantification_group_id = (self._quantificationGroupId)
        protocolElement.raw_read_count = self._rawReadCount
        protocolElement.expression = self._expression
        protocolElement.is_normalized = self._isNormalized
        protocolElement.units = self._units
        protocolElement.score = self._score
        protocolElement.conf_interval = self._confInterval
        return protocolElement


class ExpressionLevel(AbstractExpressionLevel):
    """
    Class representing a single ExpressionLevel in the GA4GH data model.
    """

    def __init__(self, parentContainer, record):
        super(ExpressionLevel, self).__init__(parentContainer, record["id"])
        self._expression = record["expression"]
        self._quantificationGroupId = record["quantification_group_id"]
        # sqlite stores booleans as int (False = 0, True = 1)
        self._isNormalized = bool(record["is_normalized"])
        self._rawReadCount = record["raw_read_count"]
        self._score = record["score"]
        self._units = record["units"]
        self._name = record["name"]
        self._confInterval = [record["conf_low"], record["conf_hi"]]

    def getName(self):
        return self._name

    def getQuantificationGroup(self):
        return self._quantificationGroupId


class AbstractQuantificationGroup(datamodel.DatamodelObject):
    """
    Class representing a single QuantificationGroup in the GA4GH model.
    """
    compoundIdClass = datamodel.QuantificationGroupCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractQuantificationGroup, self).__init__(
            parentContainer, localId)
        self._analysisId = ""
        self._name = localId
        self._description = ""
        self._info = []

    def toProtocolElement(self):
        protocolElement = protocol.QuantificationGroup()
        protocolElement.id = self.getId()
        protocolElement.analysis_id = self._analysisId
        protocolElement.name = self._name
        protocolElement.description = self._description
        protocolElement.info = self._info
        return protocolElement


class QuantificationGroup(AbstractQuantificationGroup):
    """
    Class representing a single QuantificationGroup in the GA4GH model.
    """

    def __init__(self, parentContainer, record):
        super(QuantificationGroup, self).__init__(
            parentContainer, record["quantification_group_id"])
        self._analysisId = record["rna_quantification_id"]
        self._name = record["quantification_group_id"]


class AbstractRNAQuantification(datamodel.DatamodelObject):
    """
    An abstract base class of a RNA quantification
    """
    compoundIdClass = datamodel.RnaQuantificationCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractRNAQuantification, self).__init__(
            parentContainer, localId)
        self._annotationIds = []
        self._description = ""
        self._name = localId
        self._readGroupId = ""
        self._referenceSet = None
        self._programIds = []

    def toProtocolElement(self):
        """
        Converts this rnaQuant into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.RnaQuantification()
        protocolElement.id = self.getId()
        protocolElement.name = self._name
        protocolElement.description = self._description
        protocolElement.read_group_id = self._readGroupId
        protocolElement.program_ids = self._programIds
        protocolElement.annotation_ids = self._annotationIds
        return protocolElement

    def addRnaQuantMetadata(self, fields):
        """
        data elements are:
        Id, annotations, description, name, readGroupId
        where annotations is a comma separated list
        """
        self._annotationIds = fields["annotation_ids"].split(',')
        self._description = fields["description"]
        self._name = fields["name"]
        self._readGroupId = fields["read_group_id"]

    def getReferenceSet(self):
        """
        Returns the reference set associated with this RnaQuantification.
        """
        return self._referenceSet

    def setReferenceSet(self, referenceSet):
        """
        Sets the reference set associated with this RnaQuantification to the
        specified value.
        """
        self._referenceSet = referenceSet


class RNASeqResult(AbstractRNAQuantification):
    """
    Class representing a single RnaQuantification in the GA4GH data model.
    """

    def __init__(self, parentContainer, localId, rnaQuantDataPath=None):
        super(RNASeqResult, self).__init__(parentContainer, localId)
        self._dbFilePath = None
        self._db = None

    def getRnaQuantMetadata(self):
        """
        input is tab file with no header.  Columns are:
        Id, annotations, description, name, readGroupId
        where annotation is a comma separated list
        """
        rnaQuantId = self.getLocalId()
        with self._db as dataSource:
            rnaQuantReturned = dataSource.getRnaQuantificationById(
                rnaQuantId)
        self.addRnaQuantMetadata(rnaQuantReturned)

    def populateFromFile(self, dataUrl):
        """
        Populates the instance variables of this FeatureSet from the specified
        data URL.
        """
        self._dbFilePath = dataUrl
        self._db = SqliteRNABackend(self._dbFilePath)
        self.getRnaQuantMetadata()

    def populateFromRow(self, row):
        """
        Populates the instance variables of this FeatureSet from the specified
        DB row.
        """
        self._dbFilePath = row[b'dataUrl']
        self._db = SqliteRNABackend(self._dbFilePath)
        self.getRnaQuantMetadata()

    def getDataUrl(self):
        """
        Returns the URL providing the data source for this FeatureSet.
        """
        return self._dbFilePath

    def getExpressionLevels(
            self, rnaQuantID, pageToken=0, pageSize=None, expressionId=None,
            quantificationGroupId=None, threshold=0.0):
        """
        Returns the list of ExpressionLevels in this RNA Quantification.
        """
        if expressionId is not None:
            parsedId = datamodel.ExpressionLevelCompoundId.parse(expressionId)
            expressionId = parsedId.expressionLevelId
        if quantificationGroupId is not None:
            parsedId = datamodel.QuantificationGroupCompoundId.parse(
                quantificationGroupId)
            quantificationGroupId = parsedId.quantificationGroupId
        with self._db as dataSource:
            expressionsReturned = dataSource.searchExpressionLevelsInDb(
                rnaQuantID, pageToken=pageToken,
                pageSize=pageSize, expressionId=expressionId,
                quantificationGroupId=quantificationGroupId,
                threshold=threshold)
        return [ExpressionLevel(self, expressionEntry) for
                expressionEntry in expressionsReturned]

    def getExpressionLevel(self, compoundId):
        expressionId = compoundId.expressionLevelId
        with self._db as dataSource:
            expressionReturned = dataSource.getExpressionLevelById(
                expressionId)

        return ExpressionLevel(self, expressionReturned)

    def getQuantificationGroups(
            self, rnaQuantID, pageToken=0, pageSize=None,
            quantificationGroupId=None):
        """
        Returns the list of QuantificationGroups in this RNA Quantification.

        for now the quantification group data is autogenerated by examining the
        relevant expression data file
        """
        if quantificationGroupId is not None:
            parsedId = datamodel.QuantificationGroupCompoundId.parse(
                quantificationGroupId)
            quantificationGroupId = parsedId.quantificationGroupId
        with self._db as dataSource:
            quantificationGroupsReturned = \
                dataSource.searchQuantificationGroupsInDb(
                    rnaQuantID, pageToken=pageToken,
                    pageSize=pageSize,
                    quantificationGroupId=quantificationGroupId)
        return [QuantificationGroup(self, quantificationGroupEntry) for
                quantificationGroupEntry in quantificationGroupsReturned]

    def getQuantificationGroup(self, compoundId):
        """
        for now the quantification group data is autogenerated by examining the
        relevant expression data file
        """
        quantificationGroupId = compoundId.quantificationGroupId
        with self._db as dataSource:
            quantificationGroupReturned = \
                dataSource.getQuantificationGroupById(quantificationGroupId)

        return QuantificationGroup(self, quantificationGroupReturned)


class SqliteRNABackend(sqliteBackend.SqliteBackedDataSource):
    """
    Defines an interface to a sqlite DB which stores all RNA quantifications
    in the dataset.
    """
    def __init__(self, rnaQuantSqlFile="ga4gh-rnaQuant.db"):
        super(SqliteRNABackend, self).__init__(rnaQuantSqlFile)

    def searchRnaQuantificationsInDb(
            self, pageToken=0, pageSize=None, rnaQuantificationId=None):
        """
        :param pageToken: int representing first record to return
        :param pageSize: int representing number of records to return
        :param rnaQuantificationId: string restrict search by id
        :return an array of dictionaries, representing the returned data.
        """
        sql = ("SELECT * FROM RNAQUANTIFICATION")
        sql_args = ()
        if rnaQuantificationId is not None:
            sql += " WHERE id = ? "
            sql_args += (rnaQuantificationId,)
        sql += sqliteBackend.limitsSql(pageToken, pageSize)
        query = self._dbconn.execute(sql, sql_args)
        try:
            return sqliteBackend.sqliteRows2dicts(query.fetchall())
        except AttributeError:
            raise exceptions.RnaQuantificationNotFoundException(
                rnaQuantificationId)

    def getRnaQuantificationById(self, rnaQuantificationId):
        """
        :param rnaQuantificationId: the RNA Quantification ID
        :return: dictionary representing an RnaQuantification object,
            or None if no match is found.
        """
        sql = ("SELECT * FROM RNAQUANTIFICATION WHERE id = ?")
        query = self._dbconn.execute(sql, (rnaQuantificationId,))
        try:
            return sqliteBackend.sqliteRow2Dict(query.fetchone())
        except AttributeError:
            raise exceptions.RnaQuantificationNotFoundException(
                rnaQuantificationId)

    def searchExpressionLevelsInDb(
            self, rnaQuantId, pageToken=0, pageSize=None, expressionId=None,
            quantificationGroupId=None, threshold=0.0):
        """
        :param rnaQuantId: string restrict search by quantification id
        :param pageToken: int representing first record to return
        :param pageSize: int representing number of records to return
        :param expressionId: string restrict search by expression id
        :param quantificationGroupId: string restrict search by quantification
            group id
        :param threshold: float minimum expression values to return
        :return an array of dictionaries, representing the returned data.
        """
        sql = ("SELECT * FROM EXPRESSION WHERE "
               "rna_quantification_id = ? "
               "AND expression >= ? ")
        sql_args = (rnaQuantId, threshold)
        if expressionId is not None:
            sql += "AND id = ? "
            sql_args += (expressionId,)
        if quantificationGroupId is not None:
            sql += "AND quantification_group_id = ? "
            sql_args += (quantificationGroupId,)
        sql += sqliteBackend.limitsSql(pageToken, pageSize)
        query = self._dbconn.execute(sql, sql_args)
        try:
            return sqliteBackend.sqliteRows2dicts(query.fetchall())
        except AttributeError:
            raise exceptions.ExpressionLevelNotFoundException(
                expressionId)

    def getExpressionLevelById(self, expressionId):
        """
        :param expressionId: the ExpressionLevel ID
        :return: dictionary representing an ExpressionLevel object,
            or None if no match is found.
        """
        sql = ("SELECT * FROM EXPRESSION WHERE id = ?")
        query = self._dbconn.execute(sql, (expressionId,))
        try:
            return sqliteBackend.sqliteRow2Dict(query.fetchone())
        except AttributeError:
            raise exceptions.ExpressionLevelNotFoundException(
                expressionId)

    def searchQuantificationGroupsInDb(
            self, rnaQuantId, pageToken=0, pageSize=None,
            quantificationGroupId=None):
        """
        :param rnaQuantId: string restrict search by quantification id
        :param pageToken: int representing first record to return
        :param pageSize: int representing number of records to return
        :param quantificationGroupId: string restrict search by quantification
            group id
        :return an array of dictionaries, representing the returned data.
        """
        sql = ("SELECT DISTINCT quantification_group_id,rna_quantification_id "
               "FROM EXPRESSION WHERE rna_quantification_id = ? ")
        sql_args = (rnaQuantId,)
        if quantificationGroupId is not None:
            sql += "AND quantification_group_id = ? "
            sql_args += (quantificationGroupId,)
        sql += sqliteBackend.limitsSql(pageToken, pageSize)
        query = self._dbconn.execute(sql, sql_args)
        try:
            return sqliteBackend.sqliteRows2dicts(query.fetchall())
        except AttributeError:
            raise exceptions.QuantificationGroupNotFoundException(
                quantificationGroupId)

    def getQuantificationGroupById(self, quantificationGroupId):
        """
        :param quantificationGroupId: the QuantificationGroup ID
        :return: dictionary representing a QuantificationGroup object,
            or None if no match is found.
        """
        sql = ("SELECT quantification_group_id,rna_quantification_id FROM "
               "EXPRESSION WHERE quantification_group_id = ?")
        query = self._dbconn.execute(sql, (quantificationGroupId,))
        try:
            return sqliteBackend.sqliteRow2Dict(query.fetchone())
        except AttributeError:
            raise exceptions.QuantificationGroupNotFoundException(
                quantificationGroupId)


class SimulatedRNASeqResult(AbstractRNAQuantification):
    """
    An RNA Quantification that doesn't derive from a data store.
    Used mostly for testing.
    """
    # TODO: this needs to be updated/rewritten
    def __init__(self, parentContainer, localId, rnaQuantDataPath=""):
        super(SimulatedRNASeqResult, self).__init__(parentContainer, localId)
        self._rnaQuantDataPath = rnaQuantDataPath
        self.generateRnaQuantMetadata()

    def generateRnaQuantMetadata(self):
        """
            Currently just sets to default values.
        """
        fields = {"annotation_ids": "",
                  "description": "",
                  "name": "",
                  "read_group_id": ""}
        self.addRnaQuantMetadata(fields)
