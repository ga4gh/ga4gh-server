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
        self._featureId = ""
        self._isNormalized = ""
        self._rawReadCount = 0.0
        self._score = 0.0
        self._units = 0
        self._name = localId
        self._confIntervalLow = 0.0
        self._confIntervalHigh = 0.0
        self._featureGroupIds = []

    def toProtocolElement(self):
        protocolElement = protocol.ExpressionLevel()
        protocolElement.id = self.getId()
        protocolElement.name = self._name
        protocolElement.feature_id = self._featureId
        protocolElement.rna_quantification_id = parentContainer.getId()
        protocolElement.raw_read_count = self._rawReadCount
        protocolElement.expression = self._expression
        protocolElement.is_normalized = self._isNormalized
        protocolElement.units = self._units
        protocolElement.score = self._score
        protocolElement.conf_interval_low = self._confIntervalLow
        protocolElement.conf_interval_high = self._confIntervalHigh
        protocolElement.feature_group_ids.extend(self._featureGroupIds)
        return protocolElement


class ExpressionLevel(AbstractExpressionLevel):
    """
    Class representing a single ExpressionLevel in the GA4GH data model.
    """

    def __init__(self, parentContainer, record):
        super(ExpressionLevel, self).__init__(parentContainer, record["id"])
        self._expression = record["expression"]
        self._featureId = record["feature_id"]
        # sqlite stores booleans as int (False = 0, True = 1)
        self._isNormalized = bool(record["is_normalized"])
        self._rawReadCount = record["raw_read_count"]
        self._score = record["score"]
        self._units = record["units"]
        self._name = record["name"]
        self._confIntervalLow = record["conf_low"]
        self._confIntervalHigh = record["conf_hi"]
        self._featureGroupIds = record["feature_group_ids"].split(",")

    def getName(self):
        return self._name

    def getFeatureGroups(self):
        return self._featureGroupIds


class AbstractFeatureGroup(datamodel.DatamodelObject):
    """
    Class representing a single FeatureGroup in the GA4GH model.
    """
    compoundIdClass = datamodel.FeatureGroupCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractFeatureGroup, self).__init__(
            parentContainer, localId)
        self._name = localId
        self._description = ""
        self._feature_ids = []

    def toProtocolElement(self):
        protocolElement = protocol.QuantificationGroup()
        protocolElement.id = self.getId()
        protocolElement.name = self._name
        protocolElement.description = self._description
        protocolElement.feature_ids.extend(self._feature_ids)

        return protocolElement


class FeatureGroup(AbstractFeatureGroup):
    """
    Class representing a single FeatureGroup in the GA4GH model.
    """

    def __init__(self, parentContainer, record):
        super(FeatureGroup, self).__init__(
            parentContainer, record["name"])
        self._featureIds = record["feature_ids"].split(",")
        self._description = record["description"]


class AbstractRNAQuantificationSet(datamodel.DatamodelObject):
    """
    An abstract base class of a RNA quantification set
    """
    compoundIdClass = datamodel.RnaQuantificationSetCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractRNAQuantificationSet, self).__init__(
            parentContainer, localId)
        self._name = localId
        self._referenceSet = None

    def toProtocolElement(self):
        """
        Converts this rnaQuant into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.RnaQuantificationSet()
        protocolElement.id = self.getId()
        protocolElement.dataset_id = parentContainer.getId()
        protocolElement.name = self._name


class RnaQuantificationSet(AbstractRNAQuantificationSet):
    """
    Class representing a single RnaQuantificationSet in the GA4GH model.
    """

    def __init__(self, parentContainer, name):
        super(RnaQuantificationSet, self).__init__(
            parentContainer, name)
        self._dbFilePath = None
        self._db = None

    def getReferenceSet(self):
        """
        Returns the reference set associated with this RnaQuantificationSet.
        """
        return self._referenceSet

    def setReferenceSet(self, referenceSet):
        """
        Sets the reference set associated with this RnaQuantificationSet to the
        specified value.
        """
        self._referenceSet = referenceSet

    def populateFromFile(self, dataUrl):
        """
        Populates the instance variables of this RnaQuantificationSet from the
        specified data URL.
        """
        self._dbFilePath = dataUrl
        self._db = SqliteRNABackend(self._dbFilePath)

    def populateFromRow(self, row):
        """
        Populates the instance variables of this RnaQuantificationSet from the
        specified DB row.
        """
        self._dbFilePath = row[b'dataUrl']
        self._db = SqliteRNABackend(self._dbFilePath)

    def getDataUrl(self):
        """
        Returns the URL providing the data source for this
        RnaQuantificationSet.
        """
        return self._dbFilePath


class AbstractRNAQuantification(datamodel.DatamodelObject):
    """
    An abstract base class of a RNA quantification
    """
    compoundIdClass = datamodel.RnaQuantificationCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractRNAQuantification, self).__init__(
            parentContainer, localId)
        self._featureSetIds = []
        self._description = ""
        self._name = localId
        self._readGroupIds = []
        self._referenceSet = None
        self._programs = []

    def toProtocolElement(self):
        """
        Converts this rnaQuant into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.RnaQuantification()
        protocolElement.id = self.getId()
        protocolElement.name = self._name
        protocolElement.description = self._description
        protocolElement.read_group_ids.extend(self._readGroupIds)
        protocolElement.programs.extend(self._programs)
        protocolElement.feature_set_ids.extend(self._featureSetIds)
        protocolElement.rna_quantification_set_id = parentContainer.getId()
        return protocolElement

    def addRnaQuantMetadata(self, fields):
        """
        data elements are:
        Id, annotations, description, name, readGroupId
        where annotations is a comma separated list
        """
        self._annotationIds = fields["feature_set_ids"].split(',')
        self._description = fields["description"]
        self._name = fields["name"]
        self._readGroupId = fields["read_group_ids"].split(',')
        self._programs = fields["programs"].split(',')

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
            self, rnaQuantID, pageToken=0, pageSize=None, expressionId="",
            quantificationGroupId="", threshold=0.0):
        """
        Returns the list of ExpressionLevels in this RNA Quantification.
        """
        if len(expressionId) > 0:
            parsedId = datamodel.ExpressionLevelCompoundId.parse(expressionId)
            expressionId = parsedId.expressionLevelId
        if len(quantificationGroupId) > 0:
            parsedId = datamodel.QuantificationGroupCompoundId.parse(
                quantificationGroupId)
            quantificationGroupId = parsedId.quantification_group_id
        with self._db as dataSource:
            expressionsReturned = dataSource.searchExpressionLevelsInDb(
                rnaQuantID, pageToken=pageToken,
                pageSize=pageSize, expressionId=expressionId,
                quantificationGroupId=quantificationGroupId,
                threshold=threshold)
        return [ExpressionLevel(self, expressionEntry) for
                expressionEntry in expressionsReturned]

    def getExpressionLevel(self, compoundId):
        expressionId = compoundId.expression_level_id
        with self._db as dataSource:
            expressionReturned = dataSource.getExpressionLevelById(
                expressionId)

        return ExpressionLevel(self, expressionReturned)

    def getQuantificationGroups(
            self, rnaQuantID, pageToken=0, pageSize=None,
            quantificationGroupId=""):
        """
        Returns the list of QuantificationGroups in this RNA Quantification.

        for now the quantification group data is autogenerated by examining the
        relevant expression data file
        """
        if len(quantificationGroupId) > 0:
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
        quantificationGroupId = compoundId.quantification_group_id
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
            self, pageToken=0, pageSize=None, rnaQuantificationId=""):
        """
        :param pageToken: int representing first record to return
        :param pageSize: int representing number of records to return
        :param rnaQuantificationId: string restrict search by id
        :return an array of dictionaries, representing the returned data.
        """
        sql = ("SELECT * FROM RNAQUANTIFICATION")
        sql_args = ()
        if len(rnaQuantificationId) > 0:
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
            self, rnaQuantId, pageToken=0, pageSize=None, expressionId="",
            quantificationGroupId="", threshold=0.0):
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
        if len(expressionId) > 0:
            sql += "AND id = ? "
            sql_args += (expressionId,)
        if len(quantificationGroupId) > 0:
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
            quantificationGroupId=""):
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
        if len(quantificationGroupId) > 0:
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
