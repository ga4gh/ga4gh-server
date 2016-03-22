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

# TODO: update
"""
    The RNA Quantification data model:
    Top level is the RNAQuantification which will sit right under the dataset
    Next level sits under each RNAQuant and hold expressionLevel, readCount,
    Characterization

    So to get all the RNAQuants will need to do a search over all datasets for
    the ID/name of the desired Quant

    To get all expressionLevels will have to search over all RNAQuants for all
    Datasets.

    ----

    sqlite file is going to sit directly under the dataset.  It will have
    entries for all of the quants in the dataset.  Desired objects will be
    generated on the fly by the queries and returned to the backend.  After
    that will be treated the same as now w.r.t. toProtocolElement, etc.
"""


class ExpressionLevel(datamodel.DatamodelObject):
    """
    Class representing a single ExpressionLevel in the GA4GH data model.
    """
    compoundIdClass = datamodel.ExpressionLevelCompoundId

    def __init__(self, parentContainer, record):
        # TODO this is still using the human readableId as the _id value
        # (localID)
        super(ExpressionLevel, self).__init__(parentContainer, record["id"])
        self._id = record["id"]
        self._annotationId = record["annotation_id"]
        self._expression = record["expression"]
        self._featureGroupId = record["feature_group_id"]
        self._isNormalized = record["is_normalized"]
        self._rawReadCount = record["raw_read_count"]
        self._score = record["score"]
        self._units = record["units"]

    def getName(self):
        return self._id

    def getFeatureGroup(self):
        return self._featureGroupId

    def toProtocolElement(self):
        protocolElement = protocol.ExpressionLevel()
        protocolElement.annotationId = self._annotationId
        protocolElement.expression = self._expression
        protocolElement.featureGroupId = (self._featureGroupId)
        protocolElement.id = self._id
        protocolElement.isNormalized = self._isNormalized
        protocolElement.rawReadCount = self._rawReadCount
        protocolElement.score = self._score
        protocolElement.units = self._units
        return protocolElement


class FeatureGroup(datamodel.DatamodelObject):
    """
    Class representing a single FeatureGroup in the GA4GH model.
    """
    compoundIdClass = datamodel.FeatureGroupCompoundId

    # TODO: this is just a first pass stub to get working
    # - need to formalize input data
    def __init__(self, parentContainer, record):
        super(FeatureGroup, self).__init__(parentContainer,
                                           record["feature_group_id"])
        self._id = record["feature_group_id"]
        self._analysisId = record["rna_quantification_id"]
        self.name = record["feature_group_id"]

    def toProtocolElement(self):
        protocolElement = protocol.FeatureGroup()
        protocolElement.id = self._id
        protocolElement.analysisId = self._analysisId
        protocolElement.name = self._id
        return protocolElement


class AbstractRNAQuantification(datamodel.DatamodelObject):
    """
    An abstract base class of a RNA quantification
    """
    compoundIdClass = datamodel.RnaQuantificationCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractRNAQuantification, self).__init__(parentContainer,
                                                        localId)
        self._annotationIds = []
        self._description = ""
        self._name = localId
        self._readGroupId = ""

    def toProtocolElement(self):
        """
        Converts this rnaQuant into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.RnaQuantification()
        protocolElement.annotationIds = self._annotationIds
        protocolElement.description = self._description
        protocolElement.id = self.getId()
        protocolElement.name = self._name
        protocolElement.readGroupId = self._readGroupId
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


class RNASeqResult(AbstractRNAQuantification):
    """
    Class representing a single RnaQuantification in the GA4GH data model.
    """

    def __init__(self, parentContainer, localId, rnaQuantDataPath,
                 dataRepository):
        super(RNASeqResult, self).__init__(parentContainer, localId)
        self._dbFilePath = rnaQuantDataPath  # the full path of the db file
        self._dataRepository = dataRepository
        self._db = SqliteRNABackend(self._dbFilePath)
        self.getRnaQuantMetadata()

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
        if rnaQuantReturned is not None:
            self.addRnaQuantMetadata(rnaQuantReturned)
        else:
            raise exceptions.RnaQuantificationNotFoundException(rnaQuantId)

    def getExpressionLevels(self, rnaQuantID, pageToken=0, pageSize=None,
                            expressionId=None, featureGroupId=None,
                            threshold=0.0):
        """
        Returns the list of ExpressionLevels in this RNA Quantification.
        """
        with self._db as dataSource:
            expressionsReturned = dataSource.searchExpressionLevelsInDb(
                rnaQuantID, pageToken=pageToken,
                pageSize=pageSize, expressionId=expressionId,
                featureGroupId=featureGroupId, threshold=threshold)
        return [ExpressionLevel(self, expressionEntry) for
                expressionEntry in expressionsReturned]

    def getExpressionLevel(self, compoundId):
        expressionId = compoundId.expressionLevelId
        with self._db as dataSource:
            expressionReturned = dataSource.getExpressionLevelById(
                expressionId)

        if expressionReturned is not None:
            return ExpressionLevel(self, expressionReturned)
        else:
            raise exceptions.ExpressionLevelNotFoundException(compoundId)

    def getFeatureGroups(self, rnaQuantID, pageToken=0, pageSize=None,
                         featureGroupId=None):
        """
        Returns the list of FeatureGroups in this RNA Quantification.

        for now the feature group data is autogenerated by examining the
        relevant expression data file
        """
        with self._db as dataSource:
            expressionsReturned = dataSource.searchFeatureGroupsInDb(
                rnaQuantID, pageToken=pageToken,
                pageSize=pageSize, featureGroupId=featureGroupId)
        return [FeatureGroup(self, expressionEntry) for
                expressionEntry in expressionsReturned]

    def getFeatureGroup(self, compoundId):
        """
        for now the feature group data is autogenerated by examining the
        relevant expression data file
        """
        featureGroupId = compoundId.featureGroupId
        with self._db as dataSource:
            expressionReturned = dataSource.getFeatureGroupById(featureGroupId)

        if expressionReturned is not None:
            return FeatureGroup(self, expressionReturned)
        else:
            raise exceptions.FeatureGroupNotFoundException(compoundId)


_rnaQuantColumns = [('id', 'TEXT'),
                    ('annotation_ids', 'TEXT'),
                    ('description', 'TEXT'),
                    ('name', 'TEXT'),
                    ('read_group_id', 'TEXT')]


_expressionColumns = [('id', 'TEXT'),
                      ('name', 'TEXT'),
                      ('annotation_id', 'TEXT'),
                      ('expression', 'REAL'),
                      ('feature_group_id', 'TEXT'),
                      ('is_normalized', 'BOOLEAN'),
                      ('raw_read_count', 'REAL'),
                      ('score', 'REAL'),
                      ('units', 'TEXT')]


class SqliteRNABackend(sqliteBackend.SqliteBackedDataSource):
    """
    Defines an interface to a sqlite DB which stores all RNA quantifications
    in the dataset.
    """
    def __init__(self, rnaQuantSqlFile="ga4gh-rnaQuant.db"):
        super(SqliteRNABackend, self).__init__(rnaQuantSqlFile)
        self.rnaQuantColumnNames = [f[0] for f in _rnaQuantColumns]
        self.rnaQuantColumnTypes = [f[1] for f in _rnaQuantColumns]
        self.expressionColumnNames = [f[0] for f in _expressionColumns]
        self.expressionColumnTypes = [f[1] for f in _expressionColumns]

    def searchRnaQuantificationsInDb(self, pageToken=0, pageSize=None,
                                     rnaQuantificationId=None):
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
        return sqliteBackend.sqliteRows2dicts(query.fetchall())

    def getRnaQuantificationById(self, rnaQuantificationId):
        """
        :param rnaQuantificationId: the RNA Quantification ID
        :return: dictionary representing an RnaQuantification object,
            or None if no match is found.
        """
        sql = ("SELECT * FROM RNAQUANTIFICATION WHERE id = ?")
        query = self._dbconn.execute(sql, (rnaQuantificationId,))
        return sqliteBackend.sqliteRow2Dict(query.fetchone())

    def searchExpressionLevelsInDb(self, rnaQuantId, pageToken=0,
                                   pageSize=None, expressionId=None,
                                   featureGroupId=None, threshold=0.0):
        """
        :param rnaQuantId: string restrict search by quantification id
        :param pageToken: int representing first record to return
        :param pageSize: int representing number of records to return
        :param expressionId: string restrict search by expression id
        :param featureGroupId: string restrict search by feature group id
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
        if featureGroupId is not None:
            sql += "AND feature_group_id = ? "
            sql_args += (featureGroupId,)
        sql += sqliteBackend.limitsSql(pageToken, pageSize)
        query = self._dbconn.execute(sql, sql_args)
        return sqliteBackend.sqliteRows2dicts(query.fetchall())

    def getExpressionLevelById(self, expressionId):
        """
        :param expressionId: the ExpressionLevel ID
        :return: dictionary representing an ExpressionLevel object,
            or None if no match is found.
        """
        sql = ("SELECT * FROM EXPRESSION WHERE id = ?")
        query = self._dbconn.execute(sql, (expressionId,))
        return sqliteBackend.sqliteRow2Dict(query.fetchone())

    def searchFeatureGroupsInDb(self, rnaQuantId, pageToken=0,
                                pageSize=None, featureGroupId=None):
        """
        :param rnaQuantId: string restrict search by quantification id
        :param pageToken: int representing first record to return
        :param pageSize: int representing number of records to return
        :param featureGroupId: string restrict search by feature group id
        :return an array of dictionaries, representing the returned data.
        """
        sql = ("SELECT DISTINCT feature_group_id,rna_quantification_id FROM "
               "EXPRESSION WHERE rna_quantification_id = ? ")
        sql_args = (rnaQuantId)
        if featureGroupId is not None:
            sql += "AND feature_group_id = ? "
            sql_args += (featureGroupId,)
        sql += sqliteBackend.limitsSql(pageToken, pageSize)
        query = self._dbconn.execute(sql, sql_args)
        return sqliteBackend.sqliteRows2dicts(query.fetchall())

    def getFeatureGroupById(self, featureGroupId):
        """
        :param featureGroupId: the FeatureGroup ID
        :return: dictionary representing a FeatureGroup object,
            or None if no match is found.
        """
        sql = ("SELECT feature_group_id,rna_quantification_id FROM EXPRESSION "
               "WHERE feature_group_id = ?")
        query = self._dbconn.execute(sql, (featureGroupId,))
        return sqliteBackend.sqliteRow2Dict(query.fetchone())


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
