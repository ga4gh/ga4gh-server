"""
Module responsible for translating feature expression data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.server.datamodel as datamodel
import ga4gh.server.exceptions as exceptions
import ga4gh.server.sqlite_backend as sqlite_backend

import ga4gh.schemas.protocol as protocol


"""
    The RNA Quantifications associated with a GA4GH dataset reside in a sqlite
    database which is contained in the rnaQuant subdirectory of the dataset
    directory.

    The sqlite .db file has 2 tables:
    RnaQuantification : contains rnaQuantification data
    Expression : contains feature level expression data

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
        self._isNormalized = ""
        self._rawReadCount = 0.0
        self._score = 0.0
        self._units = 0
        self._name = localId
        self._confIntervalLow = 0.0
        self._confIntervalHigh = 0.0

    def toProtocolElement(self):
        protocolElement = protocol.ExpressionLevel()
        protocolElement.id = self.getId()
        protocolElement.name = self._name
        protocolElement.rna_quantification_id = self._parentContainer.getId()
        protocolElement.raw_read_count = self._rawReadCount
        protocolElement.expression = self._expression
        protocolElement.is_normalized = self._isNormalized
        protocolElement.units = self._units
        protocolElement.score = self._score
        protocolElement.conf_interval_low = self._confIntervalLow
        protocolElement.conf_interval_high = self._confIntervalHigh
        self.serializeAttributes(protocolElement)
        return protocolElement


class SqliteExpressionLevel(AbstractExpressionLevel):
    """
    Class representing a single ExpressionLevel in the GA4GH data model.
    """
    def __init__(self, parentContainer, record):
        super(SqliteExpressionLevel, self).__init__(
            parentContainer, str(record["id"]))
        self._expression = record["expression"]
        # sqlite stores booleans as int (False = 0, True = 1)
        self._isNormalized = bool(record["is_normalized"])
        self._rawReadCount = record["raw_read_count"]
        self._score = record["score"]
        self._units = record["units"]
        self._name = record["name"]
        self._confIntervalLow = record["conf_low"]
        self._confIntervalHigh = record["conf_hi"]

    def getName(self):
        return self._name


class AbstractRnaQuantificationSet(datamodel.DatamodelObject):
    """
    An abstract base class of a RNA quantification set
    """
    compoundIdClass = datamodel.RnaQuantificationSetCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractRnaQuantificationSet, self).__init__(
            parentContainer, localId)
        self._name = localId
        self._referenceSet = None
        self._rnaQuantificationIdMap = {}
        self._rnaQuantificationIds = []

    def getRnaQuantificationByIndex(self, index):
        """
        Returns the rna quantification at the specified index in this set.
        """
        return self._rnaQuantificationIdMap[
            self._rnaQuantificationIds[index]]

    def getRnaQuantification(self, rnaQuantificationId):
        try:
            return self._rnaQuantificationIdMap[rnaQuantificationId]
        except KeyError:
            raise exceptions.RnaQuantificationNotFoundException(
                rnaQuantificationId)

    def getRnaQuantifications(self):
        return [self._rnaQuantificationIdMap[id_] for
                id_ in self._rnaQuantificationIds]

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

    def addRnaQuantification(self, rnaQuantification):
        """
        Add an rnaQuantification to this rnaQuantificationSet
        """
        id_ = rnaQuantification.getId()
        self._rnaQuantificationIdMap[id_] = rnaQuantification
        self._rnaQuantificationIds.append(id_)

    def toProtocolElement(self):
        """
        Converts this rnaQuant into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.RnaQuantificationSet()
        protocolElement.id = self.getId()
        protocolElement.dataset_id = self._parentContainer.getId()
        protocolElement.name = self._name
        self.serializeAttributes(protocolElement)
        return protocolElement


class SqliteRnaQuantificationSet(AbstractRnaQuantificationSet):
    """
    Class representing a single RnaQuantificationSet in the GA4GH model.
    """
    def __init__(self, parentContainer, name):
        super(SqliteRnaQuantificationSet, self).__init__(
            parentContainer, name)
        self._dbFilePath = None
        self._db = None

    def getDataUrl(self):
        """
        Returns the URL providing the data source for this
        RnaQuantificationSet.
        """
        return self._dbFilePath

    def populateFromFile(self, dataUrl):
        """
        Populates the instance variables of this RnaQuantificationSet from the
        specified data URL.
        """
        self._dbFilePath = dataUrl
        self._db = SqliteRnaBackend(self._dbFilePath)
        self.addRnaQuants()

    def populateFromRow(self, quantificationSetRecord):
        """
        Populates the instance variables of this RnaQuantificationSet from the
        specified DB row.
        """
        self._dbFilePath = quantificationSetRecord.dataurl
        self.setAttributesJson(quantificationSetRecord.attributes)
        self._db = SqliteRnaBackend(self._dbFilePath)
        self.addRnaQuants()

    def addRnaQuants(self):
        with self._db as dataSource:
            rnaQuantsReturned = dataSource.searchRnaQuantificationsInDb()
            for rnaQuant in rnaQuantsReturned:
                rnaQuantification = SqliteRnaQuantification(
                    self, rnaQuant["name"])
                rnaQuantification.populateFromFile(self._dbFilePath)
                self.addRnaQuantification(rnaQuantification)


class AbstractRnaQuantification(datamodel.DatamodelObject):
    """
    An abstract base class of a RNA quantification
    """
    compoundIdClass = datamodel.RnaQuantificationCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractRnaQuantification, self).__init__(
            parentContainer, localId)
        self._featureSetIds = []
        self._description = ""
        self._name = localId
        self._readGroupIds = []
        self._referenceSet = None
        self._programs = []
        self._biosampleId = ""

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
        protocolElement.biosample_id = self._biosampleId
        protocolElement.feature_set_ids.extend(self._featureSetIds)
        protocolElement.rna_quantification_set_id = \
            self._parentContainer.getId()
        self.serializeAttributes(protocolElement)
        return protocolElement

    def addRnaQuantMetadata(self, fields):
        """
        data elements are:
        Id, annotations, description, name, readGroupId
        where annotations is a comma separated list
        """
        self._featureSetIds = fields["feature_set_ids"].split(',')
        self._description = fields["description"]
        self._name = fields["name"]
        self._biosampleId = fields.get("biosample_id", "")
        if fields["read_group_ids"] == "":
            self._readGroupIds = []
        else:
            self._readGroupIds = fields["read_group_ids"].split(',')
        if fields["programs"] == "":
            self._programs = []
        else:
            # Need to use program Id's here to generate a list of Programs
            # for now set to empty
            self._programs = []

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

    def setBiosampleId(self, biosampleId):
        """
        Associates this quantification with a sample.
        """
        self._biosampleId = biosampleId

    def getBiosampleId(self):
        return self._biosampleId


class SqliteRnaQuantification(AbstractRnaQuantification):
    """
    Class representing a single RnaQuantification in the GA4GH data model.
    """
    def __init__(self, parentContainer, localId):
        super(SqliteRnaQuantification, self).__init__(parentContainer, localId)
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
        self._db = SqliteRnaBackend(self._dbFilePath)
        self.getRnaQuantMetadata()

    def populateFromRow(self, row):
        """
        Populates the instance variables of this FeatureSet from the specified
        DB row.
        """
        self._dbFilePath = row[b'dataUrl']
        self._db = SqliteRnaBackend(self._dbFilePath)
        self.getRnaQuantMetadata()

    def getDataUrl(self):
        """
        Returns the URL providing the data source for this FeatureSet.
        """
        return self._dbFilePath

    def getExpressionLevels(
            self, threshold=0.0, names=[], startIndex=0, maxResults=0):
        """
        Returns the list of ExpressionLevels in this RNA Quantification.
        """
        rnaQuantificationId = self.getLocalId()
        with self._db as dataSource:
            expressionsReturned = dataSource.searchExpressionLevelsInDb(
                rnaQuantificationId,
                names=names,
                threshold=threshold,
                startIndex=startIndex,
                maxResults=maxResults)
            expressionLevels = [
                SqliteExpressionLevel(self, expressionEntry) for
                expressionEntry in expressionsReturned]
            return expressionLevels

    def getExpressionLevel(self, compoundId):
        expressionId = compoundId.expression_level_id
        with self._db as dataSource:
            expressionReturned = dataSource.getExpressionLevelById(
                expressionId)
        return SqliteExpressionLevel(self, expressionReturned)


class SqliteRnaBackend(sqlite_backend.SqliteBackedDataSource):
    """
    Defines an interface to a sqlite DB which stores all RNA quantifications
    in the dataset.
    """
    def __init__(self, rnaQuantSqlFile="ga4gh-rnaQuant.db"):
        super(SqliteRnaBackend, self).__init__(rnaQuantSqlFile)

    def searchRnaQuantificationsInDb(
            self, rnaQuantificationId=""):
        """
        :param rnaQuantificationId: string restrict search by id
        :return an array of dictionaries, representing the returned data.
        """
        sql = ("SELECT * FROM RnaQuantification")
        sql_args = ()
        if len(rnaQuantificationId) > 0:
            sql += " WHERE id = ? "
            sql_args += (rnaQuantificationId,)
        query = self._dbconn.execute(sql, sql_args)
        try:
            return sqlite_backend.iterativeFetch(query)
        except AttributeError:
            raise exceptions.RnaQuantificationNotFoundException(
                rnaQuantificationId)

    def getRnaQuantificationById(self, rnaQuantificationId):
        """
        :param rnaQuantificationId: the RNA Quantification ID
        :return: dictionary representing an RnaQuantification object,
            or None if no match is found.
        """
        sql = ("SELECT * FROM RnaQuantification WHERE id = ?")
        query = self._dbconn.execute(sql, (rnaQuantificationId,))
        try:
            return sqlite_backend.fetchOne(query)
        except AttributeError:
            raise exceptions.RnaQuantificationNotFoundException(
                rnaQuantificationId)

    def searchExpressionLevelsInDb(
            self, rnaQuantId, names=[], threshold=0.0, startIndex=0,
            maxResults=0):
        """
        :param rnaQuantId: string restrict search by quantification id
        :param threshold: float minimum expression values to return
        :return an array of dictionaries, representing the returned data.
        """
        sql = ("SELECT * FROM Expression WHERE "
               "rna_quantification_id = ? "
               "AND expression > ? ")
        sql_args = (rnaQuantId, threshold)
        if len(names) > 0:
            sql += "AND name in ("
            sql += ",".join(['?' for name in names])
            sql += ") "
            for name in names:
                sql_args += (name,)
        sql += sqlite_backend.limitsSql(
            startIndex=startIndex, maxResults=maxResults)
        query = self._dbconn.execute(sql, sql_args)
        return sqlite_backend.iterativeFetch(query)

    def getExpressionLevelById(self, expressionId):
        """
        :param expressionId: the ExpressionLevel ID
        :return: dictionary representing an ExpressionLevel object,
            or None if no match is found.
        """
        sql = ("SELECT * FROM Expression WHERE id = ?")
        query = self._dbconn.execute(sql, (expressionId,))
        try:
            return sqlite_backend.fetchOne(query)
        except AttributeError:
            raise exceptions.ExpressionLevelNotFoundException(
                expressionId)


class SimulatedRnaQuantificationSet(AbstractRnaQuantificationSet):
    """
    An RNA Quantification set that doesn't derive from a data store.
    Used mostly for testing.
    """
    def __init__(
            self, parentContainer, localId, numRnaQuantifications=2,
            numExpressionLevels=2):
        super(SimulatedRnaQuantificationSet, self).__init__(
            parentContainer, localId)
        for i in range(numRnaQuantifications):
            localId = "simRnaQ{}".format(i)
            rnaQuantification = SimulatedRnaQuantification(
                self, localId, numExpressionLevels)
            self.addRnaQuantification(rnaQuantification)


class SimulatedRnaQuantification(AbstractRnaQuantification):
    """
    A simulated RNA Quantification
    """
    def __init__(
            self, parentContainer, localId, numExpressionLevels=2):
        super(SimulatedRnaQuantification, self).__init__(
            parentContainer, localId)
        self._expressionLevelIds = []
        self._expressionLevelIdMap = {}
        for i in range(numExpressionLevels):
            localId = "simExpLvl{}".format(i)
            expressionLevel = SimulatedExpressionLevel(self, localId)
            self.addExpressionLevel(expressionLevel)

    def addExpressionLevel(self, expressionLevel):
        id_ = expressionLevel.getId()
        self._expressionLevelIds.append(id_)
        self._expressionLevelIdMap[id_] = expressionLevel

    # TODO this makes very little sense
    def getExpressionLevels(
            self, threshold=0.0, names=[],
            startIndex=0, maxResults=0):  # NOQA
        return [self._expressionLevelIdMap[id_] for
                id_ in self._expressionLevelIds]

    def getExpressionLevel(self, compoundId):
        expressionId = str(compoundId)
        return self._expressionLevelIdMap[expressionId]


class SimulatedExpressionLevel(AbstractExpressionLevel):
    """
    A simulated expression level
    """
    def __init__(self, parentContainer, localId):
        super(SimulatedExpressionLevel, self).__init__(
            parentContainer, localId)
        self._isNormalized = False
