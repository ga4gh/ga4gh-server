"""
Module responsible for translating feature expression data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import string

import ga4gh.datamodel as datamodel
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
import ga4gh.sqliteBackend as sqliteBackend


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
        super(ExpressionLevel, self).__init__(parentContainer, record[0])
        self._id = record[0]
        self._annotationId = record[1]
        self._expression = record[2]
        self._featureGroupId = record[3]
        self._isNormalized = record[4]
        self._rawReadCount = record[5]
        self._score = record[6]
        self._units = record[7]

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


class Characterization(datamodel.DatamodelObject):
    """
    Class representing a single Characterization in the GA4GH data model.
    """
    compoundIdClass = datamodel.CharacterizationCompoundId

    def __init__(self, parentContainer, record):
        # TODO this is still using the human readableId as the _id value
        # (localID)
        super(Characterization, self).__init__(parentContainer, record[0])
        self._analysisId = record[0]
        self._complexity = float(record[1])
        self._exonicFraction = float(record[2])
        self._fractionMapped = float(record[3])
        self._intergenicFraction = float(record[4])
        self._intronicFraction = float(record[5])

    def toProtocolElement(self):
        protocolElement = protocol.Characterization
        protocolElement.analysisId = self._analysisId
        protocolElement.complexity = self._complexity
        protocolElement.exonicFraction = self._exonicFraction
        protocolElement.fractionMapped = self._fractionMapped
        protocolElement.intergenicFraction = self._intergenicFraction
        protocolElement.intronicFraction = self._intronicFraction
        return protocolElement


class ReadCounts(datamodel.DatamodelObject):
    """
    Class representing a single ReadCounts in the GA4GH data model.
    """
    compoundIdClass = datamodel.ReadCountsCompoundId

    def __init__(self, parentContainer, record):
        # TODO this is still using the human readableId as the _id value
        # (localID)
        super(ReadCounts, self).__init__(parentContainer, record[0])
        self._analysisId = record[0]
        self._multiCount = int(record[1])
        self._multiSpliceCount = int(record[2])
        self._totalReadCount = int(record[3])
        self._uniqueCount = int(record[4])
        self._uniqueSpliceCount = int(record[5])

    def toProtocolElement(self):
        protocolElement = protocol.ReadCounts
        protocolElement.analysisId = self._analysisId
        protocolElement.multiCount = self._multiCount
        protocolElement.multiSpliceCount = self._multiSpliceCount
        protocolElement.totalReadCount = self._totalReadCount
        protocolElement.uniqueCount = self._uniqueCount
        protocolElement.uniqueSpliceCount = self._uniqueSpliceCount
        return protocolElement


class FeatureGroup(datamodel.DatamodelObject):
    """
    Class representing a single FeatureGroup in the GA4GH model.
    """
    compoundIdClass = datamodel.FeatureGroupCompoundId

    # TODO: this is just a first pass stub to get working
    # - need to formalize input data
    def __init__(self, parentContainer, localId, rnaQuantId):
        super(FeatureGroup, self).__init__(parentContainer, localId)
        self._id = localId
        self._analysisId = rnaQuantId
        self.name = localId

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
        # self._rnaQuantificationId = localId
        self._name = localId
        self._expressionLevelIdMap = {}
        self._expressionLevelIds = []
        self._featureGroupNames = []
        self._featureGroupIdMap = {}
        self._featureGroupIds = []
        self._characterization = None
        self._readCounts = None

    def toProtocolElement(self):
        """
        Converts this rnaQuant into its GA4GH protocol equivalent.
        """
        protocolElement = protocol.RnaQuantification()
        protocolElement.annotationIds = self._annotationIds
        protocolElement.description = self._description
        protocolElement.id = self._id
        protocolElement.name = self._name
        protocolElement.readGroupId = self._readGroupId
        return protocolElement

    def addRnaQuantMetadata(self, fields):
        """
        input list elements are:
        Id, annotations, description, name, readGroupId
        where annotation is a comma separated list
        """
        self._id = fields[0]
        self._annotationIds = fields[1].split(',')
        self._description = fields[2]
        self._name = fields[3]
        self._readGroupId = fields[4]

    def addExpressionLevel(self, record):
        """
        Adds a ExpressionLevel for the specified sample name.
        """
        expressionLevel = ExpressionLevel(self, record)
        expressionLevelId = expressionLevel.getId()
        self._expressionLevelIdMap[expressionLevelId] = expressionLevel
        self._expressionLevelIds.append(expressionLevelId)

    def addFeatureGroup(self, localId):
        """
        Adds a FeatureGroup for the specified name.
        """
        if localId not in self._featureGroupNames:
            featureGroup = FeatureGroup(self, localId, self.getId())
            featureGroupId = featureGroup.getId()
            self._featureGroupIdMap[featureGroupId] = featureGroup
            self._featureGroupIds.append(featureGroupId)
            self._featureGroupNames.append(localId)

    def getCharacterization(self):
        return self._characterization

    def getReadCounts(self):
        return self._readCounts


class RNASeqResult(AbstractRNAQuantification):
    """
    Class representing a single RnaQuantification in the GA4GH data model.
    """

    def __init__(self, parentContainer, localId, rnaQuantDataPath):
        super(RNASeqResult, self).__init__(parentContainer, localId)
        self._rnaQuantDataPath = rnaQuantDataPath
        self.parseRnaQuantMetadataFile()
        self.addCharacterization()
        self.addReadCounts()
        self._expressionLevelFile = os.path.join(
            rnaQuantDataPath, "expression.table")

    def parseRnaQuantMetadataFile(self):
        """
        input is tab file with no header.  Columns are:
        Id, annotations, description, name, readGroupId
        where annotation is a comma separated list
        """
        rnaQuantificationFile = os.path.join(self._rnaQuantDataPath,
                                             "rnaseq.table")
        with open(rnaQuantificationFile, "r") as rnaQuantificationData:
            quantData = rnaQuantificationData.readline()
            fields = quantData.strip().split("\t")
        self.addRnaQuantMetadata(fields)

    def parseExpressionData(self):
        """
        input is tab file with no header.  Columns are:
        id, annotationId, expression, featureGroupId,
        isNormalized, rawReadCount, score, units
        expressionLevelId is not None: return only the specific expressionLevel
        object
        featureGroupId is not None: return all in that group

        Parsing the expression data file at server startup took too long so
        the procedure was moves here.  Will have overhead for these "get them
        requests" but that is on user side and not server startup.

        So we check to see if we've done it and if not do it on the first
        request.
        """
        if self._expressionLevelIds == []:
            with open(self._expressionLevelFile, "r") as expressionLevelData:
                for expressionData in expressionLevelData.readlines():
                    fields = expressionData.strip().split("\t")
                    self.addExpressionLevel(fields)
                    self.addFeatureGroup(fields[3])

    def getExpressionLevels(self, featureGroupId=None):
        """
        Returns the list of ExpressionLevels in this RNA Quantification.
        """
        self.parseExpressionData()
        if featureGroupId is None:
            return [self._expressionLevelIdMap[id_] for
                    id_ in self._expressionLevelIds]
        else:
            # TODO: need to clean up the name/ID thing as this is using the id
            # from the csv and not the compoundId of the feature group
            return [self._expressionLevelIdMap[id_] for
                    id_ in self._expressionLevelIds if
                    (self._expressionLevelIdMap[id_].getFeatureGroup() ==
                     featureGroupId)]

    def getExpressionLevel(self, id_, featureGroupId=None):
        """
        Returns a ExpressionLevel with the specified id, or raises a
        ExpressionLevelNotFoundException if it does not exist.
        """
        self.parseExpressionData()
        if id_ not in self._expressionLevelIdMap:
            raise exceptions.ExpressionLevelNotFoundException(id_)
        if featureGroupId is not None:
            if (self._expressionLevelIdMap[id_].getFeatureGroup() !=
                    featureGroupId):
                raise exceptions.ExpressionLevelNotFoundException(id_)
        return self._expressionLevelIdMap[id_]

    def addCharacterization(self):
        """
        input is tab file with no header.  Columns are:
        analysisId, complexity, exonicFraction, fractionMapped,
        intergenicFraction, intronicFraction
        """
        characterizationFile = os.path.join(self._rnaQuantDataPath,
                                            "dist.table")
        try:
            with open(characterizationFile, "r") as characterizationData:
                quantCharacterization = characterizationData.readline()
                fields = quantCharacterization.split("\t")
        except IOError:
            self._characterization = None
            return
        self._characterization = Characterization(self, fields)

    def addReadCounts(self):
        """
        input is tab file with no header.  Columns are:
        analysisId, multiCount, multiSpliceCount, totalReadCount, uniqueCount,
        uniqueSpliceCount
        """
        readCountFile = os.path.join(self._rnaQuantDataPath, "counts.table")
        try:
            with open(readCountFile, "r") as readCountData:
                countData = readCountData.readline()
                fields = countData.split("\t")
        except IOError:
            self._readCounts = None
            return
        self._readCounts = ReadCounts(self, fields)

    def getFeatureGroups(self):
        """
        Returns the list of FeatureGroups in this RNA Quantification.

        for now the feature group data is autogenerated by examining the
        relevant expression data file
        """
        self.parseExpressionData()
        return [self._featureGroupIdMap[id_] for id_ in self._featureGroupIds]

    def getFeatureGroup(self, id_):
        """
        for now the feature group data is autogenerated by examining the
        relevant expression data file
        """
        self.parseExpressionData()
        if id_ not in self._featureGroupIdMap.keys():
            raise exceptions.FeatureGroupNotFoundException(id_)
        return self._featureGroupIdMap[id_]


_rnaQuantColumns = [('id', 'TEXT'),
                    ('annotationIds', 'TEXT'),
                    ('description', 'TEXT'),
                    ('name', 'TEXT'),
                    ('readGroupId', 'TEXT')
                   ]


_expressionColumns = [('id', 'TEXT'),
                      ('name', 'TEXT'),
                      ('annotationId', 'TEXT'),
                      ('expression', 'REAL'),
                      ('featureGroupId', 'TEXT'),
                      ('isNormalized', 'BOOLEAN'),
                      ('rawReadCount', 'REAL'),
                      ('score', 'REAL'),
                      ('units', 'TEXT')
                     ]


class SqliteRNASeqResult(sqliteBackend.SqliteBackedDataSource):
    """
    Defines an interface to a sqlite DB which stores all RNA quantifications
    in the dataset.
    """
    def __init__(self, rnaQuantSqlFile="ga4gh-rnaQuant.db"):
        super(SimulatedRNASeqResult, self).__init__(rnaQuantSqlFile)
        self.rnaQuantColumnNames = [f[0] for f in _rnaQuantColumns]
        self.rnaQuantColumnTypes = [f[1] for f in _rnaQuantColumns]
        self.expressionColumnNames = [f[0] for f in _expressionColumns]
        self.expressionColumnTypes = [f[1] for f in _expressionColumns]

    def getRnaQuantifications(self, pageToken=0, pageSize=None, **query):
        """
        Returns the list of RnaQuantifications in this DB.  Raises a
        RnaQuantificationNotFoundException if id is specified and does not
        exist.
        """
        sql = "SELECT * FROM rnaQuantification "
        whereClauses = []
        for col in query:
            if col in self.featureColumnNames:
                colIdx = self.featureColumnNames.index(col)
                colType = self.featureColumnTypes[colIdx]
                colVal = query[col]
                # simple input sanity check
                if "'" in colVal:
                    throw(ga4gh.exceptions.BadRequestException)
                whereClauses.append("{} = '{}'".format(col, colVal))

        if len(whereClauses) > 0:
            sql += "WHERE {} ".format(" AND ".join(whereClauses))
        sql += "ORDER BY start, id "
        sql += sqliteBackend.limitsSql(pageToken, pageSize)
        query = self._dbconn.execute(sql)
        return sqliteBackend.sqliteRows2dicts(query.fetchall())

    def getExpressionLevels(self, pageToken=0, pageSize=None, **query):
        """
        Returns the list of ExpressionLevels in this RNA Quantification.
        """
        sql = "SELECT * FROM expression "
        whereClauses = []
        for col in query:
            if col in self.featureColumnNames:
                colIdx = self.featureColumnNames.index(col)
                colType = self.featureColumnTypes[colIdx]
                colVal = query[col]
                # simple input sanity check
                if colType is "REAL":
                    colVal = float(colVal)
                     # expression query value is a minimum threshold
                    if col is "expression":
                        whereClauses.append("expression >= {}".format(colVal))
                    else:
                        whereClauses.append("{} = {}".format(col, colVal))
                else:  # TEXT of some sort
                    if "'" in colVal:
                        throw(ga4gh.exceptions.BadRequestException)
                    whereClauses.append("{} = '{}'".format(col, colVal))
        if len(whereClauses) > 0:
            sql += "WHERE {} ".format(" AND ".join(whereClauses))
        sql += "ORDER BY start, id "
        sql += sqliteBackend.limitsSql(pageToken, pageSize)
        query = self._dbconn.execute(sql)
        return sqliteBackend.sqliteRows2dicts(query.fetchall())

    def getExpressionLevel(self, id_, featureGroupId=None):
        """
        Returns a ExpressionLevel with the specified id, or raises a
        ExpressionLevelNotFoundException if it does not exist.
        """
        results = self.getExpressionLevels(featureGroupId=featureGroupId,
                                           expressionId=id_)
        if len(results) == 0:
            raise exceptions.FeatureGroupNotFoundException(id_)
        return results


class SimulatedRNASeqResult(AbstractRNAQuantification):
    """
    An RNA Quantification that doesn't derive from a data store.
    Used mostly for testing.
    """
    def __init__(self, parentContainer, localId, rnaQuantDataPath=""):
        super(SimulatedRNASeqResult, self).__init__(parentContainer, localId)
        self._rnaQuantDataPath = rnaQuantDataPath
        self.generateRnaQuantMetadata()
        self.generateCharacterization()
        self.generateReadCounts()

    def generateCharacterization(self):
        """
            Currently just returns default values.
        """
        fields = ["", 0, 0, 0, 0, 0]
        self._characterization = Characterization(self, fields)

    def generateReadCounts(self):
        """
            Currently just returns default values.
        """
        fields = ["", 0, 0, 0, 0, 0]
        self._readCounts = ReadCounts(self, fields)

    def generateRnaQuantMetadata(self):
        """
            Currently just sets to default values.
        """
        fields = ["", "", "", "", ""]
        self.addRnaQuantMetadata(fields)


class RNASqliteStore(object):
    """
    Defines a sqlite store for RNA data as well as methods for loading and
    modifying the tables.
    """
    def __init__(self, rnaQuantDataPath, sqliteFileName=None):
        if sqliteFileName is not None:
            sqlFilePath = os.path.join(rnaQuantDataPath, sqliteFileName)
            if sqliteFileName in os.listdir(rnaQuantDataPath):
                self._dbConn = sqlite3.connect(sqlFilePath)
                self._cursor = self._dbConn.cursor()
            else:
                self.createNewRepo(sqlFilePath)

    def createNewRepo(self, sqlFilePath):
        self._dbConn = sqlite3.connect(sqlFilePath)
        self._cursor = self._dbConn.cursor()
        self.createTables(self._cursor)
        self._dbConn.commit()

    def createTables(self, cursor):
        # annotationIds is a comma separated list
        cursor.execute('''CREATE TABLE rnaQuantification (
                       id text,
                       annotationIds text,
                       description text,
                       name text,
                       readGroupId text)''')
        cursor.execute('''CREATE TABLE expression (
                       id text,
                       name text,
                       rnaQuantificationId text,
                       annotationId text,
                       expression real,
                       featureGroupId text,
                       isNormalized boolean,
                       rawReadCount real,
                       score real,
                       units text)''')

    def addRNAQuantification(self, datafields):
        """
        Adds an RNAQuantification to the db.  Datafields is a tuple in the order:
        id, annotationIds, description, name, readGroupId
        """
        self._cursor.execute('INSERT INTO rnaQuantification VALUES (?, ?, ?, ?, ?)', datafields)
        self._dbConn.commit()

    def addExpression(self, datafields):
        """
        Adds an Expression to the db.  Datafields is a tuple in the order:
        id, name, rnaQuantificationId, annotationId, expression, featureGroupId, isNormalized, rawReadCount, score, units
        """
        self._cursor.execute('INSERT INTO expression VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', datafields)
        self._dbConn.commit()
