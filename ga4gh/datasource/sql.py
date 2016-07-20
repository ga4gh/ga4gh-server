"""
Module defining concrete instances of top level containers backed
by SQL databases.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sqlite3
import json

import sqlalchemy
import sqlalchemy.orm as orm

import ga4gh.protocol as protocol
import ga4gh.registry as registry
import ga4gh.exceptions as exceptions


def sqliteRows2dicts(sqliteRows):
    """
    Unpacks sqlite rows as returned by fetchall
    into an array of simple dicts.

    :param sqliteRows: array of rows returned from fetchall DB call
    :return:  array of dicts, keyed by the column names.
    """
    return map(lambda r: dict(zip(r.keys(), r)), sqliteRows)


def sqliteRow2Dict(sqliteRow):
    """
    Unpacks a single sqlite row as returned by fetchone
    into a simple dict.

    :param sqliteRow: single row returned from fetchone DB call
    :return: dictionary corresponding to this row
    """
    return dict(zip(sqliteRow.keys(), sqliteRow))


def limitsSql(pageToken=0, pageSize=None):
    """
    Takes parsed pagination data, spits out equivalent SQL 'limit' statement.

    :param pageToken: starting row position,
        can be an int or string containing an int
    :param pageSize: number of records requested for this transaciton,
        can be an int or string containing an int
    :return: SQL 'limit' statement string to append to query.
    """
    if not pageToken:
        pageToken = 0
    if pageSize is not None or pageToken > 0:
        start = int(pageToken)
        end = start + int(pageSize)
        return " LIMIT {}, {}".format(start, end)
    else:
        return ""


def _whereClauseSql(**whereClauses):
    """
    Takes parsed search query parameters,
    produces equivalent SQL 'where' clause.

    :param whereClauses: key-value pairs of column names to limit by,
        and values to limit them to.
    :return: corresponding SQLite 'where' clause string ready to paste
        into query.
    """
    if whereClauses is not None:
        # wc is an array of "key ='value'" strings from whereClauses,
        # with all entries where value = None removed.
        wc = ["{} = '{}'".format(k, whereClauses[k])
              for k in whereClauses.keys()
              if whereClauses[k] is not None]
        if len(wc) > 0:
            return " WHERE " + " AND ".join(wc)
    else:
        return ""


class SqliteBackedDataSource(object):
    """
    Abstract class that sets up a SQLite database source
    as a context-managed data source.
    Client code of a subclass can then look something as follows:

    def search<DataModel>(self, queryParam1=val1, queryParam2=val2,
                        assemblyId=None, pageToken=tok, pageSize=size):

        limits = _makeLimits(tok, size)
        with <dataSourceModule>.<DataSource>(self._dbFile) as dataSource:
            count = dataSource.<customSQLCountMethod>()
            datasetDicts = dataSource.<customSQLSearchMethod>(limits)
        outData = []
        for dict in datasetDicts:
            outDataItem = protocol.<DataRecord>()
            outDataItem.id = dict['ID']
            outDataItem.<somethingElse> = dict['<somethingElse_columnName>']
            outData.append(outDataItem)
        return count, outData
    """
    def __init__(self, dbFile):
        """
        :param dbFile: string holding the full path to the database file.
        """
        self._dbFile = dbFile

    def __enter__(self):
        self._dbconn = sqlite3.connect(self._dbFile)
        # row_factory setting is magic pixie dust to retrieve rows
        # as dictionaries. sqliteRows2dict relies on this.
        self._dbconn.row_factory = sqlite3.Row
        return self

    def __exit__(self, type, value, traceback):
        self._dbconn.close()


"""
For this implementation, `featureSetId` is required, while `parentId`
is optional, and filters the features within the requested `featureSetId`
by their parent.

Only return features on the reference with this name. Genomic positions
are non-negative integers less than reference length.
Requests spanning the join of circular genomes are represented as two
requests one on each side of the join (position 0) end is also required
If specified, this query matches only annotations which match one of the
provided feature types.
For now do not use the features array in search

GFF3 data is represented by rows in a single table, named FEATURE.
The columns of the FEATURE table correspond to the columns of a GFF3,
with three additional columns prepended representing the ID
of this feature, the ID of its parent (if any), and a whitespace
separated array of its child IDs.

_featureColumns pairs represent the ordered (column_name, column_type).
"""
_featureColumns = [
    ('id', 'TEXT'),  # a synthetic principal key generated on ETL
    ('parent_id', 'TEXT'),
    ('child_ids', 'TEXT'),
    ('reference_name', 'TEXT'),
    ('source', 'TEXT'),
    ('type', 'TEXT'),  # corresponds to featureType, an ontology term
    ('start', 'INT'),
    ('end', 'INT'),
    ('score', 'REAL'),
    ('strand', 'TEXT'),  # limited to one of '+'/'-' or none
    ('name', 'TEXT'),  # the "ID" as found in GFF3, or '' if none
    ('gene_name', 'TEXT'),  # as found in GFF3 attributes
    ('transcript_name', 'TEXT'),  # as found in GFF3 attributes
    ('attributes', 'TEXT')]  # JSON encoding of attributes dict


class Gff3DbBackend(SqliteBackedDataSource):
    """
    Notes about the current implementation:
    For this implementation, `featureSetId` is required, while `parentId`
    is optional, and filters the features within the requested `featureSetId`
    by their parent.

    Genomic positions are non-negative integers less than reference length.
    Requests spanning the join of circular genomes are represented as two
    requests one on each side of the join (position 0)
    """

    def __init__(self, dbFile):
        super(Gff3DbBackend, self).__init__(dbFile)
        self.featureColumnNames = [f[0] for f in _featureColumns]
        self.featureColumnTypes = [f[1] for f in _featureColumns]

    def countFeaturesSearchInDb(self, request):
        _, sql, sql_args = self.featuresQuery(request)
        query = self._dbconn.execute(sql, sql_args)
        return (query.fetchone())[0]

    def searchFeaturesInDb(self, request):
        # TODO: Refactor out common bits of this and the above count query.
        sql, _, sql_args = self.featuresQuery(request)
        sql += limitsSql(request.page_token, request.page_size)
        query = self._dbconn.execute(sql, sql_args)
        return sqliteRows2dicts(query.fetchall())

    def featuresQuery(self, request):
        # TODO: Optimize by refactoring out string concatenation
        sql = ""
        sql_rows = "SELECT * FROM FEATURE WHERE id > 1 "
        sql_count = "SELECT COUNT(*) FROM FEATURE WHERE id > 1 "
        sql_args = ()
        if request.name:
            sql += "AND name = ? "  # compare this to query start
            sql_args += (request.name,)
        if request.gene_symbol:
            sql += "AND gene_name = ? "  # compare this to query start
            sql_args += (request.gene_symbol,)
        if request.start:
            sql += "AND end > ? "  # compare this to query start
            sql_args += (request.start,)
        if request.end:
            sql += "AND start < ? "  # and this to query end
            sql_args += (request.end,)
        if request.reference_name:
            sql += "AND reference_name = ?"
            sql_args += (request.reference_name,)
        if request.parent_id:
            sql += "AND parent_id = ? "
            sql_args += (request.parent_id,)
        if len(request.feature_types) > 0:
            sql += "AND type IN ("
            sql += ", ".join(["?", ] * len(request.feature_types))
            sql += ") "
            sql_args += tuple(request.feature_types)
        sql_rows += sql
        sql_rows += "ORDER BY reference_name, start, end ASC "
        sql_count += sql
        return sql_rows, sql_count, sql_args

    def getFeatureById(self, featureId):
        """
        Fetch feature by featureID.

        :param featureId: the FeatureID as found in GFF3 records
        :return: dictionary representing a feature object,
            or None if no match is found.
        """
        sql = "SELECT * FROM FEATURE WHERE id = ?"
        query = self._dbconn.execute(sql, (featureId,))
        ret = query.fetchone()
        if ret is None:
            return None
        return sqliteRow2Dict(ret)


#####################################################################
#
# Features
#
#####################################################################


class Gff3DbFeatureSet(registry.FeatureSet):
    """
    Class to directly read sequence annotation features from GFF3 files.
    Tests basic access, not to be used in production.
    """
    __tablename__ = 'Gff3DbFeatureSet'

    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('FeatureSet.id'),
        primary_key=True)

    ontology_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Ontology.id"),
        nullable=False)
    ontology = orm.relationship("Ontology", single_parent=True)

    data_url = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'Gff3DbFeatureSet',
    }

    def __init__(self, name, data_url):
        super(Gff3DbFeatureSet, self).__init__(name)
        self.data_url = data_url
        # TODO verify something about this DB???

    def getFeature(self, compoundId):
        """
        Returns a protocol.Feature object corresponding to a compoundId
        :param compoundId: a datamodel.FeatureCompoundId object
        :return: a Feature object.
        :raises: exceptions.ObjectWithIdNotFoundException if invalid
            compoundId is provided.
        """
        featureId = long(compoundId.featureId)
        with self._db as dataSource:
            featureReturned = dataSource.getFeatureById(featureId)

        if featureReturned is None:
            raise exceptions.ObjectWithIdNotFoundException(compoundId)
        else:
            gaFeature = self._gaFeatureForFeatureDbRecord(featureReturned)
            return gaFeature

    def getCompoundIdForFeatureId(self, featureId):
        """
        Returns server-style compound ID for an internal featureId.

        :param long featureId: id of feature in database
        :return: string representing ID for the specified GA4GH protocol
            Feature object in this FeatureSet.
        """
        return str(featureId)
        # if featureId is not None and featureId != "":
        #     compoundId = datamodel.FeatureCompoundId(
        #         self.getCompoundId(), str(featureId))
        # else:
        #     compoundId = ""
        # return str(compoundId)

    def _gaFeatureForFeatureDbRecord(self, feature):
        """
        :param feature: The DB Row representing a feature
        :return: the corresponding GA4GH protocol.Feature object
        """
        gaFeature = protocol.Feature()
        gaFeature.id = self.getCompoundIdForFeatureId(feature['id'])
        if feature.get('parent_id'):
            gaFeature.parent_id = str(feature['parent_id'])
        else:
            gaFeature.parent_id = ""
        gaFeature.feature_set_id = str(self.id)
        gaFeature.reference_name = feature.get('reference_name')
        gaFeature.start = feature.get('start')
        gaFeature.end = feature.get('end')
        gaFeature.name = feature.get('name')
        if feature.get('strand', '') == '-':
            gaFeature.strand = protocol.NEG_STRAND
        else:
            # default to positive strand
            gaFeature.strand = protocol.POS_STRAND
        gaFeature.child_ids.extend(
            map(str, json.loads(feature['child_ids'])))
        gaFeature.feature_type.CopyFrom(
            self.ontology.get_ga_term_by_name(feature['type']))
        attributes = json.loads(feature['attributes'])
        # TODO: Identify which values are ExternalIdentifiers and OntologyTerms
        for key in attributes:
            for v in attributes[key]:
                gaFeature.attributes.vals[key].values.add().string_value = v
        if 'gene_name' in attributes and len(attributes['gene_name']) > 0:
            gaFeature.gene_symbol = attributes['gene_name'][0]
        return gaFeature

    def run_search(self, request, response_builder):
        db = Gff3DbBackend(self.data_url)
        with db as dataSource:
            featuresCount = dataSource.countFeaturesSearchInDb(request)
            featuresReturned = dataSource.searchFeaturesInDb(request)
        offset = 0
        # TODO the pagetoken should only be parsed once. We should follow the
        # logic used in the backend for SQL queries.
        if request.page_token:
            offset = int(request.page_token)
        for featureRecord in featuresReturned:
            offset += 1
            gaFeature = self._gaFeatureForFeatureDbRecord(featureRecord)
            response_builder.addValue(gaFeature)
            if response_builder.isFull():
                break
        if offset < featuresCount:
            response_builder.setNextPageToken(str(offset))
