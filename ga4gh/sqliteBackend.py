"""
Library for accessing SQLite-backed data
optimized for typical REST queries as defined by the GA4GH server.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sqlite3


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
