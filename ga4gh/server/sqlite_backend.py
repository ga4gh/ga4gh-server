"""
Library for accessing SQLite-backed data
optimized for typical REST queries as defined by the GA4GH server.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sqlite3


def sqliteRowsToDicts(sqliteRows):
    """
    Unpacks sqlite rows as returned by fetchall
    into an array of simple dicts.

    :param sqliteRows: array of rows returned from fetchall DB call
    :return:  array of dicts, keyed by the column names.
    """
    return map(lambda r: dict(zip(r.keys(), r)), sqliteRows)


def sqliteRowToDict(sqliteRow):
    """
    Unpacks a single sqlite row as returned by fetchone
    into a simple dict.

    :param sqliteRow: single row returned from fetchone DB call
    :return: dictionary corresponding to this row
    """
    return dict(zip(sqliteRow.keys(), sqliteRow))


def limitsSql(startIndex=0, maxResults=0):
    """
    Construct a SQL LIMIT clause
    """
    if startIndex and maxResults:
        return " LIMIT {}, {}".format(startIndex, maxResults)
    elif startIndex:
        raise Exception("startIndex was provided, but maxResults was not")
    elif maxResults:
        return " LIMIT {}".format(maxResults)
    else:
        return ""


default_batch_size = 5  # arbitrary; can tweak later


def iterativeFetch(query, batchSize=default_batch_size):
    """
    Returns rows of a sql fetch query on demand
    """
    while True:
        rows = query.fetchmany(batchSize)
        if not rows:
            break
        rowDicts = sqliteRowsToDicts(rows)
        for rowDict in rowDicts:
            yield rowDict


def fetchOne(query):
    """
    Returns a dict result from the fetch of one query row
    """
    return sqliteRowToDict(query.fetchone())


class SqliteBackedDataSource(object):
    """
    Abstract class that sets up a SQLite database source
    as a context-managed data source.
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
