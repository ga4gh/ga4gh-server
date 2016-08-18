"""
Unit tests for the sql backend
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.sqliteBackend as sqliteBackend
import tests.paths as paths


class SqliteDB(sqliteBackend.SqliteBackedDataSource):

    def __init__(self, dbPath=paths.testDataRepo):
        super(SqliteDB, self).__init__(dbPath)
        self._readGroupSql = "SELECT id, name FROM ReadGroup"

    def ping(self):
        sql = "SELECT 1"
        query = self._dbconn.execute(sql)
        row = query.fetchone()
        result = int(row.keys()[0])
        return result

    def fetchOneMethod(self):
        sql = "SELECT id, name FROM ReadGroup LIMIT 1"
        query = self._dbconn.execute(sql)
        rowDict = sqliteBackend.fetchOne(query)
        return rowDict

    def getReadGroupRows(self):
        sql = self._readGroupSql
        query = self._dbconn.execute(sql)
        rows = query.fetchall()
        return rows

    def iterativeFetchMethod(self):
        sql = self._readGroupSql
        query = self._dbconn.execute(sql)
        iterator = sqliteBackend.iterativeFetch(query, 2)
        return iterator


class TestSqlBackend(unittest.TestCase):

    def setUp(self):
        self._db = SqliteDB()

    def testPing(self):
        result = None
        with self._db as db:
            result = db.ping()
        self.assertEqual(result, 1)

    def testLimitClause(self):
        noArgs = sqliteBackend.limitsSql()
        zeroArgs = sqliteBackend.limitsSql(0, 0)
        self.assertEqual(noArgs, zeroArgs)

        with self.assertRaises(Exception):
            sqliteBackend.limitsSql(startIndex=5)

        limit = sqliteBackend.limitsSql(startIndex=1, maxResults=2)
        self.assertEqual(limit, " LIMIT 1, 2")

        limit = sqliteBackend.limitsSql(maxResults=3)
        self.assertEqual(limit, " LIMIT 3")

    def _testRowDict(self, rowDict):
        self.assertEqual(len(rowDict.keys()), 2)
        self.assertIn("id", rowDict.keys())
        self.assertIn("name", rowDict.keys())
        self.assertIsInstance(rowDict["id"], unicode)
        self.assertIsInstance(rowDict["name"], unicode)

    def testRowToDict(self):
        rows = None
        with self._db as db:
            rows = db.getReadGroupRows()
        row = rows[0]
        rowDict = sqliteBackend.sqliteRowToDict(row)
        self._testRowDict(rowDict)

    def testRowsToDicts(self):
        rows = None
        with self._db as db:
            rows = db.getReadGroupRows()
        rowDicts = sqliteBackend.sqliteRowsToDicts(rows)
        for rowDict in rowDicts:
            self._testRowDict(rowDict)

    def testIterativeFetch(self):
        iterator = None
        with self._db as db:
            iterator = db.iterativeFetchMethod()
            for rowDict in iterator:
                self._testRowDict(rowDict)
            iteratorLen = len(list(db.iterativeFetchMethod()))
            regularLen = len(db.getReadGroupRows())
            self.assertEqual(iteratorLen, regularLen)

    def testFetchOne(self):
        rowDict = None
        with self._db as db:
            rowDict = db.fetchOneMethod()
        self._testRowDict(rowDict)
