"""
Tests the datarepo module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import tempfile
import unittest

import ga4gh.server.datarepo as datarepo
import ga4gh.server.exceptions as exceptions


prefix = "ga4gh_datarepo_test"


def makeTempFile():
    return tempfile.mkstemp(prefix=prefix)


def makeTempDir():
    return tempfile.mkdtemp(prefix=prefix)


class AbstractDataRepoTest(unittest.TestCase):
    """
    Parent class for data repo tests
    """
    def setUp(self):
        _, self._repoPath = makeTempFile()

    def tearDown(self):
        os.unlink(self._repoPath)


class TestDataRepoVersion(AbstractDataRepoTest):
    """
    Tests the repo schema version is written and read correctly
    and throws an error when there is a version mismatch.
    """
    def testRightVersion(self):
        repo = datarepo.SqlDataRepository(self._repoPath)
        repo.open(datarepo.MODE_WRITE)
        repo.initialise()
        anotherRepo = datarepo.SqlDataRepository(self._repoPath)
        anotherRepo.open(datarepo.MODE_READ)
        self.assertEquals(anotherRepo._schemaVersion, str(repo.version))

    def testWrongVersion(self):
        repo = datarepo.SqlDataRepository(self._repoPath)
        repo.version = datarepo.SqlDataRepository.SchemaVersion(
            "wrong.version")
        repo.open(datarepo.MODE_WRITE)
        repo.initialise()
        anotherRepo = datarepo.SqlDataRepository(self._repoPath)
        with self.assertRaises(
                exceptions.RepoSchemaVersionMismatchException):
            anotherRepo.open(datarepo.MODE_READ)


class TestBadDatabase(AbstractDataRepoTest):
    """
    Tests that errors are thrown when an invalid database is used
    """
    def testDbFileWithoutTables(self):
        repo = datarepo.SqlDataRepository(self._repoPath)
        with self.assertRaises(exceptions.RepoInvalidDatabaseException):
            repo.open(datarepo.MODE_READ)

    def testTextFile(self):
        with open(self._repoPath, 'w') as textFile:
            textFile.write('This is now a text file')
        repo = datarepo.SqlDataRepository(self._repoPath)
        with self.assertRaises(exceptions.RepoInvalidDatabaseException):
            repo.open(datarepo.MODE_READ)


class TestBadDatabaseNoSetup(unittest.TestCase):
    """
    Tests that errors are thrown when an invalid database is used
    (does not use setup/teardown functions)
    """
    def testDirectory(self):
        repoPath = makeTempDir()
        repo = datarepo.SqlDataRepository(repoPath)
        with self.assertRaises(exceptions.RepoInvalidDatabaseException):
            repo.open(datarepo.MODE_READ)

    def testNonexistentFile(self):
        repo = datarepo.SqlDataRepository("aFilePathThatDoesNotExist")
        with self.assertRaises(exceptions.RepoNotFoundException):
            repo.open(datarepo.MODE_READ)
