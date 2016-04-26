"""
Tests the datarepo module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import tempfile
import unittest

import ga4gh.datarepo as datarepo
import ga4gh.exceptions as exceptions


class TestDataRepoVersion(unittest.TestCase):
    """
    Tests the repo schema version is written and read correctly
    and throws an error when there is a version mismatch.
    """
    def setUp(self):
        fd, self._repoPath = tempfile.mkstemp(prefix="ga4gh_datarepo_test")

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

    def tearDown(self):
        os.unlink(self._repoPath)
