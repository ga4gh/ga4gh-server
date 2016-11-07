"""
Tests the file handle cache
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import shutil
import tempfile
import unittest
import uuid

import ga4gh.server.datamodel as datamodel


class TestFileHandleCache(datamodel.PysamFileHandleCache, unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestFileHandleCache, self).__init__()
        unittest.TestCase.__init__(self, *args, **kwargs)

    def setUp(self):
        self._tempdir = tempfile.mkdtemp(prefix="ga4gh_file_cache",
                                         dir=tempfile.gettempdir())

    def _getFileHandle(self, dataFile):
        def openMethod(dataFile):
            return open(dataFile, 'w')
        return self.getFileHandle(dataFile, openMethod)

    def testGetFileHandle(self):
        def genFileName(x):
            return os.path.join(self._tempdir, str(uuid.uuid4()))

        # Set cache size to 9 files max
        self.setMaxCacheSize(9)

        # Build a list of 10 files and add their handles to the cache
        fileList = map(genFileName, range(0, 10))

        for f in fileList:
            handle = self._getFileHandle(f)
            self.assertEquals(self._cache.count((f, handle)), 1)

        self.assertEquals(len(self._memoTable), len(self._cache))

        # Ensure that the first added file has been removed from the cache
        self.assertEquals(filter(lambda x: x[0] == fileList[0], self._cache),
                          [])

        topIndex = len(self._cache) - 1

        # Update priority of this file and ensure it's no longer the
        # least recently used
        self.assertEquals(self._cache[topIndex][0], fileList[1])
        self._getFileHandle(fileList[1])
        self.assertNotEqual(self._cache[topIndex][0], fileList[1])
        self.assertEquals(self._cache[0][0], fileList[1])

    def testSetCacheMaxSize(self):
        self.assertRaises(ValueError, self.setMaxCacheSize, 0)
        self.assertRaises(ValueError, self.setMaxCacheSize, -1)

    def tearDown(self):
        shutil.rmtree(self._tempdir)
