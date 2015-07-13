"""
The GA4GH data model. Defines all the methods required to translate
data in existing formats into GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import glob
import os
import json
import tempfile
import shutil
import atexit
import collections

import ga4gh.exceptions as exceptions


def _cleanupHtslibsMess(indexDir):
    """
    Cleanup the mess that htslib has left behind with the index files.
    This is a temporary measure until we get a good interface for
    dealing with indexes for remote files.
    """
    if os.path.exists(indexDir):
        shutil.rmtree(indexDir)


class PysamFileHandleCache(object):
    """
    Cache for opened file handles. We use a deque which has the
    advantage to have push/pop operations in O(1) We always add
    elements on the left of the deque and pop elements from the right.
    When a file is accessed via getFileHandle, its priority gets
    updated, it is put at the "top" of the deque.
    """

    def __init__(self):
        self._cache = collections.deque()
        self._memoTable = dict()
        # Initialize the value even if it will be set up by the config
        self._maxCacheSize = 50

    def setMaxCacheSize(self, size):
        """
        Sets the maximum size of the cache
        """
        if size <= 0:
            raise ValueError(
                "The size of the cache must be a strictly positive value")
        self._maxCacheSize = size

    def _add(self, dataFile, handle):
        """
        Add a file handle to the left of the deque
        """
        self._cache.appendleft((dataFile, handle))

    def _update(self, dataFile, handle):
        """
        Update the priority of the file handle. The element is first
        removed and then added to the left of the deque.
        """
        self._cache.remove((dataFile, handle))
        self._add(dataFile, handle)

    def _removeLru(self):
        """
        Remove the least recently used file handle from the cache.
        The pop method removes an element from the right of the deque.
        Returns the name of the file that has been removed.
        """
        (dataFile, handle) = self._cache.pop()
        handle.close()
        return dataFile

    def getCachedFiles(self):
        """
        Returns all file names stored in the cache.
        """
        return self._memoTable.keys()

    def getFileHandle(self, dataFile, openMethod):
        """
        Returns handle associated to the filename. If the file is
        already opened, update its priority in the cache and return
        its handle. Otherwise, open the file using openMethod, store
        it in the cache and return the corresponding handle.
        """
        if dataFile in self._memoTable:
            handle = self._memoTable[dataFile]
            self._update(dataFile, handle)
            return handle
        else:
            try:
                handle = openMethod(dataFile)
            except ValueError:
                raise exceptions.FileOpenFailedException(dataFile)

            self._memoTable[dataFile] = handle
            self._add(dataFile, handle)
            if len(self._memoTable) > self._maxCacheSize:
                dataFile = self._removeLru()
                del self._memoTable[dataFile]
            return handle


# LRU cache of open file handles
fileHandleCache = PysamFileHandleCache()


class DatamodelObject(object):
    """
    Superclass of all datamodel types
    """
    def __init__(self):
        # TODO move common functionality into this class from subclasses
        pass


class PysamDatamodelMixin(object):
    """
    A mixin class to simplify working with DatamodelObjects based on
    directories of files interpreted using pysam. This mixin is designed
    to work within the DatamodelObject hierarchy.
    """
    samMin = 0
    samMaxStart = 2**30 - 1
    samMaxEnd = 2**30

    vcfMin = -2**31
    vcfMax = 2**31 - 1

    fastaMin = 0
    fastaMax = 2**30 - 1

    rNameMin = 0
    rNameMax = 85

    maxStringLength = 2**10  # arbitrary

    @classmethod
    def sanitizeVariantFileFetch(cls, contig=None, start=None, stop=None):
        if contig is not None:
            contig = cls.sanitizeString(contig, 'contig')
        if start is not None:
            start = cls.sanitizeInt(start, cls.vcfMin, cls.vcfMax, 'start')
        if stop is not None:
            stop = cls.sanitizeInt(stop, cls.vcfMin, cls.vcfMax, 'stop')
        if start is not None and stop is not None:
            cls.assertValidRange(start, stop, 'start', 'stop')
        return contig, start, stop

    @classmethod
    def sanitizeAlignmentFileFetch(
            cls, referenceName=None, start=None, end=None):
        if referenceName is not None:
            referenceName = cls.sanitizeString(referenceName, 'referenceName')
        if start is not None:
            start = cls.sanitizeInt(
                start, cls.samMin, cls.samMaxStart, 'start')
        if end is not None:
            end = cls.sanitizeInt(end, cls.samMin, cls.samMaxEnd, 'end')
        if start is not None and end is not None:
            cls.assertValidRange(start, end, 'start', 'end')
        return referenceName, start, end

    @classmethod
    def sanitizeFastaFileFetch(cls, start=None, end=None):
        if start is not None:
            start = cls.sanitizeInt(
                start, cls.fastaMin, cls.fastaMax, 'start')
        if end is not None:
            end = cls.sanitizeInt(
                end, cls.fastaMin, cls.fastaMax, 'end')
        if start is not None and end is not None:
            cls.assertValidRange(start, end, 'start', 'end')
        return start, end

    @classmethod
    def sanitizeGetRName(cls, referenceId):
        cls.assertInt(referenceId, 'referenceId')
        cls.assertInRange(
            referenceId, cls.rNameMin, cls.rNameMax, 'referenceId')

    @classmethod
    def assertValidRange(cls, start, end, startName, endName):
        if start > end:
            message = "invalid coordinates: {} ({}) " \
                "greater than {} ({})".format(startName, start, endName, end)
            raise exceptions.DatamodelValidationException(message)

    @classmethod
    def assertInRange(cls, attr, minVal, maxVal, attrName):
        message = "invalid {} '{}' outside of range [{}, {}]"
        if attr < minVal:
            raise exceptions.DatamodelValidationException(message.format(
                attrName, attr, minVal, maxVal))
        if attr > maxVal:
            raise exceptions.DatamodelValidationException(message.format(
                attrName, attr, minVal, maxVal))

    @classmethod
    def assertInt(cls, attr, attrName):
        if not isinstance(attr, int):
            message = "invalid {} '{}' not an int".format(attrName, attr)
            raise exceptions.DatamodelValidationException(message)

    @classmethod
    def sanitizeInt(cls, attr, minVal, maxVal, attrName):
        cls.assertInt(attr, attrName)
        if attr < minVal:
            attr = minVal
        if attr > maxVal:
            attr = maxVal
        return attr

    @classmethod
    def sanitizeString(cls, attr, attrName):
        if not isinstance(attr, basestring):
            message = "invalid {} '{}' not a string".format(
                attrName, attr)
            raise exceptions.DatamodelValidationException(message)
        if isinstance(attr, unicode):
            attr = attr.encode('utf8')
        if len(attr) > cls.maxStringLength:
            attr = attr[:cls.maxStringLength]
        return attr

    def _setAccessTimes(self, directoryPath):
        """
        Sets the creationTime and accessTime for this file system based
        DatamodelObject. This is derived from the ctime of the specified
        directoryPath.
        """
        ctimeInMillis = int(os.path.getctime(directoryPath) * 1000)
        self._creationTime = ctimeInMillis
        self._updatedTime = ctimeInMillis

    def _scanDataFiles(self, dataDir, patterns):
        """
        Scans the specified directory for files with the specified globbing
        pattern and calls self._addDataFile for each. Raises an
        EmptyDirException if no data files are found.
        """
        numDataFiles = 0
        for pattern in patterns:
            scanPath = os.path.join(dataDir, pattern)
            for filename in glob.glob(scanPath):
                self._addDataFile(filename)
                numDataFiles += 1
        # This is a temporary workaround to allow us to use htslib's
        # facility for working with remote files. The urls.json is
        # definitely not a good idea and will be replaced later.
        # We make a temporary file for each process so that it
        # downloads its own copy and we are sure it's not overwriting
        # the copy of another process. We then register a cleanup
        # handler to get rid of these files on exit.
        urlSource = os.path.join(dataDir, "urls.json")
        if os.path.exists(urlSource):
            with open(urlSource) as jsonFile:
                urls = json.load(jsonFile)["urls"]
            indexDir = tempfile.mkdtemp(prefix="htslib_mess.")
            cwd = os.getcwd()
            os.chdir(indexDir)
            for url in urls:
                self._addDataFile(url)
                numDataFiles += 1
            os.chdir(cwd)
            atexit.register(_cleanupHtslibsMess, indexDir)
        if numDataFiles == 0:
            raise exceptions.EmptyDirException(dataDir, patterns)

    def getFileHandle(self, dataFile):
        return fileHandleCache.getFileHandle(dataFile, self.openFile)
