"""
The GA4GH data model. Defines all the methods required to translate
data in existing formats into GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import base64
import collections
import glob
import os

import ga4gh.exceptions as exceptions


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


class CompoundId(object):
    """
    Base class for an id composed of several different parts, separated
    by a separator. Each compound ID consists of a set of fields, each
    of which corresponds to a local ID in the data hierarchy. For example,
    we might have fields like ["dataset", "variantSet"] for a variantSet.
    These are available as cid.dataset, and cid.variantSet.  The actual IDs
    of the containing objects can be obtained using the corresponding
    like cid.datasetId and cid.variantSetId.
    """
    separator = ':'
    fields = []
    """
    The fields that the compound ID is composed of. These are parsed and
    made available as attributes on the object.
    """
    containerIds = []
    """
    The fields of the ID form a breadcrumb trail through the data
    hierarchy, and successive prefixes provide the IDs for objects
    further up the tree. This list is a set of tuples giving the
    name and length of a given prefix forming an identifier.
    """

    def __init__(self, parentCompoundId, *localIds):
        """
        Allocates a new CompoundId for the specified parentCompoundId and
        local identifiers. This compoundId inherits all of the fields and
        values from the parent compound ID, and must have localIds
        corresponding to its fields. If no parent id is present,
        parentCompoundId should be set to None.
        """
        index = 0
        if parentCompoundId is not None:
            for field in parentCompoundId.fields:
                setattr(self, field, getattr(parentCompoundId, field))
                index += 1
        for field, localId in zip(self.fields[index:], localIds):
            setattr(self, field, str(localId))
        if len(localIds) != len(self.fields) - index:
            raise ValueError(
                "Incorrect number of fields provided to instantiate ID")
        for idFieldName, prefix in self.containerIds:
            values = [getattr(self, f) for f in self.fields[:prefix + 1]]
            containerId = self.separator.join(values)
            obfuscated = self.obfuscate(containerId)
            setattr(self, idFieldName, obfuscated)

    def __str__(self):
        values = [getattr(self, f) for f in self.fields]
        compoundIdStr = self.separator.join(values)
        return self.obfuscate(compoundIdStr)

    @classmethod
    def parse(cls, compoundIdStr):
        """
        Parses the specified compoundId string and returns an instance
        of this CompoundId class.

        :raises: An ObjectWithIdNotFoundException if parsing fails. This is
        because this method is a client-facing method, and if a malformed
        identifier (under our internal rules) is provided, the response should
        be that the identifier does not exist.
        """
        if not isinstance(compoundIdStr, basestring):
            raise exceptions.BadIdentifierException(compoundIdStr)
        try:
            deobfuscated = cls.deobfuscate(compoundIdStr)
        except TypeError:
            # When a string that cannot be converted to base64 is passed
            # as an argument, b64decode raises a TypeError. We must treat
            # this as an ID not found error.
            raise exceptions.ObjectWithIdNotFoundException(compoundIdStr)
        try:
            splits = deobfuscated.split(cls.separator)
        except UnicodeDecodeError:
            # Sometimes base64 decoding succeeds but we're left with
            # unicode gibberish. This is also and IdNotFound.
            raise exceptions.ObjectWithIdNotFoundException(compoundIdStr)
        if len(splits) != len(cls.fields):
            raise exceptions.ObjectWithIdNotFoundException(compoundIdStr)
        return cls(None, *splits)

    @classmethod
    def obfuscate(cls, idStr):
        """
        Mildly obfuscates the specified ID string in an easily reversible
        fashion. This is not intended for security purposes, but rather to
        dissuade users from depending on our internal ID structures.
        """
        return base64.b64encode(idStr)

    @classmethod
    def deobfuscate(cls, idStr):
        """
        Reverses the obfuscation done by the :meth:`obfuscate` method.
        """
        return base64.b64decode(idStr)


class ReferenceSetCompoundId(CompoundId):
    """
    The compound ID for reference sets.
    """
    fields = ['referenceSet']
    containerIds = [('referenceSetId', 0)]


class ReferenceCompoundId(ReferenceSetCompoundId):
    """
    The compound id for a reference
    """
    fields = ReferenceSetCompoundId.fields + ['reference']


class DatasetCompoundId(CompoundId):
    """
    The compound id for a data set
    """
    fields = ['dataset']
    containerIds = [('datasetId', 0)]


class VariantSetCompoundId(DatasetCompoundId):
    """
    The compound id for a variant set
    """
    fields = DatasetCompoundId.fields + ['variantSet']
    containerIds = DatasetCompoundId.containerIds + [('variantSetId', 1)]


class VariantSetMetadataCompoundId(VariantSetCompoundId):
    """
    The compound id for a variant set
    """
    fields = VariantSetCompoundId.fields + ['key']
    containerIds = VariantSetCompoundId.containerIds + [
        ('variantSetMetadataId', 1)]


class VariantCompoundId(VariantSetCompoundId):
    """
    The compound id for a variant
    """
    fields = VariantSetCompoundId.fields + ['referenceName', 'start', 'md5']


class CallSetCompoundId(VariantSetCompoundId):
    """
    The compound id for a callset
    """
    fields = VariantSetCompoundId.fields + ['name']


class ReadGroupSetCompoundId(DatasetCompoundId):
    """
    The compound id for a read group set
    """
    fields = DatasetCompoundId.fields + ['readGroupSet']
    containerIds = DatasetCompoundId.containerIds + [('readGroupSetId', 1)]


class ReadGroupCompoundId(ReadGroupSetCompoundId):
    """
    The compound id for a read group
    """
    fields = ReadGroupSetCompoundId.fields + ['readGroup']
    containerIds = ReadGroupSetCompoundId.containerIds + [('readGroupId', 2)]


class ExperimentCompoundId(ReadGroupCompoundId):
    """
    The compound id for an experiment
    """
    fields = ReadGroupCompoundId.fields + ['experiment']
    containerIds = ReadGroupCompoundId.containerIds + [('experimentId', 3)]


class ReadAlignmentCompoundId(ReadGroupCompoundId):
    """
    The compound id for a read alignment
    """
    fields = ReadGroupCompoundId.fields + ['readAlignment']


class DatamodelObject(object):
    """
    Superclass of all datamodel types. A datamodel object is a concrete
    representation of some data, either a single observation (such as a
    read) or an aggregated set of related observations (such as a dataset).
    Every datamodel object has an ID and a localId. The ID is an identifier
    which uniquely idenfifies the object within a server instance. The
    localId is a name that identifies the object with a given its
    parent container.
    """

    compoundIdClass = None
    """ The class for compoundIds. Must be set in concrete subclasses.  """

    def __init__(self, parentContainer, localId):
        self._parentContainer = parentContainer
        self._localId = localId
        parentId = None
        if parentContainer is not None:
            parentId = parentContainer.getCompoundId()
        self._compoundId = self.compoundIdClass(parentId, localId)

    def getId(self):
        """
        Returns the string identifying this DatamodelObject within the
        server.
        """
        return str(self._compoundId)

    def getCompoundId(self):
        """
        Returns the CompoundId instance that identifies this object
        within the server.
        """
        return self._compoundId

    def getLocalId(self):
        """
        Returns the localId of this DatamodelObject. The localId of a
        DatamodelObject is a name that identifies it within its parent
        container.
        """
        return self._localId

    def getParentContainer(self):
        """
        Returns the parent container for this DatamodelObject. This the
        object that is one-level above this object in the data hierarchy.
        For example, for a Variant this is the VariantSet that it belongs
        to.
        """
        return self._parentContainer


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
    def sanitizeAlignmentFileFetch(cls, start=None, end=None):
        if start is not None:
            start = cls.sanitizeInt(
                start, cls.samMin, cls.samMaxStart, 'start')
        if end is not None:
            end = cls.sanitizeInt(end, cls.samMin, cls.samMaxEnd, 'end')
        if start is not None and end is not None:
            cls.assertValidRange(start, end, 'start', 'end')
        return start, end

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
        if not isinstance(attr, (int, long)):
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
        if numDataFiles == 0:
            raise exceptions.EmptyDirException(dataDir, patterns)

    def getFileHandle(self, dataFile):
        return fileHandleCache.getFileHandle(dataFile, self.openFile)
