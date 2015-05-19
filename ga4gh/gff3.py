"""
Object representation of GFF3 data and parser for GFF3 files.

See: http://www.sequenceontology.org/gff3.shtml
"""

from __future__ import print_function
from __future__ import unicode_literals

import urllib
import copy
import re
import gzip
import bz2
import collections


class GFF3Exception(Exception):
    "exception associated with GFF3 data"
    def __init__(self, message, fileName=None, lineNumber=None):
        if fileName is not None:
            message = str(fileName) + ":" + str(lineNumber) + ": " + message
        Exception.__init__(self, message)


class Feature(object):
    """One feature notation from a GFF3, missing attribute, as code by `.' in
    the file are stored as None """

    __slots__ = ("seqname", "source", "type", "start", "end", "score",
                 "strand", "frame", "attributes", "gff3Set", "lineNumber",
                 "parents", "children")

    def __init__(self, seqname, source, type, start, end, score, strand,
                 frame, attributes, gff3Set, lineNumber=None):
        """The attributes field is a dict of lists/tupes as attributes maybe
        multi-valued.
        """
        self.seqname = seqname
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = copy.deepcopy(attributes)
        self.gff3Set = gff3Set
        self.lineNumber = lineNumber
        self.parents = set()
        self.children = set()

    @staticmethod
    def __dotIfNone(val):
        return val if val is not None else '.'

    def __attributeStr(self, name):
        "return name=value for a single attributed"
        return name + "=" + ",".join(urllib.quote(v)
                                     for v in self.attributes[name])

    def __attributeStrs(self):
        """return name=value, semi-colon-separated string for attributes, including
        url-style quoting"""
        return ";".join([self.__attributeStr(name)
                         for name in self.attributes.attributes.iterkeys()])

    def __str__(self):
        "return the object as a valid GFF3 record line"
        return "\t".join([self.seqname, self.source, self.type,
                          self.start, self.end, self.__dotIfNone(self.score),
                          self.__dotIfNone(self.strand),
                          self.__dotIfNone(self.frame)] +
                         self.__attributeStrs())
    @property
    def featureId(self):
        return self.getRequiredAttribute("ID")[0]
    
    def getRequiredAttribute(self, name):
        "get a require attributes list of values, error if not found"
        values = self.attributes.get(name)
        if values is None:
            raise GFF3Exception("required attribute not found: " + name,
                                self.gff3Set.fileName, self.lineNumber)
        return values

class Gff3Set(object):
    "A set of GFF3 sequence annotations"
    def __init__(self, fileName=None):
        self.fileName = fileName
        self.roots = set()     # root nodes (those with out parents)
        # index of features by id. GFF3 allows disjoint features with
        # the same id, although this is rarely used.
        self.byFeatureId = collections.defaultdict(list) 

    def add(self, feature):
        "add a feature record"
        self.byFeatureId[feature.featureId].append(feature)

    def __linkFeature(self, feature):
        parentIds = feature.attributes.get("Parent")
        if parentIds is None:
            self.roots.add(feature)
        else:
            for parentId in parentIds:
                self.__linkToParent(feature, parentId)

    def __linkToParent(self, feature, parentId):
        parentParts = self.byFeatureId.get(parentId)
        if parentParts is None:
            raise GFF3Exception("Parent feature does not exist: " + parentId,
                                self.fileName, feature.lineNumber)
        # parent maybe disjoint
        for parentPart in parentParts:
            feature.parents.add(parentPart)
            parentPart.children.add(feature)

    def finish(self):
        "finish loading the set, constructing the tree"
        # features maybe disjoint
        for featureParts in self.byFeatureId.itervalues():
            for feature in featureParts:
                self.__linkFeature(feature)


class Gff3Parser(object):
    """Parser for GFF3 files.  This parse does basic validation, but does not
    fully test for conformance."""

    def __init__(self, fileName):
        """If fileName ends with .gz or .bz2, it will decompressed.  If fh is
        specified, then parse the already opened stream, with file name still
        used for error message, otherwise the file be opened and parsed"""
        self.fileName = fileName
        self.lineNumber = 0

    def __open(self):
        "open input file, optionally with decompression"
        if self.fileName.endswith(".gz"):
            return gzip.open(self.fileName)
        elif self.fileName.endswith(".bz2"):
            return bz2.BZ2File(self.fileName)
        else:
            return open(self.fileName)

    SPLIT_ATTR_RE = re.compile("^([a-zA-Z_]+)=(.+)$")  # parses `attr=val'

    def __parseAttrVal(self, attrStr):
        """returns tuple of tuple of (attr, value), multiple are returned to
        handle multi-value attributes"""
        m = self.SPLIT_ATTR_RE.match(attrStr)
        if m is None:
            raise GFF3Exception("can't parse attribute/value: '" + attrStr +
                                "'", self.fileName, self.lineNumber)
        name = m.group(1)
        val = m.group(2)
        # FIXME: parsing of value is ambiguous. Unquote then comma separate,
        # or coma separate and then unquote? Target attribute space separation
        # is also ambiguous as it has space separated arguments.  It appears
        # like special per attribute name handling is actually required.  Also
        # when to split by comma seems attribute-specific, which is
        # problematic for user-defined. attributes.
        return (name, urllib.unquote(val).split(','))

    SPLIT_ATTR_COL_RE = re.compile("; *")

    def __parseAttrs(self, attrsStr):
        "parse the attributes and values"
        attributes = dict()
        for attrStr in self.SPLIT_ATTR_COL_RE.split(attrsStr):
            name, vals = self.__parseAttrVal(attrStr)
            if name in attributes:
                raise GFF3Exception("duplicated attribute name: " + name,
                                    self.fileName, self.lineNumber)
            attributes[name] = vals
        return attributes

    GFF3_NUM_COLS = 9

    def __parseRecord(self, gff3Set, line):
        row = line.split("\t")
        if len(row) != self.GFF3_NUM_COLS:
            raise GFF3Exception("Wrong number of columns, expected " +
                                str(self.GFF3_NUM_COLS) + ", got " +
                                str(len(row)), self.fileName, self.lineNumber)
        feature = Feature(row[0], row[1], row[2], int(row[3]), int(row[4]),
                          row[5], row[6], row[7], self.__parseAttrs(row[8]),
                          gff3Set, self.lineNumber)
        gff3Set.add(feature)

    # spaces or comment line
    IGNORED_LINE_RE = re.compile("(^[ ]*$)|(^[ ]*#.*$)")

    @staticmethod
    def __ignoredLine(line):
        return Gff3Parser.IGNORED_LINE_RE.search(line) is not None

    def __checkHeader(self, line):
        GFF3_HEADER = "##gff-version 3"
        if line.rstrip() != GFF3_HEADER:
            raise GFF3Exception("First line is not GFF3 header (%s), got: %s"
                                % (GFF3_HEADER, line), self.fileName,
                                self.lineNumber)

    def __parseLine(self, gff3Set, line):
        if self.lineNumber == 1:
            self.__checkHeader(line)
        elif not self.__ignoredLine(line):
            self.__parseRecord(gff3Set, line)

    def parse(self):
        "parse and return a Gff3Set object"
        fh = self.__open()
        try:
            gff3Set = Gff3Set(self.fileName)
            for line in fh:
                self.lineNumber += 1
                self.__parseLine(gff3Set, line[0:-1])
        finally:
            fh.close()
        gff3Set.finish()
        return gff3Set
