"""
Object tree representation of GFF3 data and parser for GFF3 files.

See: http://www.sequenceontology.org/gff3.shtml
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import urllib
import copy
import re
import gzip
import bz2
import collections

GFF3_HEADER = "##gff-version 3"


class GFF3Exception(Exception):
    """
    Exception associated with GFF3 data.
    """
    def __init__(self, message, fileName=None, lineNumber=None):
        """
        :param message: The error message to emit
        :param fileName: GFF3 file being processed
        :param lineNumber: line in GFF3 file, can be int or string
        """
        if fileName is not None:
            if lineNumber is not None:
                message = "{}:{}: {}".format(fileName, lineNumber, message)
            else:
                message = "{}: {}".format(fileName, message)
        super(GFF3Exception, self).__init__(message)


# characters forcing encode for columns
_encodeColReStr = "\t|\n|\r|%|[\x00-\x1F]|\x7f"
_encodeColRe = re.compile(_encodeColReStr)

# characters forcing encode for attribute names or values
_encodeAttrReStr = _encodeColReStr + "|;|=|&|,"
_encodeAttrRe = re.compile(_encodeAttrReStr)


def _encodeAttr(v):
    """
    Encode a attribute name or value if is has special characters.

    :param str v: The attribute string with possible special characters.
    :return str: Encoded attribute string with special characters escaped.
    """
    if _encodeAttrRe.search(v):
        return urllib.quote(v)
    else:
        return v


class Feature(object):
    """
    Feature as parsed from a GFF3.
    """
    def __init__(self, seqname, source, type, start, end, score, strand,
                 frame, attributes):
        """
        :param str seqname: name of reference sequence containing annotation
        :param str source: description of algorithm or procedure
            that generated the feature. Ex: ("HAVANA")
        :param str type: ontology name of feature type (ex: "gene")
        :param int start: starting genomic coordinate of feature
        :param int end: ending genomic coordinate of feature
        :param double score: quality score assigned to feature
        :param str strand: DNA strand that feature is on
            ("+"/"-", or None if neither)
        :param frame: corresponds to the "phase" column in GFF3
        :param attributes: a dict of lists/tuples. Missing attributes,
            coded in a GFF3 file as ".", are represented here as None.
        :return:
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
        self.parents = set()
        self.children = set()

    @staticmethod
    def _dotIfNone(val):
        return val if val is not None else '.'

    def _attributeStr(self, name):
        """
        Return name=value for a single attribute
        """
        return "{}={}".format(
            _encodeAttr(name),
            ",".join([_encodeAttr(v) for v in self.attributes[name]]))

    def _attributeStrs(self):
        """
        Return name=value, semi-colon-separated string for attributes,
        including url-style quoting
        """
        return ";".join([self._attributeStr(name)
                         for name in self.attributes.iterkeys()])

    def __str__(self):
        """
        Return the object as a valid GFF3 record line.
        """
        return "\t".join([self.seqname, self.source, self.type,
                          str(self.start), str(self.end),
                          self._dotIfNone(self.score),
                          self._dotIfNone(self.strand),
                          self._dotIfNone(self.frame),
                          self._attributeStrs()])

    @property
    def featureName(self):
        """
        ID attribute from GFF3 or None if record doesn't have it.
        Called "Name" rather than "Id" within GA4GH, as there is
        no guarantee of either uniqueness or existence.
        """
        featId = self.attributes.get("ID")
        if featId is not None:
            featId = featId[0]
        return featId

    @property
    def uniqueId(self):
        """
        Integer ID for this feature, guaranteed unique within a GFF3Set
        being parsed, suitable for use as a DB primary key.
        """
        return id(self)


class Gff3Set(object):
    """
    A set of GFF3 sequence annotations
    """
    def __init__(self, fileName=None):
        self.fileName = fileName
        self.roots = set()     # root nodes (those with out parents)
        # index of features by id. GFF3 allows disjoint features with
        # the same id.  None is used to store features without ids
        self.byFeatureName = collections.defaultdict(list)

    def add(self, feature):
        """
        Add a feature record by featureName (which may be None)

        :param feature: Feature object being added.
        """
        self.byFeatureName[feature.featureName].append(feature)

    def _linkFeature(self, feature):
        """
        Link a feature with its parents.
        """
        parentNames = feature.attributes.get("Parent")
        if parentNames is None:
            self.roots.add(feature)
        else:
            for parentName in parentNames:
                self._linkToParent(feature, parentName)

    def _linkToParent(self, feature, parentName):
        """
        Link a feature with its children
        """
        parentParts = self.byFeatureName.get(parentName)
        if parentParts is None:
            raise GFF3Exception(
                "Parent feature does not exist: {}".format(parentName),
                self.fileName)
        # parent maybe disjoint
        for parentPart in parentParts:
            feature.parents.add(parentPart)
            parentPart.children.add(feature)

    def linkChildFeaturesToParents(self):
        """
        finish loading the set, constructing the tree
        """
        # features maybe disjoint
        for featureParts in self.byFeatureName.itervalues():
            for feature in featureParts:
                self._linkFeature(feature)

    @staticmethod
    def _recSortKey(r):
        """
        Sort order for Features, by genomic coordinate,
        disambiguated by feature type (alphabetically).
        """
        return r.seqname, r.start, -r.end, r.type

    def _writeRec(self, fh, rec):
        """
        Writes a single record to a file provided by the filehandle fh.
        """
        fh.write(str(rec) + "\n")
        for child in sorted(rec.children, key=self._recSortKey):
            self._writeRec(fh, child)

    def write(self, fh):
        """
        Write set to a GFF3 format file.

        :param file fh: file handle for file to write to
        """
        fh.write(GFF3_HEADER+"\n")
        for root in sorted(self.roots, key=self._recSortKey):
            self._writeRec(fh, root)


class Gff3Parser(object):
    """
    Parses a GFF3 file into a Gff3Set. Performs basic validation,
    but does not fully test for conformance to GFF3 spec.
    """

    def __init__(self, fileName):
        """
        :param str fileName: Name of GFF3 file to parse,
            compressed files (.gz or .bz2) will be automatically decompressed.
        """
        self.fileName = fileName
        self.lineNumber = 0

    def _open(self):
        """
        open input file, optionally with decompression
        """
        if self.fileName.endswith(".gz"):
            return gzip.open(self.fileName)
        elif self.fileName.endswith(".bz2"):
            return bz2.BZ2File(self.fileName)
        else:
            return open(self.fileName)

    # parses `attr=val'; GFF3 spec is not very specific on the allowed values.
    SPLIT_ATTR_RE = re.compile("^([a-zA-Z][^=]*)=(.*)$")

    def _parseAttrVal(self, attrStr):
        """
        Returns tuple of tuple of (attr, value), multiple are returned to
        handle multi-value attributes.
        """
        m = self.SPLIT_ATTR_RE.match(attrStr)
        if m is None:
            raise GFF3Exception(
                "can't parse attribute/value: '" + attrStr +
                "'", self.fileName, self.lineNumber)
        name = urllib.unquote(m.group(1))
        val = m.group(2)
        # Split by comma to separate then unquote.
        # Commas in values must be url encoded.
        return name, [urllib.unquote(v) for v in val.split(',')]

    SPLIT_ATTR_COL_RE = re.compile("; *")

    def _parseAttrs(self, attrsStr):
        """
        Parse the attributes and values
        """
        attributes = dict()
        for attrStr in self.SPLIT_ATTR_COL_RE.split(attrsStr):
            name, vals = self._parseAttrVal(attrStr)
            if name in attributes:
                raise GFF3Exception(
                    "duplicated attribute name: {}".format(name),
                    self.fileName, self.lineNumber)
            attributes[name] = vals
        return attributes

    GFF3_NUM_COLS = 9

    def _parseRecord(self, gff3Set, line):
        """
        Parse one record.
        """
        row = line.split("\t")
        if len(row) != self.GFF3_NUM_COLS:
            raise GFF3Exception(
                "Wrong number of columns, expected {}, got {}".format(
                    self.GFF3_NUM_COLS, len(row)),
                self.fileName, self.lineNumber)
        feature = Feature(
            urllib.unquote(row[0]),
            urllib.unquote(row[1]),
            urllib.unquote(row[2]),
            int(row[3]), int(row[4]),
            row[5], row[6], row[7],
            self._parseAttrs(row[8]))
        gff3Set.add(feature)

    # spaces or comment line
    IGNORED_LINE_RE = re.compile("(^[ ]*$)|(^[ ]*#.*$)")

    @staticmethod
    def _isIgnoredLine(line):
        return Gff3Parser.IGNORED_LINE_RE.search(line) is not None

    def _checkHeader(self, line):
        # split to allow multiple spaces and tabs
        if line.split() != GFF3_HEADER.split():
            raise GFF3Exception(
                "First line is not GFF3 header ({}), got: {}".format(
                    GFF3_HEADER, line), self.fileName, self.lineNumber)

    def _parseLine(self, gff3Set, line):
        if self.lineNumber == 1:
            self._checkHeader(line)
        elif not self._isIgnoredLine(line):
            self._parseRecord(gff3Set, line)

    def parse(self):
        """
        Run the parse and return the resulting Gff3Set object.
        """
        fh = self._open()
        try:
            gff3Set = Gff3Set(self.fileName)
            for line in fh:
                self.lineNumber += 1
                self._parseLine(gff3Set, line[0:-1])
        finally:
            fh.close()
        gff3Set.linkChildFeaturesToParents()
        return gff3Set
