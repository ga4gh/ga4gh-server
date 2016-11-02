"""
Generate a SQLite DB file from a sequence annotation file in GFF3 format
(http://sequenceontology.org/resources/gff3.html)
The resulting database file will be mapped to a FeatureSet in the current
GA4GH server implementation, and each row in the 'feature' table will be
mapped to a Feature.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import sys
import json
import sqlite3

import ga4gh.common.utils as utils
import glue

glue.ga4ghImportGlue()
import ga4gh.server.gff3 as gff3  # NOQA

# TODO: Shift this to use the Gff3DbBackend class.

# The columns of the FEATURE table correspond to the columns of a GFF3,
# with three additional columns prepended representing the ID of this feature,
# the ID of its parent (if any), and a whitespace separated array
# of its child IDs.

_dbTableSQL = (
    "CREATE TABLE FEATURE( "
    "id INTEGER PRIMARY KEY NOT NULL, "
    "parent_id INTEGER, "
    "child_ids TEXT, "
    "reference_name TEXT, "
    "source TEXT, "
    "type TEXT, "
    "start INT, "
    "end INT, "
    "score REAL, "
    "strand TEXT, "
    "name TEXT,"
    "gene_name TEXT,"
    "transcript_name TEXT,"
    "attributes TEXT);")


def _db_serialize(pyData):
    return json.dumps(pyData, separators=(',', ':'))


class Gff32Db(object):
    """
    Represents a unit of work for this script: Parse a GFF3 file
    using an external GFF3 parser, and create a corresponding SQLite DB file
    by iterating through the resulting parsed dictionary object
    (gff3Data.byfeatureName).
    """
    def __init__(self, inputFile, outputFile):
        """
        :param inputFile: source GFF3 filename (can be a full path)
        :param outputFile: destination sqlite filename (ditto)
        """
        self.gff3File = inputFile
        self.dbFile = outputFile
        self.valueList = []
        self.batchSize = 100
        if os.path.exists(outputFile):
            print("DB output file already exists, please remove or rename.",
                  file=sys.stderr)
            exit()

    def _batchInsertValues(self, values, dbcur, dbconn):
        self.valueList.append(values)
        if len(self.valueList) >= self.batchSize:
            self._insertValues(dbcur, dbconn)

    def _insertValues(self, dbcur, dbconn):
        if len(self.valueList) > 0:
            sql = "INSERT INTO Feature VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
            dbcur.executemany(sql, self.valueList)
            dbconn.commit()
            self.valueList = []

    def run(self):
        gff3Data = gff3.Gff3Parser(self.gff3File).parse()
        dbconn = sqlite3.connect(self.dbFile)
        dbcur = dbconn.cursor()
        dbcur.execute(_dbTableSQL)  # create table
        dbcur.execute("PRAGMA journal_mode=WAL")
        dbconn.commit()

        for i, featureName in enumerate(gff3Data.byFeatureName.keys()):
            for feature in gff3Data.byFeatureName[featureName]:
                # Ignores any parent IDs besides the first one.
                parentIds = [parent.uniqueId for parent in feature.parents]
                parentId = parentIds[0] if len(parentIds) > 0 else ''
                # FIXME: No current code to ensure childIDs in correct order.
                childIds = [child.uniqueId for child in feature.children]
                values = (
                    feature.uniqueId,
                    parentId,
                    _db_serialize(childIds),
                    feature.seqname,
                    feature.source,
                    feature.type,
                    feature.start,
                    feature.end,
                    feature.score,
                    feature.strand,
                    feature.featureName,
                    feature.attributes.get("gene_name", [None])[0],
                    feature.attributes.get("transcript_name", [None])[0],
                    _db_serialize(feature.attributes))
                self._batchInsertValues(values, dbcur, dbconn)
        self._insertValues(dbcur, dbconn)
        dbcur.execute((
            "create INDEX idx1 "
            "on feature(start, end, reference_name)"))
        dbcur.execute("PRAGMA INDEX_LIST('feature')")

        dbcur.close()
        dbconn.close()


@utils.Timed()
def main():
    parser = argparse.ArgumentParser(
        description="Script to generate SQLite database corresponding to "
        "an input GFF3 sequence annotations file.")
    parser.add_argument(
        "--outputFile", "-o", default="annotations.db",
        help="The file to output the server-ready database to.")
    parser.add_argument(
        "--inputFile", "-i",
        help="Path to input GFF3 file.",
        default='.')
    parser.add_argument('--verbose', '-v', action='count', default=0)
    args = parser.parse_args()
    g2d = Gff32Db(args.inputFile, args.outputFile)
    g2d.run()


if __name__ == "__main__":
    main()
