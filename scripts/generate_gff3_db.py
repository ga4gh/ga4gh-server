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

import utils
import ga4gh.gff3Parser as gff3

# TODO: Shift this to use the Gff3DbBackend class.

# The columns of the FEATURE table correspond to the columns of a GFF3,
# with three additional columns prepended representing the ID of this feature,
# the ID of its parent (if any), and a whitespace separated array
# of its child IDs.
_dbTableSQL = (
                "CREATE TABLE FEATURE( "
                "id TEXT PRIMARY KEY NOT NULL, "
                "parent_id TEXT, "
                "child_ids TEXT, "
                "reference_name TEXT, "
                "source TEXT, "
                "ontology_term TEXT, "
                "start INT, "
                "end INT, "
                "score REAL, "
                "strand TEXT, "
                "attributes TEXT);"
                )


def _db_serialize(pyData):
    return json.dumps(pyData, separators=(',', ':'))


class Gff32Db(object):
    """
    Represents a unit of work for this script: Parse a GFF3 file
    using an external GFF3 parser, and create a corresponding SQLite DB file
    by iterating through the resulting parsed dictionary object
    (gff3Data.byFeatureId).
    """
    def __init__(self, args):
        self.gff3File = args.inputFile
        self.dbFile = args.outputFile
        if os.path.exists(args.outputFile):
            print("DB output file already exists, please remove or rename.",
                  file=sys.stderr)
            exit()

    def run(self):
        print("Running GFF3 parser...", file=sys.stderr)
        gff3Data = gff3.Gff3Parser(self.gff3File).parse()

        print("Generating database...", file=sys.stderr)
        dbconn = sqlite3.connect(self.dbFile)
        dbcur = dbconn.cursor()
        dbcur.execute(_dbTableSQL)  # create table
        for featureId in gff3Data.byFeatureId.keys():
            feature = gff3Data.byFeatureId[featureId][0]
            rowInsertSQL = "INSERT INTO FEATURE VALUES("
            rowInsertSQL += "'{}', ".format(featureId)

            # FIXME: Ignores any parent IDs besides the first one.
            parentIds = [parent.featureId for parent in feature.parents]
            if len(parentIds) == 0:
                rowInsertSQL += "'',"
            else:
                rowInsertSQL += "'{}', ".format(parentIds[0])

            # FIXME: No current way to get childIDs in correct order.
            childIds = [child.featureId for child in feature.children]
            rowInsertSQL += "'{}', ".format(_db_serialize(childIds))

            rowInsertSQL += "'{}', ".format(feature.seqname)
            rowInsertSQL += "'{}', ".format(feature.source)
            rowInsertSQL += "'{}', ".format(feature.type)
            rowInsertSQL += "'{}', ".format(feature.start)
            rowInsertSQL += "'{}', ".format(feature.end)
            rowInsertSQL += "'{}', ".format(feature.score)
            rowInsertSQL += "'{}', ".format(feature.strand)
            rowInsertSQL += "'{}');".format(_db_serialize(feature.attributes))

            # print(rowInsertSQL)  #DEBUG
            dbcur.execute(rowInsertSQL)
            dbconn.commit()
        dbcur.close()
        dbconn.close()
        print("Done.", file=sys.stderr)


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
    g2d = Gff32Db(args)
    g2d.run()


if __name__ == "__main__":
    main()
