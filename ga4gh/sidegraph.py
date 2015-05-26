"""
Objects and methods encapsulating graph reference data whose graph
topology is stored in a Sqlite database and whose actual sequences
are found in an associated FASTA file.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import logging

import pyfasta
import sqlite3
import re

SIDEGRAPH_TRUE = 'TRUE'
SIDEGRAPH_FALSE = 'FALSE'

# Use the following regular expression to check
# if findBadChars.search(input) is not None:
# In which case, handle error - invalid characters found
findBadChars = re.compile('[\W]')


def _sqliteRows2dicts(sqliteRows):
    """
    unpacks sqlite rows as returned by fetchall
    into an array of simple dicts.
    """
    return map(lambda r: dict(zip(r.keys(), map(str, r))), sqliteRows)


def _limitsSql(limits):
    if limits is not None and len(limits) > 1:
        start = int(limits[0])
        end = int(limits[1])
        return " LIMIT {}, {}".format(start, end)
    else:
        return ""


def _whereClauseSql(**whereClauses):
    if whereClauses is not None:
        # wc is an array of "key ='value'" strings from whereClauses,
        # with all entries where value = None removed.
        wc = ["{} = '{}'".format(k, whereClauses[k])
              for k in whereClauses.keys()
              if whereClauses[k] is not None]
    if len(wc) > 0:
        return " WHERE " + " AND ".join(wc)
    else:
        return ""


class SideGraph(object):
    """
    Graph reference set based on a Sqlite file.

    Any named FASTA files without a fully qualified URL need to be
    located in the fastaDir.
    """

    def __init__(self, sqlFile, fastaDir="."):
        logging.basicConfig()
        self._logger = logging.getLogger(__name__)
        # NOTE: Set logging level to something reasonable
        # for semi-production
        self._sqlFile = sqlFile
        self._logger.setLevel(logging.DEBUG)
        self._fastaDir = fastaDir

    def __enter__(self):
        self._graphDb = sqlite3.connect(self._sqlFile)
        # row_factory setting is magic pixie dust to retrieve rows
        # as dictionaries. sqliteRows2dict relies on this.
        self._graphDb.row_factory = sqlite3.Row
        return self

    def __exit__(self, type, value, traceback):
        self._graphDb.close()

    def _countRows(self, tableName, **whereClauses):
        """
        Returns the total number of rows in the named table.
        """
        if findBadChars.search(tableName):
            raise Exception("table name contains invalid characters")

        sql = "SELECT count(*) from {}".format(tableName)
        sql += _whereClauseSql(**whereClauses)
        return int(self._graphDb.execute(sql).fetchone()[0])

    def _getRowsAsDicts(self, tableName, limits=None, **whereClauses):
        """
        Returns an array of dictionaries, each one encoding key-value
        information about each row of the requested table.

        Limits, if provided, need to be a pair of ints: (start,end)
        end can be a value past the end of the resultset.

        It is assumed that the table has an ID field that provides
        a canonical ordering on the table's rows.
        """
        if findBadChars.search(tableName):
            raise Exception("table name contains invalid characters")

        sql = "SELECT * FROM {}".format(tableName)
        sql += _whereClauseSql(**whereClauses)
        sql += " ORDER BY ID"
        sql += _limitsSql(limits)
        query = self._graphDb.execute(sql)
        return _sqliteRows2dicts(query.fetchall())

    def _getRowByIdAsDict(self, tableName, id):
        """
        Returns a single row of the table as a dictionary.

        It's assumed that the table has ID as its primary key.
        """
        if findBadChars.search(tableName):
            raise Exception("table name contains invalid characters")

        sql = """SELECT * FROM {} WHERE ID='{}'""".format(
            tableName, id)
        query = self._graphDb.execute(sql)
        return _sqliteRows2dicts(query.fetchall())[0]

    def searchReferenceSetsCount(self):
        return self._countRows("ReferenceSet")

    def searchReferenceSets(self, limits=None):
        return self._getRowsAsDicts("ReferenceSet", limits)

    def searchReferencesCount(self):
        return self._countRows("Reference")

    def searchReferences(self, limits=None):
        return self._getRowsAsDicts("Reference", limits)

    def searchVariantSetsCount(self):
        return self._countRows("VariantSet")

    def searchVariantSets(self, limits=None):
        return self._getRowsAsDicts("VariantSet", limits)

    def searchAllelesCount(self, variantSetId=None):
        return self._countRows("Allele",
                               variantSetID=variantSetId)

    def searchAlleles(self, limits=None, variantSetId=None):
        return self._getRowsAsDicts("Allele", limits,
                                    variantSetID=variantSetId)

    def searchCallSetsCount(self):
        return self._countRows("CallSet")

    def searchCallSets(self, limits=None):
        callSets = self._getRowsAsDicts("CallSet", limits)
        for cs in callSets:
            sql = """SELECT variantSetID FROM VariantSet_CallSet_Join
                WHERE callSetID='{}'""".format(cs["ID"])
            cs["variantSetIds"] = \
                [str(vs[0]) for vs in self._graphDb.execute(sql).fetchall()]
        return callSets

    def searchAlleleCallsCount(self, alleleId=None,
                               callSetId=None, variantSetId=None):
        return self._countRows("AlleleCall",
                               alleleID=alleleId, callSetID=callSetId)

    def searchAlleleCalls(self, limits=None,
                          alleleId=None, callSetId=None, variantSetId=None):
        # Can't use the regular _getRowsAsDicts mechanism as it has two
        # primary keys, not ID. Ordering is by those, lexicographic.
        # TODO: restrict by variantSet
        sql = "SELECT * FROM AlleleCall"
        sql += _whereClauseSql(alleleID=alleleId, callSetID=callSetId)
        sql += " ORDER BY alleleID, callSetID"
        sql += _limitsSql(limits)
        query = self._graphDb.execute(sql)
        return _sqliteRows2dicts(query.fetchall())

    def searchSequencesCount(self):
        return self._countRows("Sequence")

    def searchSequences(self, limits=None):
        return self._getRowsAsDicts("Sequence", limits)

    def searchJoinsCount(self):
        return self._countRows("GraphJoin")

    def searchJoins(self, limits=None):
        return self._getRowsAsDicts("GraphJoin", limits)

    def getSequenceBases(self, id, start=0, end=None):
        """
        Returns a string composed of bases in the sequence, between
        the start and end positions (or to end of string if end isn't
        specified).
        """
        if findBadChars.search(str(id)):
            raise Exception("ID contains invalid characters")

        sql = """SELECT F.fastaURI, S.sequenceRecordName
            FROM Sequence S
            JOIN FASTA F ON S.fastaID=F.ID
            WHERE S.ID='{}'""".format(id)
        query = self._graphDb.execute(sql)
        fetched = query.fetchone()
        if fetched is None:
            raise Exception("sequence id not found")

        fastaURI, recordName = fetched
        fastaFileName = os.path.join(self._fastaDir,
                                     os.path.basename(fastaURI))
        fasta = pyfasta.Fasta(
            fastaFileName, key_fn=lambda key: key.split()[0])
        return fasta[recordName][start:end]

    def getVariantSetIdForAllele(self, alleleID):
        """
        Returns variantSet ID for a given allele ID
        """
        allele = self._getRowByIdAsDict("Allele", alleleID)
        return allele["variantSetID"]

    def getAllelePathItems(self, alleleID):
        """
        Returns n array of dicts with keys
        {sequenceID, start, length, strandIsForward},
        where isForward is a boolean indicating if segment should
        be read in the formward (from 5' to 3') direction.
        """
        if findBadChars.search(alleleID):
            raise Exception("alleleID contains invalid characters")

        sql = """SELECT sequenceID, start, length, strandIsForward
            FROM AllelePathItem WHERE alleleID = '{}'
            ORDER BY pathItemIndex
        """.format(alleleID)
        query = self._graphDb.execute(sql)
        return _sqliteRows2dicts(query.fetchall())

    def getSequence(self, id):
        return self._getRowByIdAsDict("Sequence", id)

    def getJoins(self, seqId, start=0, end=None):
        """
        Return a list of all joins in the side graph adjacent to
        the requested sequence in the specified interval.

        If no valid interval is given,
        return all joins adjacent to the sequence.
        Joins are returned as a list of dicts, as usual.
        """
        if findBadChars.search(str(seqId)):
            raise Exception("input contains invalid characters")

        joinsSQL = "SELECT * from GraphJoin"
        if type(end) is int:
            joinsSQL += """
                WHERE (side1SequenceId = '{0}'
                    AND side1Position >= {1}
                    AND side1Position <= {2})
                OR (side2SequenceId = '{0}'
                    AND side2Position >= {1}
                    AND side2Position <= {2})
                """.format(seqId, int(start), int(end))
        else:  # only constrain start, which defaults to 0
            joinsSQL += """
                WHERE (side1SequenceId = '{0}' AND side1Position >= {1})
                OR (side2SequenceId = '{0}' AND side2Position >= {1})
                """.format(seqId, int(start))

        return _sqliteRows2dicts(self._graphDb.execute(joinsSQL).fetchall())

    def getSubgraph(self, seedSequenceId, seedPosition=0, radius=0,
                    referenceSetId=None):
        """
        Returns the sepcified subgraph as a pair, (segments, joins)
        with segment represented by a dictionary with keys
        (id, start, length) and each join as a typical join dictionary.

        Note: All segments returned are assumed forward strand.
        """
        # recursive method: segments, joins are arrays of already
        # traversed elements. Assume half-open interval notation.
        # The refname, pos and rad are recomputed for each new segment
        # discovered by traversing an allowable join.
        # segments and joins are the arrays being built up,
        # and joinTaken, when not null, is the join just traversed
        # to arrive at the current position.
        def _getSubgraph(seqId, pos, fwd, rad,
                         segments, joins,
                         joinTaken=None):
            """
            First three inputs describe current search "head":
              seqId - sequence id
              pos - position on sequence
              fwd - boolean: if true, walking forward, else reverse.
            Next two describe accumulated lists of segments and joins
              found so far.
            Final argument encodes which join was taken to arrive
              at current position, if any.

            Nothing is returned: it's a tail recursion, with the result
            being built up in the segments and joins arrays passed in.
            """
            self._logger.debug("at {}:{}{}{} via {}".format(
                seqId, pos, ">" if fwd else "<", rad, joinTaken))
            if rad <= 0:
                self._logger.debug("radius reached zero")
                return
            if joinTaken is not None and joinTaken not in joins:
                joins.append(joinTaken)

            # TODO - limit to reference set logic
            # There is nothing stopping a reference from sourcing several
            # noncontiguous segments of the same sequence, but it would
            # SERIOUSLY screw up any implementation of the below. I suggest
            # we disallow that by fiat.
            sq = self.getSequence(seqId)
            sqStart = 0
            sqLength = int(sq["length"])
            segStart = segEnd = pos
            if fwd:
                segEnd = min(sqLength, pos + rad)
            else:
                segStart = max(sqStart, pos - rad)

            # the specified segment to explore may already be inside
            # another segment, or may partially intersect one or more
            # other segments. In the former case, just exit.
            # In the latter, figure out which existing segments to merge
            # into a larger segment with this one.
            intersectingSegs = []
            unionStart = segStart
            unionEnd = segEnd
            for i, segDict in enumerate(segments):
                iSeqId = segDict["sequenceID"]
                iStart = segDict["start"]
                iLen = segDict["length"]
                iEnd = iStart + iLen
                if seqId == iSeqId:
                    if iStart <= segStart and iEnd >= segEnd:
                        self._logger.debug("nothing new here")
                        return
                    if segStart <= iStart <= segEnd or\
                       segStart <= iEnd <= segEnd:
                        intersectingSegs.append(i)
                        unionStart = min(unionStart, iStart)
                        unionEnd = max(unionEnd, iEnd)
            # reverse-iterate the intersecting segment indices and prune
            for i in intersectingSegs[::-1]:
                del segments[i]
            # replace them with the unified new segment
            segments.append(dict(
                sequenceID=seqId,
                start=unionStart,
                length=unionEnd - unionStart,
                strandIsForward=SIDEGRAPH_TRUE))
            # With segments adjusted, now explore joins...
            foundJoins = self.getJoins(seqId, segStart, segEnd)
            self._logger.debug("looking for joins on seqId {} {}-{}".format(
                seqId, segStart, segEnd))

            for foundJoin in foundJoins:
                self._logger.debug("found join {}".format(foundJoin))
                seq1 = foundJoin["side1SequenceID"]
                pos1 = int(foundJoin["side1Position"])
                fwd1 = foundJoin["side1StrandIsForward"] == SIDEGRAPH_TRUE

                seq2 = foundJoin["side2SequenceID"]
                pos2 = int(foundJoin["side2Position"])
                fwd2 = foundJoin["side2StrandIsForward"] == SIDEGRAPH_TRUE
                # make recursive calls to follow all joins
                # encountered on the segment
                if seq1 == seqId and segStart <= pos1 <= segEnd:
                    # check 1st side of join
                    if fwd and not fwd1:  # walking forward
                        _getSubgraph(
                            seq2, pos2, fwd2, rad - pos1 + pos - 1,
                            segments, joins, foundJoin)
                    elif not fwd and fwd1:  # walking reverse
                        _getSubgraph(
                            seq2, pos2, fwd2, rad - pos + pos1 - 1,
                            segments, joins, foundJoin)
                if seq2 == seqId and segStart <= pos2 <= segEnd:
                    # check 2nd side of join
                    if fwd and not fwd2:  # walking forward
                        _getSubgraph(
                            seq1, pos1, fwd1, rad - pos2 + pos - 1,
                            segments, joins, foundJoin)
                    elif not fwd and fwd2:  # walking reverse
                        _getSubgraph(
                            seq1, pos1, fwd1, rad - pos + pos2 - 1,
                            segments, joins, foundJoin)
            self._logger.debug("end {}:{}~{} via {}".format(
                seqId, pos, rad, joinTaken))

        segments = []
        joins = []
        # recursively fill out subgraph walking forward
        _getSubgraph(str(seedSequenceId), seedPosition, True,
                     radius, segments, joins, None)
        # and back
        _getSubgraph(str(seedSequenceId), seedPosition, False,
                     radius, segments, joins, None)

        return (segments, joins)
