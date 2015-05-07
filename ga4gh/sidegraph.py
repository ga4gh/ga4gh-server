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


def sqliteRows2dicts(sqliteRows):
    """
    unpacks sqlite rows as returned by fetchall
    into an array of simple dicts.
    """
    return map(lambda r: dict(zip(r.keys(), map(str, r))), sqliteRows)


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

    def searchReferenceSetsCount(self):
        """
        Returns the total number of results for a ReferenceSets search.
        """
        query = self._graphDb.execute("SELECT count(*) from ReferenceSet")
        count = int(query.fetchone()[0])
        return count

    def searchReferenceSets(self, limits=None):
        """
        Returns an array of dictionaries, each one encoding key-value
        information about a reference set.
        Limits, if provided, need to be a pair of ints: (start,end)
        end can be a value past the end of the resultset.
        """
        sql = """SELECT id FROM ReferenceSet ORDER BY ID"""
        if limits is not None and len(limits) > 1:
            start = int(limits[0])
            end = int(limits[1])
            sql += " LIMIT {}, {}".format(start, end)
        query = self._graphDb.execute(sql)
        return sqliteRows2dicts(query.fetchall())

    def searchVariantSetsCount(self):
        """
        Returns the total number of results for a variant set search.
        """
        query = self._graphDb.execute("SELECT count(*) from VariantSet")
        count = int(query.fetchone()[0])
        return count

    def searchVariantSets(self, limits=None):
        """
        Returns an array of dictionaries, each one encoding key-value
        information about a variant set.
        Limits, if provided, need to be a pair of ints: (start,end)
        end can be a value past the end of the resultset.
        """
        sql = """SELECT id FROM VariantSet ORDER BY ID"""
        if limits is not None and len(limits) > 1:
            start = int(limits[0])
            end = int(limits[1])
            sql += " LIMIT {}, {}".format(start, end)
        query = self._graphDb.execute(sql)
        return sqliteRows2dicts(query.fetchall())

    def searchSequencesCount(self):
        """
        Returns the total number of results for a given sequences search.
        """
        query = self._graphDb.execute("SELECT count(*) from Sequence")
        count = int(query.fetchone()[0])
        return count

    def searchSequences(self, limits=None):
        """
        Returns an array of dictionaries, each one encoding key-value
        information about a sequence.
        Limits, if provided, need to be a pair of ints: (start,end)
        end can be a value past the end of the resultset.
        """
        sql = """SELECT id, length FROM Sequence ORDER BY ID"""
        if limits is not None and len(limits) > 1:
            start = int(limits[0])
            end = int(limits[1])
            sql += " LIMIT {}, {}".format(start, end)
        query = self._graphDb.execute(sql)
        return sqliteRows2dicts(query.fetchall())

    def searchJoinsCount(self):
        """
        Returns the total number of results for a given joins search.
        """
        cursor = self._graphDb.cursor()
        query = cursor.execute("SELECT count(*) from GraphJoin")
        count = int(query.fetchone()[0])
        return count

    def searchJoins(self, limits=None):
        """
        Returns an array of dictionaries, each one encoding key-value
        information about a sequence.
        Limits, if provided, need to be a pair of ints: (start,end)
        end can be a value past the end of the resultset.
        """
        sql = """SELECT side1SequenceID, side1Position, side1StrandIsForward,
            side2SequenceID, side2Position, side2StrandIsForward
            FROM GraphJoin ORDER BY ID"""
        if limits is not None and len(limits) > 1:
            start = int(limits[0])
            end = int(limits[1])
            sql += " LIMIT {}, {}".format(start, end)
        query = self._graphDb.execute(sql)
        return sqliteRows2dicts(query.fetchall())

    def getSequenceBases(self, id, start=0, end=None):
        """
        Returns a string composed of bases in the sequence, between
        the start and end positions (or to end of string if end isn't
        specified).
        """
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
        returns variantSet ID for a given allele ID
        """
        query = self._graphDb.execute("""SELECT variantSetID
            FROM Allele where ID='{}'""".format(alleleID))
        return str(query.fetchone()[0])

    def getAllelePathItems(self, alleleID):
        """
        Returns n array of dicts with keys
        {sequenceID, start, length, strandIsForward}
        where isForward is a boolean indicating if segment should
        be read in the formward (from 5' to 3') direction.
        """
        sql = """SELECT sequenceID, start, length, strandIsForward
            FROM AllelePathItem WHERE alleleID = '{}'
            ORDER BY pathItemIndex
        """.format(alleleID)
        query = self._graphDb.execute(sql)
        return sqliteRows2dicts(query.fetchall())

    def getRefNames(self):
        """
        Return a list of all sequence names in the side graph,
        ordered by reference ID.
        """
        query = self._graphDb.execute("SELECT name FROM Reference ORDER BY ID")
        return list(map(lambda x: x[0], query.fetchall()))

    def getReference(self, name):
        """
        Returns a dictionary containing all known info about
        a particular reference. I'd make an object, but c'mon, this
        is simple key-value data, people.
        For now start=0 (TODO: Fix this as soon as SQL is fixed)
        """
        if name not in self.getRefNames():
            return None

        refSQL = """SELECT R.ID, R.name, S.sequenceRecordName, 0, S.length,
                S.md5checksum, R.isDerived, R.sourceDivergence,
                R.ncbiTaxonID
            FROM Reference R
            JOIN Sequence S ON R.sequenceID=S.ID
            WHERE R.name = '{}'
            ORDER BY R.ID
            LIMIT {}, {}
            """.format(name)
        r = {}
        (referenceId,
            r['name'],
            r['sequenceId'],
            r['start'],
            r['length'],
            r['md5checksum'],
            r['isDerived'],
            r['sourceDivergence'],
            r['ncbiTaxonId']) = self._graphDb.execute(refSQL).fetchone()

        # add the list of source accessions from a separate table
        accSQL = """SELECT accessionID FROM ReferenceAccession
            WHERE referenceID={}
            """.format(referenceId)
        acc = self._graphDb.execute(accSQL).fetchall()
        r['sourceAccessions'] = [] if len(acc) == 0\
            else [a[0] for a in acc]
        return r

    def getJoins(self, refName=None, start=0, end=None):
        """
        Return a list of all joins in the side graph adjacent to
        reference sequence 'refName' in the specified interval.
        If no refname is given, return all joins.
        If no valid interval is given,
        return all joins adjacent to that reference.
        Specifying an interval without a refName is invalid and
        should result in an error.
        Each join is represented by a tuple:
        (seq1Name, seq1Pos, seq1Strand, seq2Name, seq2Pos, seq2Strand)
        where strand is 'F' for forward, 'R' for reverse.
        """
        joinsSQL = """SELECT L.name, J.side1position,
                (CASE J.side1StrandIsForward WHEN 'TRUE'
                    THEN 'F' ELSE 'R' END),
                   R.name, J.side2position,
                (CASE J.side2StrandIsForward WHEN 'TRUE'
                    THEN 'F' ELSE 'R' END)
            FROM GraphJoin J
            JOIN Reference L ON J.side1SequenceID = L.sequenceID
            JOIN Reference R ON J.side2SequenceID = R.sequenceID
            """
        if refName is not None and refName in self.getRefNames():
            if type(end) is int:
                joinsSQL += """
                    fWHERE (L.name = '{0}'
                        AND J.side1Position >= {1}
                        AND J.side1Position <= {2})
                    OR (R.name = '{0}'
                        AND J.side2Position >= {1}
                        AND J.side2Position <= {2})
                    """.format(refName, start, end)
            else:  # only constrain start, which defaults to 0
                joinsSQL += """
                    WHERE (L.name = '{0}' AND J.side1Position >= {1})
                    OR (R.name = '{0}' AND J.side2Position >= {1})
                    """.format(refName, start)
        joinsCursor = self._graphDb.cursor()
        return joinsCursor.execute(joinsSQL).fetchall()

    def getSubgraph(self, refName=None, seedPosition=0, radius=0):
        """
        Returns the sepcified subgraph as a pair, (segments, joins)
        with segment represented by a triple (name, start, length)
        and each join a sextuple as returned by yieldJoins above.

        TODO: Current algorithm is recursive. This may overflow the stack
        with large subgraphs.
        """
        # recursive method: segments, joins are arrays of already
        # traversed elements. Assume half-open interval notation.
        # The refname, pos and rad are recomputed for each new segment
        # discovered by traversing an allowable join.
        # segments and joins are the arrays being built up,
        # and joinTaken, when not null, is the join just traversed
        # to arrive at the current position.
        segments = []
        joins = []

        def _getSubgraph(refName, pos, rad, segments, joins, joinTaken):
            self._logger.debug("at {}:{}~{} via {}".format(
                refName, pos, rad, joinTaken))
            if rad <= 0 or refName not in self.getRefNames():
                self._logger.debug("not a valid leaf")
                return
            if joinTaken is not None and joinTaken not in joins:
                joins.append(joinTaken)

            ref = self.getReference(refName)
            refStart = ref.get("start", 0)
            refLength = ref.get("length")
            segStart = max(refStart, pos - rad)
            segEnd = min(refStart + refLength, pos + rad)
            # the specified segment to explore may already be inside
            # another segment, or may partially intersect one or more
            # other segments. In the former case, just exit.
            # In the latter, figure out which existing segments to merge
            # into a larger segment with this one.
            intersectingSegs = []
            unionStart = segStart
            unionEnd = segEnd
            for i, (irefName, iStart, iLen) in enumerate(segments):
                iEnd = iStart + iLen
                if refName == irefName:
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
            segments.append((refName, unionStart, unionEnd - unionStart))
            # With segments adjusted, now explore joins...
            foundJoins = self.getJoins(refName, segStart, segEnd)
            for foundJoin in foundJoins:
                self._logger.debug("found join {}".format(foundJoin))
                (lref, lpos, lstr, rref, rpos, rstr) = foundJoin
                # make the recursive calls to follow all joins
                # encountered on the segment
                if lref == refName:
                    if lstr == 'R' and pos <= lpos <= segEnd:
                        _getSubgraph(
                            rref, rpos, rad - lpos + pos - 1, segments, joins,
                            foundJoin)
                    elif lstr == 'F' and segStart <= lpos <= pos:
                        _getSubgraph(
                            rref, rpos, rad - pos + lpos - 1, segments, joins,
                            foundJoin)
                if rref == refName:
                    if rstr == 'R' and pos <= rpos <= segEnd:
                        _getSubgraph(
                            lref, lpos, rad - rpos + pos - 1, segments, joins,
                            foundJoin)
                    elif rstr == 'F' and segStart <= rpos <= pos:
                        _getSubgraph(
                            lref, lpos, rad - pos + rpos - 1, segments, joins,
                            foundJoin)
            self._logger.debug("end {}:{}~{} via {}".format(
                refName, pos, rad, joinTaken))

        _getSubgraph(refName, seedPosition, radius, segments, joins, None)
        return (segments, joins)
