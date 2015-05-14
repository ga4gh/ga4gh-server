"""
Unit tests for sidegraph.py, my toy side graph library.
Uses local data for now.
TODO: Generalize to take any old SQL/FASTA encoded sidegraph.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import sqlite3
import os
import glob

import ga4gh.sidegraph as sidegraph


class TestCase(unittest.TestCase):

    def setUp(self):
        self._files = "tests/data/graphs"
        self._db = os.path.join(self._files, "graph.db")

        self._conn = sqlite3.connect(self._db)

    def tearDown(self):
        self._conn.close()

    def testSearchSequencesCount(self):
        with sidegraph.SideGraph(self._db, self._files) as sg:
            self.assertEquals(sg.searchSequencesCount(), 5)

    def testSearchSequences(self):
        expected = [{'length': '1000',
                'sequenceRecordName': 'chr1',
                'ID': '1',
                'md5checksum': 'f28788f2203d71311b5d8fe81fe1bd2e',
                'fastaID': '1'},
            {'length': '1',
                'sequenceRecordName': 'chr1snp1',
                'ID': '2',
                'md5checksum': '7fc56270e7a70fa81a5935b72eacbe29',
                'fastaID': '1'},
            {'length': '1000',
                'sequenceRecordName': 'chr2',
                'ID': '3',
                'md5checksum': '2895d94c7491966ef9df7af5ecf77e9f',
                'fastaID': '1'},
            {'length': '1000',
                'sequenceRecordName': 'chr3',
                'ID': '5',
                'md5checksum': '52ae3ef016c60c5e978306e8d3334cd8',
                'fastaID': '1'},
            {'length': '30',
                'sequenceRecordName': 'chr2ins1',
                'ID': '8',
                'md5checksum': '4833b4fa1627b1ee25f83698f768f997',
                'fastaID': '1'}]
        with sidegraph.SideGraph(self._db, self._files) as sg:
            self.assertEquals(sg.searchSequences(), expected)

    def testSearchJoinsCount(self):
        with sidegraph.SideGraph(self._db, self._files) as sg:
            self.assertEquals(sg.searchJoinsCount(), 8)

    def testSearchJoins(self):
        expected = [{'side2Position': '0',
                'side2StrandIsForward': 'TRUE',
                'side2SequenceID': '2',
                'side1SequenceID': '1',
                'ID': '1',
                'side1Position': '79',
                'side1StrandIsForward': 'FALSE'},
            {'side2Position': '0',
                'side2StrandIsForward': 'FALSE',
                'side2SequenceID': '2',
                'side1SequenceID': '1',
                'ID': '2',
                'side1Position': '81',
                'side1StrandIsForward': 'TRUE'},
            {'side2Position': '85',
                'side2StrandIsForward': 'FALSE',
                'side2SequenceID': '1',
                'side1SequenceID': '1',
                'ID': '3',
                'side1Position': '72',
                'side1StrandIsForward': 'TRUE'},
            {'side2Position': '108',
                'side2StrandIsForward': 'TRUE',
                'side2SequenceID': '3',
                'side1SequenceID': '3',
                'ID': '4',
                'side1Position': '33',
                'side1StrandIsForward': 'FALSE'},
            {'side2Position': '310',
                'side2StrandIsForward': 'FALSE',
                'side2SequenceID': '3',
                'side1SequenceID': '1',
                'ID': '5',
                'side1Position': '233',
                'side1StrandIsForward': 'TRUE'},
            {'side2Position': '311',
                'side2StrandIsForward': 'TRUE',
               'side2SequenceID': '3',
                'side1SequenceID': '1',
                'ID': '6',
               'side1Position': '289',
                'side1StrandIsForward': 'FALSE'},
            {'side2Position': '0',
                'side2StrandIsForward': 'TRUE',
                'side2SequenceID': '8',
                'side1SequenceID': '3',
                'ID': '7',
             'side1Position': '300',
             'side1StrandIsForward': 'FALSE'},
            {'side2Position': '29',
                'side2StrandIsForward': 'FALSE',
                'side2SequenceID': '8',
                'side1SequenceID': '3',
                'ID': '8',
                'side1Position': '301',
                'side1StrandIsForward': 'TRUE'}]
        with sidegraph.SideGraph(self._db, self._files) as sg:
            self.assertEquals(sg.searchJoins(), expected)

    def testGetSequenceBases(self):
        expected = 'CAGTAGTAACCATAAACTTACGCTGGGGCT'
        with sidegraph.SideGraph(self._db, self._files) as sg:
            self.assertEquals(sg.getSequenceBases(8), expected)


    @unittest.skip("not yet implemented")
    def getVariantSetIdForAllele(self):
        pass

    def testGetJoins(self):
        expected = [
            {'side2Position': '0',
                'side2StrandIsForward': 'TRUE',
                'side2SequenceID': '2',
                'side1SequenceID': '1', 'ID': '1',
                'side1Position': '79',
                'side1StrandIsForward': 'FALSE'},
            {'side2Position': '0',
                'side2StrandIsForward': 'FALSE',
                'side2SequenceID': '2',
                'side1SequenceID': '1',
                'ID': '2',
                'side1Position': '81',
                'side1StrandIsForward': 'TRUE'}]
        with sidegraph.SideGraph(self._db, self._files) as sg:
            self.assertEquals(sg.getJoins(2), expected)

    def testGetSubgraph(self):
        expected = ([('2', 0, 1), ('1', 65, 20)],
                    [{'side2Position': '0',
                      'side2StrandIsForward': 'TRUE',
                      'side2SequenceID': '2',
                      'side1SequenceID': '1', 'ID': '1',
                      'side1Position': '79',
                      'side1StrandIsForward': 'FALSE'},
                     {'side2Position': '0',
                        'side2StrandIsForward': 'FALSE',
                        'side2SequenceID': '2',
                        'side1SequenceID': '1',
                        'ID': '2',
                        'side1Position': '81',
                        'side1StrandIsForward': 'TRUE'},
                     {'side2Position': '85',
                        'side2StrandIsForward': 'FALSE',
                        'side2SequenceID': '1',
                        'side1SequenceID': '1',
                        'ID': '3',
                        'side1Position': '72',
                        'side1StrandIsForward': 'TRUE'}])
        with sidegraph.SideGraph(self._db, self._files) as sg:
            self.assertEquals(sg.getSubgraph(1, 75, 10), expected)


if __name__ == '__main__':
    unittest.main()
