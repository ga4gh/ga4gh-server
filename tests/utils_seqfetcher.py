"""
Test remote sequencing fetching as provided by ga4gh.utils.seqfetcher
"""

from __future__ import unicode_literals

import unittest

from ga4gh.utils.seqfetcher import fetch_seq


class TestSeqFetcher(unittest.TestCase):
    def test_ensembl_full(self):
        "Fetch full-length sequence from Ensembl"
        s = fetch_seq('ENSP00000288602')

        self.assertEqual(766,
                         len(s))

        self.assertEqual(s[0:10],
                         u'MAALSGGGGG')

    def test_ensembl_slice(self):
        "Fetch sequence slice from Ensembl"
        self.assertEqual(u'MAALSGGGGG',
                         fetch_seq('ENSP00000288602', 0, 10))

    def test_ncbi_full(self):
        "Fetch full-length sequence from NCBI"
        s = fetch_seq('NP_056374.2')

        self.assertEqual(1596,
                         len(s))

        self.assertEqual(s[0:10],
                         u'MESRETLSSS')

    def test_ncbi_slice(self):
        "Fetch sequence slice from NCBI"
        self.assertEqual(u'MESRETLSSS',
                         fetch_seq('NP_056374.2', 0, 10))

    def test_sequence_not_found(self):
        "Raise error for non-existent sequence"
        with self.assertRaises(RuntimeError):
            fetch_seq('NP_056374.299')

    def test_no_fetcher_for_accession(self):
        "Raise error for unrecognized accession"
        with self.assertRaises(RuntimeError):
            fetch_seq('ZZ_01234.5')
