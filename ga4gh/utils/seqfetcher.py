"""Fetch sequences and sequence slices from remote sources

This code was originally derived from the hgvs package [1] (Apache
licensed). The code was changed to remove caching, which is not
required for GA4GH server purposes, and to raise standard Python
exceptions for errors (rather than HGVSError subclasses).

IMPORTANT: This code uses interbase *intervals*, which are numerically
equivalent to 0-based, right-open nucleotide coordinates and familiar
to Python programmers. The _i suffix in call signatures emphasizes
this convention.

[1] source code at bitbucket: https://goo.gl/tO7GeW

"""

import logging
import re

import requests

logger = logging.getLogger(__name__)


def fetch_seq(ac, start_i=None, end_i=None):
    """Fetches sequences and subsequences from NCBI eutils and Ensembl
    REST interfaces.

    :param string ac: accession of sequence to fetch
    :param int start_i: start position of *interbase* interval
    :param int end_i: end position of *interbase* interval

    **IMPORTANT** start_i and end_i specify 0-based interbase
    coordinates, which refer to junctions between nucleotides.  This
    is numerically equivalent to 0-based, right-open nucleotide
    coordinates.

    Without an interval, the full sequence is returned::

    >>> len(fetch_seq('NP_056374.2'))
    1596

    Therefore, it's preferable to provide the interval rather than
    using Python slicing sequence on the delivered sequence::

    >>> fetch_seq('NP_056374.2',0,10)   # This!
    'MESRETLSSS'

    >>> fetch_seq('NP_056374.2')[0:10]  # Not this!
    'MESRETLSSS'

    >>> fetch_seq('NP_056374.2',0,10) == fetch_seq('NP_056374.2')[0:10]
    True

    Providing intervals is especially important for large sequences::

    >>> fetch_seq('NC_000001.10',2000000,2000030)
    'ATCACACGTGCAGGAACCCTTTTCCAAAGG'

    This call will pull back 30 bases plus overhead; without the
    interval, one would receive 250MB of chr1 plus overhead!

    Essentially any RefSeq, Genbank, BIC, or Ensembl sequence may be
    fetched:

    >>> [(ac,fetch_seq(ac,0,25))
    ... for ac in ['NG_032072.1', 'NW_003571030.1', 'NT_113901.1',
    ... 'NC_000001.10','NP_056374.2', 'GL000191.1', 'KB663603.1',
    ... 'ENST00000288602', 'ENSP00000288602']] # doctest: +NORMALIZE_WHITESPACE
    [('NG_032072.1', 'AAAATTAAATTAAAATAAATAAAAA'),
     ('NW_003571030.1', 'TTGTGTGTTAGGGTGCTCTAAGCAA'),
     ('NT_113901.1', 'GAATTCCTCGTTCACACAGTTTCTT'),
     ('NC_000001.10', 'NNNNNNNNNNNNNNNNNNNNNNNNN'),
     ('NP_056374.2', 'MESRETLSSSRQRGGESDFLPVSSA'),
     ('GL000191.1', 'GATCCACCTGCCTCAGCCTCCCAGA'),
     ('KB663603.1', 'TTTATTTATTTTAGATACTTATCTC'),
     ('ENST00000288602', u'CGCCTCCCTTCCCCCTCCCCGCCCG'),
     ('ENSP00000288602', u'MAALSGGGGGGAEPGQALFNGDMEP')]


    RuntimeError is thrown in the case of errors::

    >>> fetch_seq('NM_9.9')
    Traceback (most recent call last):
       ...
    RuntimeError: No sequence available for NM_9.9

    >>> fetch_seq('QQ01234')
    Traceback (most recent call last):
       ...
    RuntimeError: No sequence fetcher for QQ01234

    """

    ac_dispatch = [
        {'re': re.compile('^(?:AC|N[CGMPRTW])_|^[A-L]\w\d|^U\d'),
         'fetcher': _fetch_seq_ncbi},
        {'re': re.compile('^ENS[TP]\d+'),
         'fetcher': _fetch_seq_ensembl},
        ]

    eligible_fetchers = [dr['fetcher']
                         for dr in ac_dispatch
                         if dr['re'].match(ac)]

    if len(eligible_fetchers) == 0:
        raise RuntimeError("No sequence fetcher for {ac}".format(ac=ac))

    elif len(eligible_fetchers) == 1:
        fetcher = eligible_fetchers[0]

    elif len(eligible_fetchers) >= 1:
        logger.debug("Multiple sequence fetchers found for "
                     "{ac}; using first".format(ac=ac))
        fetcher = eligible_fetchers[0]

    logger.debug("fetching {ac} with {f}".format(
        ac=ac, f=fetcher))

    try:
        return fetcher(ac, start_i, end_i)
    except requests.HTTPError:
        raise RuntimeError("No sequence available for {ac}".format(ac=ac))


# ###########################################################################
# Internal functions

def _fetch_seq_ensembl(ac, start_i=None, end_i=None):
    """Fetch the specified sequence slice from Ensembl using the public
    REST interface.

    An interbase interval may be optionally provided with start_i and
    end_i. However, the Ensembl REST interface does not currently
    accept intervals, so the entire sequence is returned and sliced
    locally.

    >>> len(_fetch_seq_ensembl('ENSP00000288602'))
    766

    >>> _fetch_seq_ensembl('ENSP00000288602',0,10)
    u'MAALSGGGGG'

    >>> _fetch_seq_ensembl('ENSP00000288602')[0:10]
    u'MAALSGGGGG'

    >>> ac = 'ENSP00000288602'
    >>> _fetch_seq_ensembl(ac ,0, 10) == _fetch_seq_ensembl(ac)[0:10]
    True

    """

    url_fmt = "http://rest.ensembl.org/sequence/id/{ac}"
    url = url_fmt.format(ac=ac)
    r = requests.get(url, headers={"Content-Type": "application/json"})
    r.raise_for_status()
    seq = r.json()['seq']
    return seq if (start_i is None or end_i is None) else seq[start_i:end_i]


def _fetch_seq_ncbi(ac, start_i=None, end_i=None):
    """Fetch sequences from NCBI using the eutils interface.

    An interbase interval may be optionally provided with start_i and
    end_i. NCBI eutils will return just the requested subsequence,
    which might greatly reduce payload sizes (especially with
    chromosome-scale sequences).

    >>> len(_fetch_seq_ncbi('NP_056374.2'))
    1596

    Pass the desired interval rather than using Python's [] slice
    operator.

    >>> _fetch_seq_ncbi('NP_056374.2',0,10)
    'MESRETLSSS'

    >>> _fetch_seq_ncbi('NP_056374.2')[0:10]
    'MESRETLSSS'

    >>> _fetch_seq_ncbi('NP_056374.2',0,10) == _fetch_seq_ncbi('NP_056374.2')[0:10]
    True

    """

    url_fmt = ("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
               "db=nucleotide&id={ac}&rettype=fasta")
    if (start_i is None or end_i is None):
        url = url_fmt.format(ac=ac)
    else:
        url_fmt += "&seq_start={start}&seq_stop={stop}"
        url = url_fmt.format(ac=ac, start=start_i + 1, stop=end_i)
    resp = requests.get(url)
    resp.raise_for_status()
    return ''.join(resp.content.splitlines()[1:])


# So that I don't forget why I didn't use ebi too:
# $ curl 'http://www.ebi.ac.uk/ena/data/view/AM407889.1&display=fasta'
# >ENA|AM407889|AM407889.2 Medicago sativa partial mRNA ...
# AACGTATCACACTTCTTCTCCATTTCTTTTTCTTACATCTTCTCTCTACAAATTCATTTC
# Note that we requested .1, got .2.  Implicit behavior bites again.


if __name__ == "__main__":
    import doctest
    doctest.testmod()
