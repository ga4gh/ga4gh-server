#!/usr/bin/env python
"""
Example of client query
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.client.client as client


def runDemo():

    httpClient = client.HttpClient("http://localhost:8000")
    iterator = httpClient.search_variants(
        "WyIxa2ctcDMtc3Vic2V0IiwidnMiLCJtdm5jYWxsIl0",
        reference_name="1", start=45000, end=50000)
    for variant in iterator:
        print(
            variant.reference_name, variant.start, variant.end,
            variant.reference_bases, variant.alternate_bases, sep="\t")


if __name__ == '__main__':
    runDemo()
