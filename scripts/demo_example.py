#!/usr/bin/env python
"""
Example of client query
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import ga4gh.client as client
import ga4gh.protocol as protocol


def runDemo():
    httpClient = client.HttpClient("http://server:80/current")
    request = protocol.SearchVariantsRequest()
    request.variant_set_ids = ["1kg-phase1"]
    request.reference_name = "2"
    request.start = 33100
    request.end = 34000
    for variant in httpClient.searchVariants(request):
        print(
            variant.reference_name, variant.start, variant.end,
            variant.reference_bases, variant.alternate_bases,
            sep="\t")


if __name__ == '__main__':
    runDemo()
