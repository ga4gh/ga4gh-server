#!/usr/bin/env python
"""
Example of client query
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import utils
utils.ga4ghImportGlue()
import ga4gh.client as client  # NOQA


def runDemo():

    httpClient = client.HttpClient("http://localhost:8000")
    iterator = httpClient.search_variants(
        variantSetId="MWtnLXAzLXN1YnNldDptdm5jYWxs",
        referenceName="1", start=45000, end=50000)
    for variant in iterator:
        print(
            variant.referenceName, variant.start, variant.end,
            variant.referenceBases, variant.alternateBases, sep="\t")


if __name__ == '__main__':
    runDemo()
