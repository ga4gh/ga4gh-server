#!/usr/bin/env python
"""
script to download gene level RNA quantifications from the ENCODE DCC site
https://www.encodeproject.org and convert to format used by the GA4GH API
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
import json
import errno
import urllib2
import optparse

import utils


class Color(object):
    """
    Color printing in terminal
    """

    RED = '\033[91m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    NORMAL = '\033[0m'
    YELLOW = '\033[93m'

    @staticmethod
    def red(astr):
        return '{}{}{}'.format(Color.RED, str(astr), Color.NORMAL)

    @staticmethod
    def green(astr):
        return '{}{}{}'.format(Color.GREEN, str(astr), Color.NORMAL)

    @staticmethod
    def yellow(astr):
        return '{}{}{}'.format(Color.YELLOW, str(astr), Color.NORMAL)

    @staticmethod
    def blue(astr):
        return '{}{}{}'.format(Color.BLUE, str(astr), Color.NORMAL)


def log(message):
    print("== {}".format(Color.green(message)), file=sys.stderr)


def getFilesFromHost(data, host, outputType, request=False, subset=None):
    """
        analysisId is the ENCODE accession number for the quantification
        result file.  There may be multiple quantifications on a single ENCODE
        project or experiment page.
    """
    fileList = [f for f in data.get('files') if
                (f.get('output_type') == outputType)]
    if subset:
        fileList = fileList[:subset]
    for f in fileList:
        analysisId = f.get('accession')
        if not request:
            yield analysisId
        else:
            href = f.get('href')
            url = "{host}{file}".format(host=host, file=href)
            yield (analysisId, urllib2.urlopen(url))


def getDataFromHost(rnaDB, url, headers, host, outputType, outputFolder,
                    description, annotationId, subset=None,
                    readGroupId=None):
    req = urllib2.Request(url, headers=headers)
    try:
        response = urllib2.urlopen(req)
    except urllib2.URLError as e:
        if hasattr(e, 'reason'):
            print('We failed to reach a server.')
            print('Reason: {}'.format(e.reason))
        elif hasattr(e, 'code'):
            print("The server couldn't fulfill the request.")
            print('Error code: {}'.format(e.code))
    else:
        jsonData = json.load(response)
        # TODO: if data already exists don't write anything - this is a loader
        # not an updater
        makeDir(outputFolder)
        writeRnaseqTables(rnaDB, getFilesFromHost(jsonData, host, outputType,
                          subset=subset), description, annotationId,
                          outputFolder, readGroupId=readGroupId)
        writeGeneExpressionTables(rnaDB, getFilesFromHost(jsonData, host,
                                  outputType, subset=subset, request=True),
                                  annotationId, outputFolder)


def makeDir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


# TODO: placeholder values need to be calculated then removed
def getCount(expressionId):
    rawCount = 0
    return "{:d}".format(rawCount)


# TODO: placeholder values need to be calculated then removed
def getScore(expressionId):
    rawScore = 0.0

    return "{:0.2f}".format(rawScore)


def writeRNAQuant(rnaDB, analysisId, description, annotationId,
                  readGroupId=None):
    if readGroupId is None:
        readGroupId = ""
    datafields = (analysisId, annotationId, description, analysisId,
                  readGroupId)
    rnaDB.addRNAQuantification(datafields)


def writeGeneExpression(rnaDB, analysisId, annotationId, quantfile,
                        tool='RSEM'):
    # RSEM gene expression table header:
    #   gene_id transcript_id(s)    length  effective_length    expected_count
    #   TPM FPKM    pme_expected_count  pme_TPM pme_FPKM    TPM_ci_lower_bound
    #   TPM_ci_upper_bound  FPKM_ci_lower_bound FPKM_ci_upper_bound
    # TODO: placeholder values need to be calculated then removed
    isNormalized = "True"
    units = "TPM"
    # log expression file header
    log(quantfile.readline())
    for expression in quantfile.readlines():
        fields = expression.strip().split("\t")
        expressionLevel = fields[5]
        # TODO: generate expressionID
        expressionId = fields[0]
        name = fields[0]
        featureGroupId = fields[0]
        rawCount = fields[4]
        score = (float(fields[10]) + float(fields[11]))/2

        datafields = (expressionId, name, analysisId, annotationId,
                      expressionLevel, featureGroupId, isNormalized, rawCount,
                      score, units)
        rnaDB.addExpression(datafields)


def writeRnaseqTables(rnaDB, analysisIds, description, annotationId,
                      outputFolder, readGroupId=None):
    log("Writing rnaseq tables")
    for analysisId in analysisIds:
        writeRNAQuant(rnaDB, analysisId, description, annotationId,
                      readGroupId=readGroupId)


def writeGeneExpressionTables(rnaDB, data, annotationId, outputFolder):
    log("Writing gene expression tables")
    for analysisId, quantfile in data:
        print("processing {}".format(analysisId))
        writeGeneExpression(rnaDB, analysisId, annotationId, quantfile)


def makeParser(usage):
    """
    Default values will download and parse 4 quantifications from the
    ENCODE Evaluation dataset as an example set of test data.
    """

    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--dataset", dest="dataset")
    parser.add_option("--description", dest="description")
    parser.add_option("--annotationId", dest="annotationId")
    parser.add_option("--file_limit", type="int", dest="subset")
    parser.add_option("--readGroupId", dest="readGroupId")
    parser.add_option("--expressionType", dest="expressionType")
    parser.set_defaults(dataset="ENCSR000AJW", annotationId="Gencodev16",
                        description="RNAseq data from ENCODE evaluation",
                        subset=4, readGroupId=None, expressionType="gene")
    return parser


def main(argv):

    usage = "Usage: {} <data-folder> <db-file>".format(argv[0])
    if len(argv) < 2:
        print(Color.red(usage))
        sys.exit(1)
    parser = makeParser(usage)
    (options, args) = parser.parse_args(argv[1:])
    host = "https://www.encodeproject.org"
    dataset = options.dataset
    url = "{host}/datasets/{dataset}/?frame=embedded".format(host=host,
                                                             dataset=dataset)
    headers = {
        'Accept': 'application/json; charset=utf-8'
    }
    if options.expressionType == "gene":
        dataType = "gene quantifications"
    dataFolder = argv[1]
    sqlFilename = argv[2]
    rnaFolder = "rnaQuant"
    outputFolder = os.path.join(dataFolder, rnaFolder)
    rnaDB = utils.RNASqliteStore(outputFolder, sqlFilename)
    subset = options.subset
    log("Downloading GA4GH test dataset - RNA Quantification API")
    print("ENCODE dataset: {}".format(Color.blue(dataset)))
    print("data type:      {}".format(Color.blue(dataType)))
    print("subset size:    {}".format(Color.blue(subset)))
    print("readGroupId  :  {}".format(Color.blue(options.readGroupId)))
    print("output folder:  {}".format(Color.blue(outputFolder)))
    # TODO: annotation name can be found at DCC...
    getDataFromHost(rnaDB, url, headers, host, dataType, outputFolder,
                    options.description, options.annotationId, subset,
                    options.readGroupId)
    log("DONE")


if __name__ == '__main__':
    main(sys.argv)
