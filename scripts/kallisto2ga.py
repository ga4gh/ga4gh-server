"""
    Script to parse the abundances.txt output file produced by kallisto.  As
    current version of the rnaQuant is using a flat file for the data this
    parser is living in server/scripts.  Once the backend store is ready for
    rnaQuant then this functionality will move to server/loaders.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
import errno
import optparse

import utils


def makeDir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


# TODO: placeholder values need to be calculated then removed
def getCount(expressionId):
    rawCount = 0
    return "{}".format(rawCount)


# TODO: placeholder values need to be calculated then removed
def getScore(expressionId):
    rawScore = 0.0

    return "{:0.2f}".format(rawScore)


def writeExpression(analysisId, annotationId, rnaQuantId, quantfile, rnaDB):
    # kallisto header
    # target_id	length	eff_length	est_counts	tpm
    isNormalized = True
    units = "TPM"
    # strip header and print - log it instead?
    print(quantfile.readline())
    for expression in quantfile:
        fields = expression.strip().split("\t")
        expressionLevel = fields[4]
        expressionId = fields[0]
        name = analysisId
        # TODO: properly handle feature group
        featureGroupId = expressionId
        rawCount = fields[3]
        score = 0.0

        datafields = (expressionId, name, rnaQuantId, annotationId, expressionLevel,
                      featureGroupId, isNormalized, rawCount,
                      str(score), units)
        rnaDB.addExpression(datafields)


def writeRnaseqTable(rnaDB, analysisIds, description, annotationId, readGroupId=None):
    if readGroupId is None:
        readGroupId = ""
    for analysisId in analysisIds:
        datafields = (analysisId, annotationId, description, analysisId, readGroupId)
        rnaDB.addRNAQuantification(datafields)


def writeExpressionTable(data, annotationId, rnaQuantId, rnaDB):
    for analysisId, quantfile in data:
        print("processing {}".format(analysisId))
        writeExpression(analysisId, annotationId, rnaQuantId, quantfile, rnaDB)


# TODO: not implemented
def writeRawBootstrap(data, annotationId, outputFolder, bootstrapId=None):
    pass


def makeParser(usage):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--readgroup", dest="readGroupId")
    parser.add_option("--annotation", dest="annotationId")
    parser.set_defaults(readGroupId=None, annotationId=None)
    return parser


# TODO: this is going to be the test of the sqlite backend.  Use this to parse and populate a
# dbfile and then see if it all connects properly
def main(argv):
    usage = "Usage: {0} <data-folder> <dbfile>".format(argv[0])
    if len(argv) < 3:
        print(usage)
        sys.exit(1)
    parser = makeParser(usage)
    (options, args) = parser.parse_args(argv[1:])
    description = "Kallisto demo data"
    # TOOD: remove default annotationId in favor of options.annotationId
    annotationId = "Homo_sapiens.GRCh38.rel79.cdna.all+ERCC"

    dataFolder = argv[1]
    sqlFilename = argv[2]
    rnaFolder = "rnaQuant"
    outputFolder = os.path.join(dataFolder, rnaFolder)
    rnaDB = utils.RNASqliteStore(outputFolder, sqlFilename)
    print("output folder:  {0}".format(outputFolder))
    rnaQuantId = "kallisto"
    # TODO: check to see if the rnaQuantId is in the db and exit if it is since this is a
    # generator and not an updater
    writeRnaseqTable(rnaDB, [rnaQuantId], description, annotationId,
                      readGroupId=options.readGroupId)
    analysisId = "kallisto"
    quantFile = open("abundance.txt", "r")
    writeExpressionTable([(analysisId, quantFile)], rnaQuantId, annotationId,
                          rnaDB)

    # TODO: add raw bootstrap stages

    print("DONE")


if __name__ == '__main__':
    main(sys.argv)
