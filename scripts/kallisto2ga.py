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


def writeRNAQuant(outfile, analysisId, description, annotationId,
                  readGroupId=None):
    if readGroupId is None:
        readGroupId = ""
    outline = "\t".join([analysisId, annotationId, description, analysisId,
                         readGroupId])
    print(outline, file=outfile)


def writeExpression(analysisId, annotationId, quantfile, quantOutfile):
    # kallisto header
    # target_id	length	eff_length	est_counts	tpm
    isNormalized = "True"
    units = "TPM"
    # strip header and print - log it instead?
    print(quantfile.readline())
    for expression in quantfile.readlines():
        fields = expression.strip().split("\t")
        expressionLevel = fields[4]
        expressionId = fields[0]
        # TODO: properly handle feature group
        featureGroupId = expressionId
        rawCount = fields[3]
        score = 0.0
        outline = "\t".join([expressionId, annotationId, expressionLevel,
                             featureGroupId, isNormalized, rawCount,
                             str(score), units])
        print(outline, file=quantOutfile)


def writeRnaseqTables(analysisIds, description, annotationId, outputFolder):
    for analysisId in analysisIds:
        makeDir(os.path.join(outputFolder, analysisId))
        rnaSeqTable = os.path.join(outputFolder, analysisId, "rnaseq.table")
        with open(rnaSeqTable, "w") as rnaQuantFile:
            writeRNAQuant(rnaQuantFile, analysisId, description, annotationId)


def writeExpressionTables(data, annotationId, outputFolder):
    for analysisId, quantfile in data:
        expTable = os.path.join(outputFolder, analysisId, "expression.table")
        print("processing {}".format(analysisId))
        with open(expTable, "w") as quantOutfile:
            writeExpression(analysisId, annotationId, quantfile, quantOutfile)


# TODO: not implemented - this is just copied from writeExpressionTables
def writeRawBootstrap(data, annotationId, outputFolder, bootstrapId=None):
    for (analysisId, quantfile) in data:
        # output table
        expTable = os.path.join(outputFolder, analysisId, "expression.table")

        # write expression table
        print("processing {}".format(analysisId))
        with open(expTable, "w") as quantOutfile:
            writeExpression(analysisId, annotationId, quantfile, quantOutfile,
                            featureGroupId=bootstrapId)


def makeParser(usage):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--readgroup", dest="readGroupId")
    parser.set_defaults(dataset=None)
    return parser


def main(argv):

    if len(argv) < 2:
        usage = "Usage: {0} <data-folder>".format(argv[0])
        print(usage)
        sys.exit(1)
    parser = makeParser(usage)
    (options, args) = parser.parse_args(argv[1:])
    description = "Kallisto demo data"
    annotationId = "Homo_sapiens.GRCh38.rel79.cdna.all+ERCC"

    dataFolder = argv[1]
    rnaFolder = "rnaQuant"
    outputFolder = os.path.join(dataFolder, rnaFolder)
    print("output folder:  {0}".format(outputFolder))
    makeDir(outputFolder)
    writeRnaseqTables(["kallisto"], description, annotationId, outputFolder,
                      readgroupId=options.readGroupId)
    analysisId = "kallisto"
    quantFile = open("abundance.txt", "r")
    writeExpressionTables([(analysisId, quantFile)], annotationId,
                          outputFolder)

    # TODO: add raw bootstrap stages

    print("DONE")


if __name__ == '__main__':
    main(sys.argv)
