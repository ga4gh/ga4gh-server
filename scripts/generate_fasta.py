"""
Generate a random FASTA file
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import math
import random

import utils


firstLine = ">Generated file\n"
fileName = "generated.fa"


def parseArgs():
    """
    Parse the command line args
    """
    parser = argparse.ArgumentParser(
        description="Generate random FASTA files")
    basesDefault = 1000
    parser.add_argument(
        "--num-bases", "-n", default=basesDefault,
        help="number of bases to include; default {}".format(basesDefault))
    args = parser.parse_args()
    return args


def writeFasta(args):
    """
    Write the random fasta file
    """
    numBases = args.num_bases
    utils.log("writing {} bases to {} ...".format(numBases, fileName))
    with open(fileName, 'w') as fastaFile:
        fastaFile.write(firstLine)
        basesPerLine = 70
        numLines = int(math.ceil(numBases / basesPerLine))
        baseChoices = ['A', 'G', 'C', 'T']
        basesRemaining = numBases
        for i in range(numLines):
            if basesRemaining < basesPerLine:
                basesToWrite = basesRemaining
            else:
                basesToWrite = basesPerLine
            bases = ''.join(
                [random.choice(baseChoices) for _ in range(basesToWrite)])
            line = "{}\n".format(bases)
            fastaFile.write(line)
            basesRemaining -= basesToWrite
        assert basesRemaining == 0


def zipFasta():
    """
    Compress the fasta file
    """
    utils.log("zipping {} ...".format(fileName))
    cmd = "bgzip {}".format(fileName)
    utils.runCommand(cmd)


def indexFasta():
    """
    Create index on the fasta file
    """
    zipFileName = "{}.gz".format(fileName)
    utils.log("indexing {} ...".format(zipFileName))
    cmd = "samtools faidx {}".format(zipFileName)
    utils.runCommand(cmd)


def main():
    args = parseArgs()
    writeFasta(args)
    zipFasta()
    indexFasta()


if __name__ == '__main__':
    main()
