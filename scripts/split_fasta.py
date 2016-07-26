"""
Split a multiple-sequence FASTA file into many single-sequence FASTA files
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os

import utils


def parseArgs():
    description = ("Split a multiple-sequence FASTA file into many "
                   "single-sequence FASTA files")
    parser = argparse.ArgumentParser(
        description=description)
    parser.add_argument("fastaFile", help="the FASTA file to split")
    args = parser.parse_args()
    return args


def decompressFasta(args):
    fastaFileName = args.fastaFile
    filename, extension = os.path.splitext(fastaFileName)
    if extension == '.gz':
        utils.log("Decompressing {}".format(fastaFileName))
        cmd = "gunzip {}".format(fastaFileName)
        utils.runCommand(cmd)
        fastaFileName = filename
    return fastaFileName


def splitFasta(fastaFileName):
    splitFileNames = []
    with open(fastaFileName) as fastaFile:
        currentFile = None
        currentFileName = None
        for line in fastaFile:
            if line[0] == '>':
                if currentFile is not None:
                    currentFile.close()
                currentFileName = line[1:].split()[0].strip() + '.fa'
                utils.log("Creating {}".format(currentFileName))
                splitFileNames.append(currentFileName)
                currentFile = open(currentFileName, 'w')
            currentFile.write(line)
        currentFile.close()
    return splitFileNames


def compressSplits(splitFileNames):
    compressedFileNames = []
    for splitFileName in splitFileNames:
        utils.log("Compressing {}".format(splitFileName))
        cmd = "bgzip {}".format(splitFileName)
        utils.runCommand(cmd)
        compressedFileName = "{}.gz".format(splitFileName)
        compressedFileNames.append(compressedFileName)
    return compressedFileNames


def indexSplits(compressedFileNames):
    for compressedFileName in compressedFileNames:
        utils.log("Indexing {}".format(compressedFileName))
        cmd = "samtools faidx {}".format(compressedFileName)
        utils.runCommand(cmd)


@utils.Timed()
def main():
    requiredExecutables = ['gunzip', 'bgzip', 'samtools']
    utils.requireExecutables(requiredExecutables)
    args = parseArgs()
    fastaFileName = decompressFasta(args)
    splitFileNames = splitFasta(fastaFileName)
    compressedFileNames = compressSplits(splitFileNames)
    indexSplits(compressedFileNames)


if __name__ == '__main__':
    main()
