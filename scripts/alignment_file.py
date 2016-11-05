"""
Tools for alignment files
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import pysam

import ga4gh.common.utils as utils


class AlignmentFileConstants(object):
    """
    A container class for constants dealing with alignment files
    """
    SAM = "SAM"
    BAM = "BAM"
    BAI = "BAI"


class AlignmentFileTool(object):
    """
    Helps with operations on BAM and SAM files
    """
    def __init__(self, inputFileFormat, outputFileFormat):
        self.inputFileFormat = inputFileFormat
        self.outputFileFormat = outputFileFormat
        self.args = None

    def parseArgs(self):
        description = "{} to {} conversion tool".format(
            self.inputFileFormat, self.outputFileFormat)
        parser = argparse.ArgumentParser(
            description=description)
        inputHelpText = "the name of the {} file to read".format(
            self.inputFileFormat)
        parser.add_argument(
            "inputFile", help=inputHelpText)
        outputHelpText = "the name of the {} file to write".format(
            self.outputFileFormat)
        defaultOutputFilePath = "out.{}".format(
            self.outputFileFormat.lower())
        parser.add_argument(
            "--outputFile", "-o", default=defaultOutputFilePath,
            help=outputHelpText)
        parser.add_argument(
            "--numLines", "-n", default=10,
            help="the number of lines to write")
        parser.add_argument(
            "--skipIndexing", default=False, action='store_true',
            help="don't create an index file")
        args = parser.parse_args()
        self.args = args

    def convert(self):
        # set flags
        if self.inputFileFormat == AlignmentFileConstants.SAM:
            inputFlags = "r"
        elif self.inputFileFormat == AlignmentFileConstants.BAM:
            inputFlags = "rb"
        if self.outputFileFormat == AlignmentFileConstants.SAM:
            outputFlags = "wh"
        elif self.outputFileFormat == AlignmentFileConstants.BAM:
            outputFlags = "wb"
        # open files
        inputFile = pysam.AlignmentFile(
            self.args.inputFile, inputFlags)
        outputFile = pysam.AlignmentFile(
            self.args.outputFile, outputFlags, header=inputFile.header)
        outputFilePath = outputFile.filename
        utils.log("Creating alignment file '{}'".format(outputFilePath))
        # write new file
        for _ in xrange(self.args.numLines):
            alignedSegment = inputFile.next()
            outputFile.write(alignedSegment)
        # clean up
        inputFile.close()
        outputFile.close()
        # create index file
        if (not self.args.skipIndexing and
                self.outputFileFormat == AlignmentFileConstants.BAM):
            indexFilePath = "{}.{}".format(
                outputFilePath, AlignmentFileConstants.BAI.lower())
            utils.log("Creating index file '{}'".format(indexFilePath))
            pysam.index(outputFilePath)
