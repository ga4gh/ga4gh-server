"""
Utilities for scripts
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import functools
import os
import shlex
import subprocess
import sys
import time

import humanize
import requests
import yaml
import pysam


def log(message):
    print(message)


class Timed(object):
    """
    Decorator that times a method, reporting runtime at finish
    """
    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            self.start = time.time()
            result = func(*args, **kwargs)
            self.end = time.time()
            self._report()
            return result
        return wrapper

    def _report(self):
        delta = self.end - self.start
        timeString = humanize.time.naturaldelta(delta)
        log("Finished in {} ({} seconds)".format(timeString, delta))


class FileDownloader(object):
    """
    Provides a wget-like file download and terminal display
    """
    defaultChunkSize = 1048576  # 1MB
    defaultStream = sys.stdout

    def __init__(self, url, path, chunkSize=defaultChunkSize,
                 stream=defaultStream):
        self.url = url
        self.path = path
        self.basename = os.path.basename(url)
        self.basenameLength = len(self.basename)
        self.chunkSize = chunkSize
        self.stream = stream
        self.bytesWritten = 0
        self.displayIndex = 0
        self.displayWindowSize = 20

    def download(self):
        self.stream.write("Downloading '{}' to '{}'\n".format(
            self.url, self.path))
        response = requests.get(self.url, stream=True)
        response.raise_for_status()
        self.contentLength = int(response.headers['content-length'])
        with open(self.path, 'wb') as outputFile:
            for chunk in response.iter_content(chunk_size=self.chunkSize):
                self.bytesWritten += self.chunkSize
                self._updateDisplay()
                outputFile.write(chunk)
        self.stream.write("\n")
        self.stream.flush()

    def _getFileNameDisplayString(self):
        if self.basenameLength <= self.displayWindowSize:
            return self.basename
        else:
            return self.basename  # TODO scrolling window here

    def _updateDisplay(self):
        fileName = self._getFileNameDisplayString()

        # TODO contentLength seems to slightly under-report how many bytes
        # we have to download... hence the min functions
        percentage = min(self.bytesWritten / self.contentLength, 1)
        numerator = humanize.filesize.naturalsize(
            min(self.bytesWritten, self.contentLength))
        denominator = humanize.filesize.naturalsize(
            self.contentLength)

        displayString = "{}   {:<6.2%} ({:>9} / {:<9})\r"
        self.stream.write(displayString.format(
            fileName, percentage, numerator, denominator))
        self.stream.flush()


def runCommandSplits(splits, silent=False):
    """
    Run a shell command given the command's parsed command line
    """
    if silent:
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(splits, stdout=devnull, stderr=devnull)
    else:
        subprocess.check_call(splits)


def runCommand(command, silent=False):
    """
    Run a shell command
    """
    splits = shlex.split(command)
    runCommandSplits(splits, silent=silent)


def getAuthValues(filePath='scripts/auth.yml'):
    """
    Return the script authentication file as a dictionary
    """
    return getYamlDocument(filePath)


def getYamlDocument(filePath):
    """
    Return a yaml file's contents as a dictionary
    """
    with open(filePath) as stream:
        doc = yaml.load(stream)
        return doc


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
        log("Creating alignment file '{}'".format(outputFilePath))
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
            log("Creating index file '{}'".format(indexFilePath))
            pysam.index(outputFilePath)
