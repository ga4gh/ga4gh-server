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

import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()


def getPathOfExecutable(executable):
    """
    Returns the full path of the executable, or None if the executable
    can not be found.
    """
    exe_paths = os.environ['PATH'].split(':')
    for exe_path in exe_paths:
        exe_file = os.path.join(exe_path, executable)
        if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
            return exe_file
    return None


def requireExecutables(executables):
    """
    Check that all of the given executables are on the path.
    If at least one of them is not, exit the script and inform
    the user of the missing requirement(s).
    """
    missingExecutables = []
    for executable in executables:
        if getPathOfExecutable(executable) is None:
            missingExecutables.append(executable)
    if len(missingExecutables) > 0:
        log("In order to run this script, the following "
            "executables need to be on the path:")
        for missingExecutable in missingExecutables:
            print(missingExecutable)
        exit(1)


def ga4ghImportGlue():
    """
    Call this method before importing a ga4gh module in the scripts dir.
    Otherwise, you will be using the installed package instead of
    the development package.
    Assumes a certain directory structure.
    """
    path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(path)


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
        log("Finished in {} ({:.2f} seconds)".format(timeString, delta))


class FileDownloader(object):
    """
    Base class for file downloaders of different protocols
    """
    defaultStream = sys.stdout

    def __init__(self, url, path, stream=defaultStream):
        self.url = url
        self.path = path
        self.basename = path
        self.basenameLength = len(self.basename)
        self.stream = stream
        self.bytesReceived = 0
        self.displayIndex = 0
        self.displayWindowSize = 20
        self.fileSize = None
        self.displayCounter = 0

    def _printStartDownloadMessage(self):
        self.stream.write("Downloading '{}' to '{}'\n".format(
            self.url, self.path))

    def _cleanUp(self):
        self.stream.write("\n")
        self.stream.flush()

    def _getFileNameDisplayString(self):
        if self.basenameLength <= self.displayWindowSize:
            return self.basename
        else:
            return self.basename  # TODO scrolling window here

    def _updateDisplay(self, modulo=1):
        self.displayCounter += 1
        if self.displayCounter % modulo != 0:
            return
        fileName = self._getFileNameDisplayString()
        if self.fileSize is None:
            displayString = "{}   bytes received: {}\r"
            bytesReceived = humanize.filesize.naturalsize(
                self.bytesReceived)
            self.stream.write(displayString.format(
                fileName, bytesReceived))
        else:
            # TODO contentlength seems to slightly under-report how many
            # bytes we have to download... hence the min functions
            percentage = min(self.bytesReceived / self.fileSize, 1)
            numerator = humanize.filesize.naturalsize(
                min(self.bytesReceived, self.fileSize))
            denominator = humanize.filesize.naturalsize(
                self.fileSize)
            displayString = "{}   {:<6.2%} ({:>9} / {:<9})\r"
            self.stream.write(displayString.format(
                fileName, percentage, numerator, denominator))
        self.stream.flush()


class HttpFileDownloader(FileDownloader):
    """
    Provides a wget-like file download and terminal display for HTTP
    """
    defaultChunkSize = 1048576  # 1MB

    def __init__(self, url, path, chunkSize=defaultChunkSize,
                 stream=FileDownloader.defaultStream):
        super(HttpFileDownloader, self).__init__(
            url, path, stream)
        self.chunkSize = chunkSize

    def download(self):
        self._printStartDownloadMessage()
        response = requests.get(self.url, stream=True)
        response.raise_for_status()
        try:
            contentLength = int(response.headers['content-length'])
            self.fileSize = contentLength
        except KeyError:
            # chunked transfer encoding
            pass
        with open(self.path, 'wb') as outputFile:
            for chunk in response.iter_content(chunk_size=self.chunkSize):
                self.bytesReceived += self.chunkSize
                self._updateDisplay()
                outputFile.write(chunk)
        self._cleanUp()


def runCommandSplits(splits, silent=False):
    """
    Run a shell command given the command's parsed command line
    """
    try:
        if silent:
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call(splits, stdout=devnull, stderr=devnull)
        else:
            subprocess.check_call(splits)
    except OSError, e:
        if e.errno == 2:  # cmd not found
            raise Exception(
                "Can't find command while trying to run {}".format(splits))
        else:
            raise


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
