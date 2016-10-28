"""
Simple tool to download files
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys

import humanize

import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()


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
