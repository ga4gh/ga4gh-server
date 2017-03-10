"""
Module responsible for translating continuous sequence annotation data
into GA4GH native objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import random
import re
import math

# no step/span; requires numpy
import pyBigWig

# for running bigwig tool externally
import subprocess

import ga4gh.server.datamodel as datamodel
import ga4gh.server.exceptions as exceptions
import ga4gh.schemas.pb as pb
import ga4gh.schemas.protocol as protocol


"""
These classes handle continuous message.
This message defines a format for exchanging continuous valued signal data,
such as those produced experimentally (e.g. ChIP-Seq data) or through
calculations (e.g. conservation scores). It can be used, for example,
to share data stored in Wiggle, BigWig, and BedGraph formats.

Assumes 0-based for everything (wiggle is 1-based).

The pyBiWig package is used to read the bigwig file it doesn't return step
and span, though. Other packages were
considered. The bx-python package is missing the bigwig reader in the
pypi verions. The ngslib package is not being maintained is difficult
to install.

"""


class WiggleReader:
    """
    Class for reading Wiggle data (from a file or a pipe) and returning
    protocol objects within the defined query region.

    Currently, only a single protocol object is returned, using getData().
    Unsampled values are assigned NaN in the object values array.
    The object start is the position of the first NaN value in the
    the sampled region. The values array extends to the last value
    that is not NaN. So, rather than covering the query region, the
    values array drops any initial NaN or final NaN. This was done
    to shorten the message and to avoid sending values containing only NaN.
    """

    # Wiggle has two modes:
    # 1. variable step, where each data line has a position and value.
    # 2. fixed step, where each data line has only a value. Positions
    #    are calculated as a fixed length from a start position, which
    #    are given in the header line (i.e. 'fixedStep ...').
    _VARIABLE_STEP = 1
    _FIXED_STEP = 2

    def __init__(self, reference, start, end):
        self._queryReference = reference
        self._queryStart = start
        self._queryEnd = end
        self._mode = None
        self._span = 1
        self._step = 1
        self._data = protocol.Continuous()
        self._position = None

    def getData(self):
        return self._data

    def parseStep(self, line):
        """
        Parse the line describing the mode.

        One of:
        variableStep chrom=<reference> [span=<window_size>]
        fixedStep chrom=<reference> start=<position> step=<step_interval>
                  [span=<window_size>]

        Span is optional, defaulting to 1. It indicates that each value
        applies to region, starting at the given position and extending
        <span> positions.
        """
        fields = dict([field.split('=') for field in line.split()[1:]])

        if 'chrom' in fields:
            self._reference = fields['chrom']
        else:
            raise ValueError("Missing chrom field in %s" % line.strip())

        if line.startswith("fixedStep"):
            if 'start' in fields:
                self._start = int(fields['start']) - 1  # to 0-based
            else:
                raise ValueError("Missing start field in %s" % line.strip())

        if 'span' in fields:
            self._span = int(fields['span'])
        if 'step' in fields:
            self._step = int(fields['step'])

    def readWiggleLine(self, line):
        """
        Read a wiggle line. If it is a data line, add values to the
        protocol object.
        """
        if(line.isspace() or line.startswith("#")
                or line.startswith("browser") or line.startswith("track")):
            return
        elif line.startswith("variableStep"):
            self._mode = self._VARIABLE_STEP
            self.parseStep(line)
            return
        elif line.startswith("fixedStep"):
            self._mode = self._FIXED_STEP
            self.parseStep(line)
            return
        elif self._mode is None:
            raise ValueError("Unexpected input line: %s" % line.strip())

        if self._queryReference != self._reference:
            return

        # read data lines
        fields = line.split()
        if self._mode == self._VARIABLE_STEP:
            start = int(fields[0])-1  # to 0-based
            val = float(fields[1])
        else:
            start = self._start
            self._start += self._step
            val = float(fields[0])

        if start < self._queryEnd and start > self._queryStart:
            if self._position is None:
                self._position = start
                self._data.start = start

            # fill gap
            while self._position < start:
                self._data.values.append(float('NaN'))
                self._position += 1
            for _ in xrange(self._span):
                self._data.values.append(val)
            self._position += self._span

    def wiggleFileHandleToProtocol(self, fileHandle):
        """
        Return a continuous protocol object satsifiying the given query
        parameters from the given wiggle file handle.
        """
        for line in fileHandle:
            self.readWiggleLine(line)
        return self._data

    def wiggleFileToProtocol(self, fileName):
        """
        Return a continuous protocol object satsifiying the given query
        parameters from the given wiggle file.
        """
        with open(fileName, 'r') as f:
            return self.wiggleFileHandleToProtocol(f)


class BigWigDataSource:
    """
    Class for reading from bigwig files.

    Two different readers are implemented:
    1. pyBigWig: a python package that wraps custom C code. It
       doesn't seem to return span and step values, limiting it's use.
    2. bigWigToWig: This is a command line tool for the Kent library.
       It must be installed separately.
    """

    def __init__(self, sourceFile):
        self._sourceFile = sourceFile
        self._INCREMENT = 10000  # max results per bw query
        self._MAX_VALUES = 1000  # max values length

    def checkReference(self, reference):
        """
        Check the reference for security. Tries to avoid any characters
        necessary for doing a script injection.
        """
        pattern = re.compile(r'[\s,;"\'&\\]')
        if pattern.findall(reference.strip()):
            return False
        return True

    def readValuesPyBigWig(self, reference, start, end):
        """
        Use pyBigWig package to read a BigWig file for the
        given range and return a protocol object.

        pyBigWig returns an array of values that fill the query range.
        Not sure if it is possible to get the step and span.

        This method trims NaN values from the start and end.

        pyBigWig throws an exception if end is outside of the
        reference range. This function checks the query range
        and throws its own exceptions to avoid the ones thrown
        by pyBigWig.
        """
        if not self.checkReference(reference):
            raise exceptions.ReferenceNameNotFoundException(reference)
        if start < 0:
            start = 0
        bw = pyBigWig.open(self._sourceFile)
        referenceLen = bw.chroms(reference)
        if referenceLen is None:
            raise exceptions.ReferenceNameNotFoundException(reference)
        if end > referenceLen:
            end = referenceLen
        if start >= end:
            raise exceptions.ReferenceRangeErrorException(
                reference, start, end)

        data = protocol.Continuous()
        curStart = start
        curEnd = curStart + self._INCREMENT
        while curStart < end:
            if curEnd > end:
                curEnd = end
            for i, val in enumerate(bw.values(reference, curStart, curEnd)):
                if not math.isnan(val):
                    if len(data.values) == 0:
                        data.start = curStart + i
                    data.values.append(val)
                    if len(data.values) == self._MAX_VALUES:
                        yield data
                        data = protocol.Continuous()
                elif len(data.values) > 0:
                    # data.values.append(float('NaN'))
                    yield data
                    data = protocol.Continuous()
            curStart = curEnd
            curEnd = curStart + self._INCREMENT

        bw.close()
        if len(data.values) > 0:
            yield data

    def readValuesBigWigToWig(self, reference, start, end):
        """
        Read a bigwig file and return a protocol object with values
        within the query range.

        This method uses the bigWigToWig command line tool from UCSC
        GoldenPath. The tool is used to return values within a query region.
        The output is in wiggle format, which is processed by the WiggleReader
        class.

        There could be memory issues if the returned results are large.

        The input reference can be a security problem (script injection).
        Ideally, it should be checked against a list of known chromosomes.
        Start and end should not be problems since they are integers.
        """
        if not self.checkReference(reference):
            raise exceptions.ReferenceNameNotFoundException(reference)
        if start < 0:
            raise exceptions.ReferenceRangeErrorException(
                reference, start, end)
            # TODO: CHECK IF QUERY IS BEYOND END

        cmd = ["bigWigToWig", self._sourceFile, "stdout", "-chrom="+reference,
               "-start="+str(start), "-end="+str(end)]
        wiggleReader = WiggleReader(reference, start, end)
        try:
            # run command and grab output simultaneously
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            while True:
                line = process.stdout.readline()
                if line == '' and process.poll() is not None:
                    break
                wiggleReader.readWiggleLine(line.strip())
        except ValueError:
            raise
        except:
            raise Exception("bigWigToWig failed to run")

        return wiggleReader.getData()

    def bigWigToProtocol(self, reference, start, end):
        # return self.readValuesBigWigToWig(reference, start, end)
        for continuousObj in self.readValuesPyBigWig(reference, start, end):
            yield continuousObj


class AbstractContinuousSet(datamodel.DatamodelObject):
    """
    A continuous sequence annotation set
    """
    compoundIdClass = datamodel.ContinuousSetCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractContinuousSet, self).__init__(parentContainer, localId)
        self._name = localId
        self._sourceUri = ""
        self._referenceSet = None

    def getReferenceSet(self):
        """
        Returns the reference set associated with this ContinuousSet.
        """
        return self._referenceSet

    def setReferenceSet(self, referenceSet):
        """
        Sets the reference set associated with this ContinuousSet to the
        specified value.
        """
        self._referenceSet = referenceSet

    def toProtocolElement(self):
        """
        Returns the representation of this ContinuousSet as the corresponding
        ProtocolElement.
        """
        gaContinuousSet = protocol.ContinuousSet()
        gaContinuousSet.id = self.getId()
        gaContinuousSet.dataset_id = self.getParentContainer().getId()
        gaContinuousSet.reference_set_id = pb.string(
                                            self._referenceSet.getId())
        gaContinuousSet.name = self._name
        gaContinuousSet.source_uri = self._sourceUri
        attributes = self.getAttributes()
        for key in attributes:
            gaContinuousSet.attributes.attr[key] \
                .values.extend(protocol.encodeValue(attributes[key]))
        return gaContinuousSet


class FileContinuousSet(AbstractContinuousSet):
    """
    Data associated with a file containing continuous data.
    """
    def __init__(self, parentContainer, localId):
        super(FileContinuousSet, self).__init__(parentContainer, localId)
        self._filePath = None

    def populateFromFile(self, dataUrl):
        """
        Populates the instance variables of this ContinuousSet from the
        specified data URL.
        """
        self._filePath = dataUrl

    def populateFromRow(self, continuousSetRecord):
        """
        Populates the instance variables of this ContinuousSet from the
        specified DB row.
        """
        self._filePath = continuousSetRecord.dataurl
        self.setAttributesJson(continuousSetRecord.attributes)

    def getDataUrl(self):
        """
        Returns the URL providing the data source for this ContinuousSet.
        """
        return self._filePath

    def getContinuous(self, referenceName=None, start=None, end=None):
        """
        Method passed to runSearchRequest to fulfill the request to
        yield continuous protocol objects that satisfy the given query.

        :param str referenceName: name of reference (ex: "chr1")
        :param start: castable to int, start position on reference
        :param end: castable to int, end position on reference
        :return: yields a protocol.Continuous at a time
        """
        bigWigReader = BigWigDataSource(self._filePath)
        for continuousObj in bigWigReader.bigWigToProtocol(
                                            referenceName, start, end):
            yield continuousObj


class SimulatedContinuousSet(AbstractContinuousSet):
    """
    Simulated data backend for ContinuousSet, used for internal testing.
    """
    def __init__(self, parentContainer, localId, randomSeed=1):
        self._randomSeed = randomSeed
        super(SimulatedContinuousSet, self).__init__(parentContainer, localId)

    def _generateSimulatedContinuous(self, randomNumberGenerator):
        continuous = protocol.Continuous()
        continuous.start = randomNumberGenerator.randint(1000, 2000)
        continuous.values = [100, 200.3, 400]

    def getContinuousData(self, referenceName=None, start=None, end=None):
        """
        Returns a set number of simulated continuous data.

        :param referenceName: name of reference to "search" on
        :param start: start coordinate of query
        :param end: end coordinate of query
        :return: Yields continuous list
        """
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(self._randomSeed)
        for i in range(100):
            gaContinuous = self._generateSimulatedContinuous(
                                    randomNumberGenerator)
            match = (
                gaContinuous.start < end and
                gaContinuous.end > start and
                gaContinuous.reference_name == referenceName)
            if match:
                yield gaContinuous
