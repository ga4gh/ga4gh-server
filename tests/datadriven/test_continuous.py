"""
Unit tests for continuous objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import math

from nose.tools import raises

import ga4gh.server.datarepo as datarepo
import ga4gh.server.datamodel.continuous as continuous
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.exceptions as exceptions

import tests.paths as paths


class TestContinuous(unittest.TestCase):
    """
    Unit tests for continuous data
    """
    def _createContinuousSet(self):
        """
        Creates a ContinuousSet from the specified directory.
        """
        self._continuousSetName = "testContinuous"
        self._repo = datarepo.SqlDataRepository(paths.testDataRepo)
        self._repo.open(datarepo.MODE_READ)
        self._dataset = datasets.Dataset("testDs")
        self._continuousSet = continuous.readSet(
            self._dataset, self._continuousSetName)

    def setUp(self):
        dataDir = "tests/data/datasets/dataset1/continuous"
        self._wiggleFile = dataDir + "/wiggle_2.txt"
        self._bigWigFile = dataDir + "/bigwig_1.bw"

    def testReadWiggle(self):
        continuousObj = continuous.WiggleReader(
                            'chr19', 49307698, 49308020)
        obj = continuousObj.wiggleFileToProtocol(self._wiggleFile)
        self.assertEqual(obj.start, 49307700)
        self.assertEqual(obj.values[0], 900)
        self.assertEqual(obj.values[300], 800)
        self.assertEqual(len(obj.values), 302)

    def getTuples(self, generator):
        """
        Convert a generator of continuous objects into tuples of
        (position,value).
        """
        tuples = []
        for obj in generator:
            for i, value in enumerate(obj.values):
                if not math.isnan(value):
                    tuples.append((obj.start+i, value))
        return tuples

    def testReadBigWig(self):
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        generator = continuousObj.bigWigToProtocol("chr19", 49305897, 49306090)
        tuples = self.getTuples(generator)
        self.assertEqual(tuples[0], (49305900, 20.0))
        self.assertEqual(tuples[4], (49305904, 20.0))
        self.assertEqual(tuples[5], (49306080, 17.5))
        self.assertEqual(tuples[9], (49306084, 17.5))
        self.assertEqual(len(tuples), 10)

    def testReadBigWigAllNan(self):
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        generator = continuousObj.bigWigToProtocol(
                                        "chr19", 49305927, 49305997)
        tuples = self.getTuples(generator)
        self.assertEqual(len(tuples), 0)

    @raises(exceptions.ReferenceRangeErrorException)
    def testReadBigWigInvalidRange(self):
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        generator = continuousObj.bigWigToProtocol(
                                        "chr19", 493059030, 49305934)
        next(generator)

    def testReadBigWigOutsideReferenceRange(self):
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        generator = continuousObj.bigWigToProtocol(
                                    "chr19", 49306897, 493059304)
        tuples = self.getTuples(generator)
        self.assertEqual(len(tuples), 5)

    def testReadBigWigNegativeReferenceRange(self):
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        generator = continuousObj.bigWigToProtocol("chr19", -1, 5)
        tuples = self.getTuples(generator)
        self.assertEqual(len(tuples), 0)

    @raises(exceptions.ReferenceNameNotFoundException)
    def testReadBigWigChromsomeException(self):
        """
        Test for catching bad chromosome names.
        """
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        generator = continuousObj.bigWigToProtocol(
                                            "chr&19", 49305602, 49308000)
        next(generator)
