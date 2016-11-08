"""
Random data based tests for the interval iterator code.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import random

import ga4gh.server.paging as paging


def setUp():
    """
    Set the random seed to ensure that the tests in this module are
    reproducable.
    """
    random.seed(1)


def intervalsIntersect(a1, b1, a2, b2):
    """
    Returns True if the specified half-closed intervals [a1, b1)
    and [a2, b2) intersect.
    """
    # Intervals do not interset if [a1, b1) is wholly to the left
    # or right of [a2, b2).
    return not (b1 <= a2 or a1 >= b2)


class TestIntervalsIntersect(unittest.TestCase):
    """
    Tests to ensure our interval intersection code works.
    """
    def testEmptyIntervals(self):
        self.assertFalse(intervalsIntersect(0, 0, 0, 0))
        self.assertFalse(intervalsIntersect(1, 1, 1, 1))
        self.assertFalse(intervalsIntersect(1, 1, 1, 2))
        self.assertFalse(intervalsIntersect(1, 10, 1, 1))

    def testNonIntersecting(self):
        tests = [
            [0, 1, 1, 2],
            [0, 1, 2, 3],
            [2, 1, 0, 1],
            [2, 3, 0, 1]]
        for test in tests:
            self.assertFalse(intervalsIntersect(*test))

    def testIntersecting(self):
        tests = [
            [0, 1, 0, 1],
            [0, 1, 0, 3],
            [0, 3, 0, 1],
            [1, 10, 9, 10],
            [1, 10, 9, 12],
            [1, 10, 5, 6],
            [1, 10, 1, 2]]
        for test in tests:
            self.assertTrue(intervalsIntersect(*test))


def randomIntervals(start, end, numIntervals):
    """
    Returns a list of numIntervals half-closed intervals. Each element
    of the returned array is a tuple a, b such that a <= start < b <= end.
    """
    intervals = []
    for _ in range(numIntervals):
        a = random.randrange(start, end - 1)
        b = random.randrange(a + 1, end)
        assert start <= a < b <= end
        intervals.append((a, b))
    return intervals


class IntervalSet(object):
    """
    Class representing a set of half-closed intervals.
    """
    def __init__(self, start, end, intervals):
        self.start = start
        self.end = end
        self.intervals = sorted(intervals, key=lambda x: x[0])

    def get(self, start, end):
        """
        Returns an iterator over all intervals in this set that intersect
        with the specified interval.
        """
        for interval in self.intervals:
            if intervalsIntersect(start, end, interval[0], interval[1]):
                yield interval


class FakeRequest(object):
    """
    A class to stand in as a Request object.
    """
    def __init__(self, start, end, pageToken):
        self.page_token = pageToken
        self.start = start
        self.end = end


class TrivialIntervalIterator(paging.IntervalIterator):
    """
    The simplest possible instance of the interval iterator
    used to test the iteration code.
    """
    def __init__(self, intervalSet, start, end, pageToken=None):
        self.intervalSet = intervalSet
        request = FakeRequest(start, end, pageToken)
        super(TrivialIntervalIterator, self).__init__(request, None)

    def _getContainer(self):
        return None

    def _search(self, start, end):
        return self.intervalSet.get(start, end)

    def _getStart(self, interval):
        return interval[0]

    def _getEnd(self, interval):
        return interval[1]


class TestIntervalIterator(unittest.TestCase):
    """
    A class to systematically test the paging code over interval search
    by randomly generating densely packed interval data, and comparing the
    results for a range of page sizes.
    """
    def setUp(self):
        self.num_random_tests = 5
        self.testIntervalSets = []
        # Set up some interval sets that provoke specific issues
        intervals = [
            (0, 1), (1, 8), (2, 9), (4, 7), (4, 8), (5, 9), (6, 7), (6, 7),
            (7, 8), (8, 9)]
        self.testIntervalSets.append(IntervalSet(0, 10, intervals))
        # Make two interval sets, and merge them. We want to have patches
        # of intervals, but with queries happening in between them.
        intervals = randomIntervals(10, 20, 10) + randomIntervals(30, 40, 10)
        self.testIntervalSets.append(IntervalSet(0, 50, intervals))
        # Add other random interval sets
        self.testIntervalSets.append(
            IntervalSet(0, 10, randomIntervals(0, 10, 10)))
        self.testIntervalSets.append(
            IntervalSet(0, 10, randomIntervals(0, 10, 10)))
        self.testIntervalSets.append(
            IntervalSet(0, 100, randomIntervals(0, 100, 10)))
        self.testIntervalSets.append(
            IntervalSet(0, 100, randomIntervals(0, 100, 100)))

    def verifyInterval(self, intervalSet, start, end):
        """
        Verify that we can pick up iteration of the interval from
        anywhere by starting a new iterator from every point.
        """
        topIterator = list(TrivialIntervalIterator(intervalSet, start, end))
        allIntervals = list(intervalSet.get(start, end))
        topIntervals = []
        for topInterval, topPageToken in topIterator[:-1]:
            topIntervals.append(topInterval)
            self.assertIsNotNone(topPageToken)
            # We should be able to pick the iteration up from here and go
            # forward, getting the same set of intervals
            subIterator = TrivialIntervalIterator(
                intervalSet, start, end, topPageToken)
            subIntervals = list(topIntervals)
            for subInterval, subPageToken in subIterator:
                subIntervals.append(subInterval)
            self.assertEqual(allIntervals, subIntervals)
            self.assertIsNone(subPageToken)
        topInterval, topPageToken = topIterator[-1]
        self.assertIsNone(topPageToken)
        topIntervals.append(topInterval)
        self.assertEqual(allIntervals, topIntervals)

    def testFullInterval(self):
        for intervalSet in self.testIntervalSets:
            self.verifyInterval(
                intervalSet, intervalSet.start, intervalSet.end)

    def verifyEmptyInterval(self, intervalSet, start, end):
        """
        Verify that we correctly return an empty iterator.
        """
        iterator = TrivialIntervalIterator(intervalSet, start, end)
        self.assertIsNone(next(iterator, None))

    def testEmptyInterval(self):
        for intervalSet in self.testIntervalSets:
            for start, end in [(-1, -1), (intervalSet.end, intervalSet.end)]:
                self.verifyEmptyInterval(intervalSet, start, end)

    def testStartInGap(self):
        for intervalSet in self.testIntervalSets:
            starts = set(start for start, _ in intervalSet.intervals)
            gaps = set(range(intervalSet.start, intervalSet.end - 1)) - starts
            start = -1
            if len(gaps) > 0:
                start = list(gaps)[0]
            self.verifyInterval(intervalSet, start, intervalSet.end)

    def testOutsideRange(self):
        for intervalSet in self.testIntervalSets:
            self.verifyInterval(
                intervalSet, intervalSet.start - 1, intervalSet.end)
            self.verifyInterval(
                intervalSet, intervalSet.start, intervalSet.end + 1)
            self.verifyInterval(
                intervalSet, intervalSet.start - 1, intervalSet.end + 1)
            self.verifyInterval(
                intervalSet, intervalSet.start - 10, intervalSet.end + 10)

    def testRandomStarts(self):
        for intervalSet in self.testIntervalSets:
            for _ in range(self.num_random_tests):
                start = random.randrange(
                    intervalSet.start, intervalSet.end - 1)
                end = intervalSet.end
                intervals = list(intervalSet.get(start, end))
                if len(intervals) == 0:
                    self.verifyEmptyInterval(intervalSet, start, end)
                else:
                    self.verifyInterval(intervalSet, start, end)

    def testRandomEnds(self):
        for intervalSet in self.testIntervalSets:
            for _ in range(self.num_random_tests):
                start = intervalSet.start
                end = random.randrange(
                    intervalSet.start + 1, intervalSet.end)
                intervals = list(intervalSet.get(start, end))
                if len(intervals) == 0:
                    self.verifyEmptyInterval(intervalSet, start, end)
                else:
                    self.verifyInterval(intervalSet, start, end)

    def testRandomIntervals(self):
        for intervalSet in self.testIntervalSets:
            for _ in range(self.num_random_tests):
                start = random.randrange(
                    intervalSet.start, intervalSet.end - 1)
                end = random.randrange(start, intervalSet.end - 1)
                intervals = list(intervalSet.get(start, end))
                if len(intervals) == 0:
                    self.verifyEmptyInterval(intervalSet, start, end)
                else:
                    self.verifyInterval(intervalSet, start, end)
