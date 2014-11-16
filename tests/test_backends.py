"""
Tests for the backend objects. We instantiate local copies of the backends
and invoke the entry points for the protocol methods. We do not set up
any server processes or communicate over sockets.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import shutil
import tempfile
import unittest
import itertools
import subprocess

import wormtable as wt

import tests
import ga4gh.server as server
import ga4gh.protocol as protocol


class WormtableTestFixture(object):
    """
    Class representing the test fixtures used for the wormtable tests.
    This is a temporary directory containing converted wormtables
    for example VCFs in the test directory.
    """
    def __init__(self):
        self.dataDir = tempfile.mkdtemp(prefix="ga4gh_wt")

    def convertVariantSet(self, vcfFile):
        """
        Converts the specified VCF file into wormtable format, storing it
        the data directory.
        """
        vsid = vcfFile.split("/")[-1].split(".")[0]
        wtDir = os.path.join(self.dataDir, vsid)
        # convert the vcf to wormtable format.
        cmd = ["vcf2wt", "-q", vcfFile, wtDir]
        subprocess.check_call(cmd)
        # Build the indexes
        # The CHROM index is used here to find the distinct reference names.
        for index in ["CHROM", "CHROM+POS", "CHROM+ID"]:
            cmd = ["wtadmin", "add", "-q", wtDir, index]
            subprocess.check_call(cmd)

    def setUp(self):
        for f in glob.glob("tests/data/*.vcf.gz"):
            self.convertVariantSet(f)

    def tearDown(self):
        shutil.rmtree(self.dataDir)

# Test fixture module variable
_wormtableTestFixture = None


def setUp():
    """
    Module level setup for nosetests.
    """
    global _wormtableTestFixture
    _wormtableTestFixture = WormtableTestFixture()
    _wormtableTestFixture.setUp()


def tearDown():
    """
    Module level teardown for nosetests.
    """
    global _wormtableTestFixture
    _wormtableTestFixture.tearDown()


def convertInfoValue(column, value):
    """
    Converts a value from the specified info column into the format
    expected by the protocol.
    """
    if column.get_num_elements() == 1 or column.get_type() == wt.WT_CHAR:
        ret = [str(value)]
    else:
        ret = [str(v) for v in value]
    return ret


class TestWormtableBackend(unittest.TestCase):

    def setUp(self):
        global _wormtableTestFixture
        self._dataDir = _wormtableTestFixture.dataDir
        self._tables = {}
        self._chromIndexes = {}
        self._chromPosIndexes = {}
        for f in os.listdir(self._dataDir):
            t = wt.open_table(os.path.join(self._dataDir, f))
            self._tables[f] = t
            self._chromIndexes[f] = t.open_index("CHROM")
            self._chromPosIndexes[f] = t.open_index("CHROM+POS")
        self._backend = server.WormtableBackend(self._dataDir)

    def tearDown(self):
        for t in self._tables.values():
            t.close()

    def resultIterator(
            self, request, pageSize, searchMethod, ResponseClass, listMember):
        """
        Returns an iterator over the list of results from the specified
        request.  All results are returned, and paging is handled
        automatically.
        """
        not_done = True
        request.pageSize = pageSize
        while not_done:
            response = searchMethod(request)
            self.assertEqual(type(response), ResponseClass)
            l = getattr(response, listMember)
            self.assertLessEqual(len(l), pageSize)
            for vs in l:
                yield vs
            not_done = response.nextPageToken is not None
            request.pageToken = response.nextPageToken

    def getVariantSets(self, pageSize=100):
        """
        Returns an iterator over the variantSets, abstracting away the details
        of the pageSize.
        """
        request = protocol.GASearchVariantsRequest()
        return self.resultIterator(
            request, pageSize, self._backend.searchVariantSets,
            protocol.GASearchVariantSetsResponse, "variantSets")

    def getVariants(
            self, variantSetIds, referenceName, start=0, end=2**32,
            pageSize=100, callSetIds=[]):
        """
        Returns an iterator over the specified list of variants, abstracting
        out paging details.
        """
        request = protocol.GASearchVariantsRequest()
        request.variantSetIds = variantSetIds
        request.referenceName = referenceName
        request.start = start
        request.end = end
        request.callSetIds = callSetIds
        return self.resultIterator(
            request, pageSize, self._backend.searchVariants,
            protocol.GASearchVariantsResponse, "variants")

    def getReferenceNames(self, variantSetId):
        """
        Returns the distinct reference names for the specified variantSetId.
        """
        i = self._chromIndexes[variantSetId]
        for k in i.counter().keys():
            yield k

    def getCommonRefNames(self, variantSetIds):
        """
        Returns the reference names common to all the specified variantSetIds.
        """
        assert len(variantSetIds) > 0
        commonNames = set(self.getReferenceNames(variantSetIds[0]))
        for vsid in variantSetIds[1:]:
            commonNames &= set(self.getReferenceNames(vsid))
        return commonNames

    def getWormtableVariants(
            self, variantSetIds, referenceName, start=0, end=2**32,
            callSetIds=[]):
        """
        Returns the rows from the table corresponding to the specified values.
        Each row is a dictionary mapping the column to its value. If callSetIds
        is empty, return all sample specific columns. Otherwise, only return
        columns sample specific columns for the callSetIds in question.
        """
        for vsid in variantSetIds:
            t = self._tables[vsid]
            i = self._chromPosIndexes[vsid]
            s = (referenceName, start)
            e = (referenceName, end)
            cols = t.columns()[1:8]
            for c in t.columns()[8:]:
                name = c.get_name()
                if name.startswith("INFO"):
                    cols.append(c)
                else:
                    if len(callSetIds) == 0:
                        cols.append(c)
                    else:
                        csid = name.split(".")[0]
                        if csid in callSetIds:
                            cols.append(c)
            for r in i.cursor(cols, s, e):
                yield dict((c, v) for c, v in zip(cols, r))

    def getStartBounds(self, variantSetId, referenceName):
        """
        Returns the first and last start coorinate for the specified
        referenceName in the specified variantSet.
        """
        i = self._chromPosIndexes[variantSetId]
        first = i.min_key(referenceName)
        last = i.max_key(referenceName)
        return first[1], last[1]

    def getCallSetIds(self, variantSetId):
        """
        Returns the set of callSetIds derived from the wormtable header.
        """
        t = self._tables[variantSetId]
        s = set()
        for c in t.columns()[8:]:  # skip VCF fixed cols
            name = c.get_name()
            if not name.startswith("INFO"):
                s.add(name.split(".")[0])
        return s


class TestVariantSets(TestWormtableBackend):
    """
    Test the searchVariantSets end point
    """

    def testPagination(self):
        vsids = []
        for ps in range(1, 10):
            s = [vs.id for vs in self.getVariantSets(pageSize=ps)]
            vsids.append(s)
        for v in vsids[1:]:
            self.assertEqual(v, vsids[0])
        # TODO repeat this for full object equality, not just ids

    def testIds(self):
        variantSets = [vs for vs in self.getVariantSets()]
        self.assertEqual(len(variantSets), len(self._tables))
        ids = set(vs.id for vs in variantSets)
        self.assertEqual(ids, set(self._tables.keys()))


class TestVariants(TestWormtableBackend):
    """
    Tests the searchVariants end point.
    """

    def verifyGenotype(self, call, vcfString):
        """
        Verifies the specified Call object has the correct interpretation of
        the specified VCF genotype value.
        """
        # TODO parse the VCF string and verify that it's correct.

    def verifyVariantsEqual(self, variant, row):
        """
        Verifies that the specified protocol Variant object is equal to the
        specified wormtable row. The row is a dictionary mapping columns
        to values.
        """
        r = {}
        columns = {}
        for c, value in row.items():
            name = c.get_name()
            columns[name] = c
            r[name] = value
            # Check the info cols.
            if name.startswith("INFO."):
                infoField = name.split(".")[1]
                if value is None:
                    # Either it is not in the map, or it maps to None
                    if infoField in variant.info:
                        self.assertNone(variant.info[infoField])
                else:
                    v = variant.info[infoField]
                    self.assertEqual(v, convertInfoValue(c, value))
        self.assertEqual(variant.referenceName, r["CHROM"])
        self.assertEqual(variant.start, r["POS"])
        ref = r["REF"]
        self.assertEqual(variant.referenceBases, ref)
        alt = r["ALT"]
        if alt is None:
            alt = []
        else:
            alt = alt.split(",")
        self.assertEqual(variant.alternateBases, alt)
        end = variant.start + len(ref)
        self.assertEqual(variant.end, end)
        names = r["ID"]
        if names is None:
            names = []
        else:
            names = names.split(",")
        self.assertEqual(variant.names, names)
        for call in variant.calls:
            csid = call.callSetId
            k = csid + ".GT"
            self.assertTrue(k in r)
            self.verifyGenotype(call, r[k])
            k = csid + ".GL"
            if k in r:
                self.assertTrue(False)  # TODO add test case and data.
            else:
                self.assertEqual(len(call.genotypeLikelihood), 0)
            for infoField, v in call.info.items():
                k = csid + "." + infoField
                self.assertTrue(k in r)
                self.assertEqual(v, convertInfoValue(columns[k], r[k]))

    def verifySearchVariants(
            self, variantSetIds, referenceName, start, end, pageSize=1000):
        """
        Verifies queries on the specified reference name and list of variant
        set ids.
        """
        l1 = [v for v in self.getVariants(
            variantSetIds, referenceName, start=start, end=end)]
        l2 = [r for r in self.getWormtableVariants(
            variantSetIds, referenceName, start=start, end=end)]
        self.assertEqual(len(l1), len(l2))
        for variant, row in zip(l1, l2):
            self.verifyVariantsEqual(variant, row)

    def verifySearchByCallSetIds(
            self, vsid, referenceName, start, end, callSetIds, pageSize=1000):
        """
        Verifies the variants we get back contain the correct callSet
        information.
        """
        l1 = [v for v in self.getVariants(
            [vsid], referenceName, start=start, end=end,
            callSetIds=callSetIds)]
        l2 = [r for r in self.getWormtableVariants(
            [vsid], referenceName, start=start, end=end,
            callSetIds=callSetIds)]
        self.assertEqual(len(l1), len(l2))
        for variant, row in zip(l1, l2):
            self.verifyVariantsEqual(variant, row)
            # Verify we've got the correct csids.
            s1 = set(cs.callSetId for cs in variant.calls)
            self.assertEqual(len(s1), len(variant.calls))
            self.assertEqual(s1, set(callSetIds))

    def testSearchAllVariants(self):
        for vs in self.getVariantSets():
            for referenceName in self.getReferenceNames(vs.id):
                self.verifySearchVariants([vs.id], referenceName, 0, 2**32)

    def testSearchVariantSlices(self):
        for vs in self.getVariantSets():
            for referenceName in self.getReferenceNames(vs.id):
                first, last = self.getStartBounds(vs.id, referenceName)
                mid = (last - first) // 2
                for pageSize in [1, 2, 3, 5, 1000]:
                    self.verifySearchVariants(
                        [vs.id], referenceName, first, last, pageSize=pageSize)
                    self.verifySearchVariants(
                        [vs.id], referenceName, first, mid, pageSize=pageSize)
                    self.verifySearchVariants(
                        [vs.id], referenceName, mid, last, pageSize=pageSize)
                    self.verifySearchVariants(
                        [vs.id], referenceName, first, first + 1,
                        pageSize=pageSize)
                    self.verifySearchVariants(
                        [vs.id], referenceName, mid, mid + 1,
                        pageSize=pageSize)
                    self.verifySearchVariants(
                        [vs.id], referenceName, mid - 1, mid + 1,
                        pageSize=pageSize)
                    self.verifySearchVariants(
                        [vs.id], referenceName, last, last + 1,
                        pageSize=pageSize)

    def testSearchByCallSetIds(self):
        for vs in self.getVariantSets():
            csids = self.getCallSetIds(vs.id)
            # limit the subsets we check over to some small value.
            for subset in tests.powerset(csids, 20):
                for referenceName in self.getReferenceNames(vs.id):
                    self.verifySearchByCallSetIds(
                        vs.id, referenceName, 0, 2**32, subset)

    def testUniqueIds(self):
        ids = set()
        for vs in self.getVariantSets():
            for referenceName in self.getReferenceNames(vs.id):
                for v in self.getVariants([vs.id], referenceName):
                    self.assertTrue(v.id not in ids)
                    ids.add(v.id)

    def testAcrossVariantSets(self):
        allVariantSets = [vs.id for vs in self.getVariantSets()]
        for permLen in range(1, len(allVariantSets) + 1):
            for variantSets in itertools.permutations(allVariantSets, permLen):
                for referenceName in self.getCommonRefNames(variantSets):
                    for pageSize in [1, 2, 3, 5, 1000]:
                        self.verifySearchVariants(
                            variantSets, referenceName, 0, 2**32,
                            pageSize=pageSize)
