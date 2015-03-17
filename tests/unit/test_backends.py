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

import tests.utils as utils
import ga4gh.backend as backend
import ga4gh.protocol as protocol


class WormtableTestFixture(object):
    """
    Class representing the test fixtures used for the wormtable tests.
    This is a temporary directory containing converted wormtables
    for example VCFs in the test directory.
    """
    def __init__(self):
        self.dataDir = tempfile.mkdtemp(prefix="ga4gh_wt")
        self.variantsDir = os.path.join(self.dataDir, "variants")
        self.referencesDir = os.path.join(self.dataDir, "references")
        self.readsDir = os.path.join(self.dataDir, "reads")
        subdirs = [self.variantsDir, self.referencesDir, self.readsDir]
        for subdir in subdirs:
            os.mkdir(subdir)

    def convertVariantSet(self, vcfFile):
        """
        Converts the specified VCF file into wormtable format, storing it
        the data directory.
        """
        variantSetid = vcfFile.split("/")[-1].split(".")[0] + '.wt'
        wtDir = os.path.join(self.variantsDir, variantSetid)
        # convert the vcf to wormtable format.
        cmd = ["vcf2wt", "-q", vcfFile, wtDir]
        subprocess.check_call(cmd)
        # Build the indexes
        # The CHROM index is used here to find the distinct reference names.
        for index in ["CHROM", "CHROM+POS", "CHROM+ID"]:
            cmd = ["wtadmin", "add", "-q", wtDir, index]
            subprocess.check_call(cmd)

    def setUp(self):
        for filename in glob.glob("tests/data/variants/example_*/*.vcf.gz"):
            self.convertVariantSet(filename)

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
        ret = [str(listElement) for listElement in value]
    return ret


class TestWormtableBackend(unittest.TestCase):

    def setUp(self):
        global _wormtableTestFixture
        self._dataDir = _wormtableTestFixture.dataDir
        self._variantsDir = _wormtableTestFixture.variantsDir
        self._tables = {}
        self._chromIndexes = {}
        self._chromPosIndexes = {}
        for relativePath in os.listdir(self._variantsDir):
            table = wt.open_table(
                os.path.join(self._variantsDir, relativePath))
            self._tables[relativePath] = table
            self._chromIndexes[relativePath] = table.open_index("CHROM")
            self._chromPosIndexes[relativePath] = table.open_index("CHROM+POS")
        self._backend = backend.FileSystemBackend(self._dataDir)

    def tearDown(self):
        for table in self._tables.values():
            table.close()

    def resultIterator(
            self, request, pageSize, searchMethod, ResponseClass, listMember):
        """
        Returns an iterator over the list of results from the specified
        request.  All results are returned, and paging is handled
        automatically.
        """
        notDone = True
        request.pageSize = pageSize
        while notDone:
            # TODO validate the response there.
            responseStr = searchMethod(request.toJsonString())
            response = ResponseClass.fromJsonString(responseStr)
            objectList = getattr(response, listMember)
            self.assertLessEqual(len(objectList), pageSize)
            for obj in objectList:
                yield obj
            notDone = response.nextPageToken is not None
            request.pageToken = response.nextPageToken

    def getVariantSets(self, pageSize=100):
        """
        Returns an iterator over the variantSets, abstracting away the details
        of the pageSize.
        """
        request = protocol.GASearchVariantSetsRequest()
        return self.resultIterator(
            request, pageSize, self._backend.searchVariantSets,
            protocol.GASearchVariantSetsResponse, "variantSets")

    def getVariantSetIds(self):
        """
        Return all variantSetIds.
        """
        return [variantSet.id for variantSet in self.getVariantSets()]

    def getVariants(
            self, variantSetIds, referenceName, start=0, end=2 ** 32,
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

    def getCallSets(self, variantSetId, pageSize=100):
        """
        Returns an iterator over the callsets in a specified variant set.
        """
        request = protocol.GASearchCallSetsRequest()
        request.variantSetIds = [variantSetId]
        return self.resultIterator(
            request, pageSize, self._backend.searchCallSets,
            protocol.GASearchCallSetsResponse, "callSets")

    def getReferenceNames(self, variantSetId):
        """
        Returns the distinct reference names for the specified variantSetId.
        """
        index = self._chromIndexes[variantSetId]
        for referenceName in index.counter().keys():
            yield referenceName

    def getCommonRefNames(self, variantSetIds):
        """
        Returns the reference names common to all the specified variantSetIds.
        """
        self.assertGreater(len(variantSetIds), 0)
        commonNames = set(self.getReferenceNames(variantSetIds[0]))
        for variantSetId in variantSetIds[1:]:
            commonNames &= set(self.getReferenceNames(variantSetId))
        return commonNames

    def getWormtableVariants(
            self, variantSetIds, referenceName, start=0, end=2 ** 32,
            callSetIds=[]):
        """
        Returns the rows from the table corresponding to the specified values.
        Each row is a dictionary mapping the column to its value. If callSetIds
        is empty, return all sample specific columns. Otherwise, only return
        columns sample specific columns for the callSetIds in question.
        """
        for variantSetId in variantSetIds:
            table = self._tables[variantSetId]
            index = self._chromPosIndexes[variantSetId]
            cols = table.columns()[1:8]
            for col in table.columns()[8:]:
                name = col.get_name()
                if name.startswith("INFO"):
                    cols.append(col)
                else:
                    if len(callSetIds) == 0:
                        cols.append(col)
                    else:
                        callSetId = name.split(".")[0]
                        if callSetId in callSetIds:
                            cols.append(col)
            cursor = index.cursor(
                cols, (referenceName, start), (referenceName, end))
            for row in cursor:
                yield dict((col, value) for col, value in zip(cols, row))

    def getStartBounds(self, variantSetId, referenceName):
        """
        Returns the first and last start coorinate for the specified
        referenceName in the specified variantSet.
        """
        index = self._chromPosIndexes[variantSetId]
        first = index.min_key(referenceName)
        last = index.max_key(referenceName)
        return first[1], last[1]

    def getWormtableCallSetIds(self, variantSetId):
        """
        Returns the set of callSetIds derived from the wormtable header.
        """
        table = self._tables[variantSetId]
        callSetIds = list()
        for col in table.columns()[8:]:  # skip VCF fixed cols
            name = col.get_name().split(".")[0]
            if not name.startswith("INFO") and name not in callSetIds:
                callSetIds.append(name.split(".")[0])
        return callSetIds


class TestVariantSets(TestWormtableBackend):
    """
    Test the searchVariantSets end point
    """

    def testPagination(self):
        results = []
        for pageSize in range(1, 10):
            variantSetIds = [
                variantSet.id for variantSet in self.getVariantSets(
                    pageSize=pageSize)]
            results.append(variantSetIds)
        for result in results[1:]:
            self.assertEqual(result, results[0])
        # TODO repeat this for full object equality, not just ids

    def testIds(self):
        variantSets = [variantSet for variantSet in self.getVariantSets()]
        self.assertEqual(len(variantSets), len(self._tables))
        ids = set(variantSet.id for variantSet in variantSets)
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
        if vcfString is None:
            genotype = [-1]
        else:
            delim = "/"
            if "|" in vcfString:
                delim = "|"
            if "." in vcfString:
                genotype = [-1]
            else:
                genotype = map(int, vcfString.split(delim))
        self.assertEqual(call.genotype, genotype)

    def verifyVariantsEqual(self, variant, columnValueMap):
        """
        Verifies that the specified protocol Variant object is equal to the
        specified wormtable row, provided in the form of a dictionary
        mapping column objects to values.
        """
        nameValueMap = {}
        nameColumnMap = {}
        for col, value in columnValueMap.items():
            name = col.get_name()
            nameColumnMap[name] = col
            nameValueMap[name] = value
            # Check the info cols.
            if name.startswith("INFO."):
                infoField = name.split(".")[1]
                if value is None:
                    # Either it is not in the map, or it maps to None
                    if infoField in variant.info:
                        self.assertNone(variant.info[infoField])
                else:
                    infoValue = variant.info[infoField]
                    self.assertEqual(infoValue, convertInfoValue(col, value))
        self.assertEqual(variant.referenceName, nameValueMap["CHROM"])
        self.assertEqual(variant.start, nameValueMap["POS"])
        ref = nameValueMap["REF"]
        self.assertEqual(variant.referenceBases, ref)
        alt = nameValueMap["ALT"]
        if alt is None:
            alt = []
        else:
            alt = alt.split(",")
        self.assertEqual(variant.alternateBases, alt)
        end = variant.start + len(ref)
        self.assertEqual(variant.end, end)
        names = nameValueMap["ID"]
        if names is None:
            names = []
        else:
            names = names.split(",")
        self.assertEqual(variant.names, names)
        for call in variant.calls:
            callSetId = call.callSetId
            colName = callSetId + ".GT"
            self.assertTrue(colName in nameValueMap)
            self.verifyGenotype(call, nameValueMap[colName])
            colName = callSetId + ".GL"
            if colName in nameValueMap:
                self.assertTrue(False)  # TODO add test case and data.
            else:
                self.assertEqual(len(call.genotypeLikelihood), 0)
            for infoField, value in call.info.items():
                colName = callSetId + "." + infoField
                self.assertTrue(colName in nameValueMap)
                self.assertEqual(
                    value, convertInfoValue(
                        nameColumnMap[colName], nameValueMap[colName]))

    def verifySearchVariants(
            self, variantSetIds, referenceName, start, end, pageSize=1000):
        """
        Verifies queries on the specified reference name and list of variant
        set ids.
        """
        ga4ghVariants = list(self.getVariants(
            variantSetIds, referenceName, start=start, end=end))
        wormtableVariants = list(self.getWormtableVariants(
            variantSetIds, referenceName, start=start, end=end))
        self.assertEqual(len(ga4ghVariants), len(wormtableVariants))
        for ga4ghVariant, wormtableVariant in zip(
                ga4ghVariants, wormtableVariants):
            self.verifyVariantsEqual(ga4ghVariant, wormtableVariant)

    def verifySearchByCallSetIds(
            self, variantSetId, referenceName, start, end, callSetIds,
            pageSize=1000):
        """
        Verifies the variants we get back contain the correct callSet
        information.
        """
        ga4ghVariants = list(self.getVariants(
            [variantSetId], referenceName, start=start, end=end,
            callSetIds=callSetIds))
        wormtableVariants = list(self.getWormtableVariants(
            [variantSetId], referenceName, start=start, end=end,
            callSetIds=callSetIds))
        self.assertEqual(len(ga4ghVariants), len(wormtableVariants))
        for ga4ghVariant, wormtableVariant in zip(
                ga4ghVariants, wormtableVariants):
            self.verifyVariantsEqual(ga4ghVariant, wormtableVariant)
            # Verify we've got the correct callSetIds.
            returnedCallSetIds = set(cs.callSetId for cs in ga4ghVariant.calls)
            self.assertEqual(len(returnedCallSetIds), len(ga4ghVariant.calls))
            self.assertEqual(returnedCallSetIds, set(callSetIds))

    def testSearchAllVariants(self):
        for variantSet in self.getVariantSets():
            for referenceName in self.getReferenceNames(variantSet.id):
                self.verifySearchVariants(
                    [variantSet.id], referenceName, 0, 2 ** 32)

    def testSearchVariantSlices(self):
        for variantSet in self.getVariantSets():
            for referenceName in self.getReferenceNames(variantSet.id):
                first, last = self.getStartBounds(variantSet.id, referenceName)
                mid = (last - first) // 2
                for pageSize in [1, 2, 3, 5, 1000]:
                    self.verifySearchVariants(
                        [variantSet.id], referenceName, first, last,
                        pageSize=pageSize)
                    self.verifySearchVariants(
                        [variantSet.id], referenceName, first, mid,
                        pageSize=pageSize)
                    self.verifySearchVariants(
                        [variantSet.id], referenceName, mid, last,
                        pageSize=pageSize)
                    self.verifySearchVariants(
                        [variantSet.id], referenceName, first, first + 1,
                        pageSize=pageSize)
                    self.verifySearchVariants(
                        [variantSet.id], referenceName, mid, mid + 1,
                        pageSize=pageSize)
                    self.verifySearchVariants(
                        [variantSet.id], referenceName, mid - 1, mid + 1,
                        pageSize=pageSize)
                    self.verifySearchVariants(
                        [variantSet.id], referenceName, last, last + 1,
                        pageSize=pageSize)

    def testSearchByCallSetIds(self):
        for variantSet in self.getVariantSets():
            callSetIds = self.getWormtableCallSetIds(variantSet.id)
            # limit the subsets we check over to some small value.
            for subset in utils.powerset(callSetIds, 20):
                for referenceName in self.getReferenceNames(variantSet.id):
                    self.verifySearchByCallSetIds(
                        variantSet.id, referenceName, 0, 2 ** 32, subset)

    def testUniqueIds(self):
        ids = set()
        for variantSet in self.getVariantSets():
            for referenceName in self.getReferenceNames(variantSet.id):
                for variant in self.getVariants(
                        [variantSet.id], referenceName):
                    self.assertTrue(variant.id not in ids)
                    ids.add(variant.id)

    def testAcrossVariantSets(self):
        allVariantSets = self.getVariantSetIds()
        for permLen in range(1, len(allVariantSets) + 1):
            for variantSets in itertools.permutations(allVariantSets, permLen):
                for referenceName in self.getCommonRefNames(variantSets):
                    for pageSize in [1, 2, 3, 5, 1000]:
                        self.verifySearchVariants(
                            variantSets, referenceName, 0, 2 ** 32,
                            pageSize=pageSize)


class TestCallSets(TestWormtableBackend):
    """
    Tests the searchCallSets end point.
    """

    def verifySearchCallSets(self, variantSetId, pageSize=1000):
        """
        Verifies callsets queries based on specified variant set ids.
        """
        ga4ghCallSets = list(self.getCallSets(variantSetId))
        wormtableCallSets = list(self.getWormtableCallSetIds(variantSetId))
        self.assertEqual(len(ga4ghCallSets), len(wormtableCallSets))
        for ga4ghCallSet, wormtableCallSet in zip(
                ga4ghCallSets, wormtableCallSets):
            self.verifyCallSetsEqual(
                variantSetId, ga4ghCallSet, wormtableCallSet)

    def verifyCallSetsEqual(
            self, variantSetId, ga4ghCallSet, wormtableCallSet):
        """
        Verifies that the ga4ghCallSet Ojbect is equal to the wormtable column.
        """
        self.assertEqual(ga4ghCallSet.id,
                         "{0}.{1}".format(variantSetId, wormtableCallSet))
        self.assertEqual(ga4ghCallSet.name, wormtableCallSet)
        self.assertEqual(ga4ghCallSet.sampleId, wormtableCallSet)

    def testSearchAllCallSets(self):
        for variantSet in self.getVariantSets():
            self.verifySearchCallSets(variantSet.id)

    def testSearchByCallSetIds(self):
        # TODO implement after schemas on name string search is determined
        pass

    def testCallSetIds(self):
        # verify that callSetIds are equal to sample columns in wormtable
        variantSetIdMap = self._backend.getVariantSetIdMap()
        for variantSetId in self.getVariantSetIds():
            callSetIds = set(self.getWormtableCallSetIds(variantSetId))
            variantSets = variantSetIdMap[variantSetId]
            sampleNames = set(variantSets.getSampleNames())
            self.assertEqual(callSetIds, sampleNames)
            self.assertIsNotNone(callSetIds)
            self.assertIsNotNone(sampleNames)

    def testUniqueCallSetIds(self):
        # verify that all callSetIds are unique
        for variantSetId in self.getVariantSetIds():
            callSetIds = set()
            for callSetId in self.getWormtableCallSetIds(variantSetId):
                self.assertTrue(callSetId not in callSetIds)
                callSetIds.add(callSetId)

    def testVariantSetIds(self):
        # test the VariantSetId field of GACallSet object is correct
        allCallSetIdMap = dict()
        for variantSetId in self.getVariantSetIds():
            for callSetId in self.getWormtableCallSetIds(variantSetId):
                if callSetId not in allCallSetIdMap:
                    allCallSetIdMap[callSetId] = []
                allCallSetIdMap[callSetId].append(variantSetId)
        wormtableCallSetIdMap = self._backend.getCallSetIdMap()
        # can't just call assertEqual on allCallSetIdMap and
        # wormtableCallSetIdMap since their values (arrays)
        # may have different orderings of the items
        self.assertEqual(
            sorted(allCallSetIdMap.keys()),
            sorted(wormtableCallSetIdMap.keys()))
        for key in allCallSetIdMap:
            for value in allCallSetIdMap[key]:
                self.assertIn(value, wormtableCallSetIdMap[key])
