"""
Tests for the backend objects. We instantiate local copies of
the backends and invoke the entry points for the protocol methods.
We do not set up any server processes or communicate over sockets.
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

import pysam
import wormtable as wt

import tests.utils as utils
import ga4gh.backend as backend
import ga4gh.protocol as protocol


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


class TestAbstractBackend(unittest.TestCase):
    """
    Provides testing harness for testing methods in AbstractBackend,
    using an instance of the mock SimulatedBackend object.
    """
    def setUp(self):
        self._backend = backend.SimulatedBackend(
            numCalls=100, numVariantSets=10)
        # TODO arbitrary values, pepper to taste

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
        Returns an iterator over the variantSets, abstracting away
        the details of the pageSize.
        """
        request = protocol.GASearchVariantSetsRequest()
        return self.resultIterator(
            request, pageSize, self._backend.searchVariantSets,
            protocol.GASearchVariantSetsResponse, "variantSets")

    def getVariants(
            self, variantSetIds, referenceName, start=0, end=2 ** 32,
            pageSize=100, callSetIds=None):
        """
        Returns an iterator over the specified list of variants,
        abstracting out paging details.
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
        Returns an iterator over the callsets in a specified
        variant set.
        """
        request = protocol.GASearchCallSetsRequest()
        request.variantSetIds = [variantSetId]
        return self.resultIterator(
            request, pageSize, self._backend.searchCallSets,
            protocol.GASearchCallSetsResponse, "callSets")

    def testGetVariantSets(self):
        sortedVariantSetsFromGetter = sorted(self._backend.getVariantSets())
        sortedVariantSetMapValues = sorted(
            self._backend._variantSetIdMap.values())
        self.assertEqual(
            sortedVariantSetMapValues, sortedVariantSetsFromGetter)

    def testGetVariantSetIdMap(self):
        variantSetMapFromGetter = self._backend.getVariantSetIdMap()
        variantSetMap = self._backend._variantSetIdMap
        self.assertEqual(variantSetMapFromGetter, variantSetMap)

    def testParsePageToken(self):
        goodPageToken = "12:34:567:8:9000"
        parsedToken = self._backend.parsePageToken(goodPageToken, 5)
        self.assertEqual(parsedToken[2], 567)

    def testRunSearchRequest(self):
        request = protocol.GASearchVariantSetsRequest()
        responseStr = self._backend.runSearchRequest(
            request.toJsonString(), protocol.GASearchVariantSetsRequest,
            protocol.GASearchVariantSetsResponse,
            self._backend.variantSetsGenerator)
        response = protocol.GASearchVariantSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.GASearchVariantSetsResponse))

    def testSearchVariantSets(self):
        request = protocol.GASearchVariantSetsRequest()
        responseStr = self._backend.searchVariantSets(request.toJsonString())
        response = protocol.GASearchVariantSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.GASearchVariantSetsResponse))

    def testSearchVariants(self):
        request = protocol.GASearchVariantsRequest()
        responseStr = self._backend.searchVariants(request.toJsonString())
        response = protocol.GASearchVariantsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.GASearchVariantsResponse))

    def testSearchCallSets(self):
        variantSetIds = [
            variantSet.id for variantSet in self.getVariantSets(pageSize=1)]
        request = protocol.GASearchCallSetsRequest()
        request.variantSetIds = variantSetIds[:1]
        responseStr = self._backend.searchCallSets(request.toJsonString())
        response = protocol.GASearchCallSetsResponse.fromJsonString(
            responseStr)
        self.assertTrue(
            isinstance(response, protocol.GASearchCallSetsResponse))

    def testVariantSetPagination(self):
        results = []
        for pageSize in range(1, 100):
            variantSetIds = [
                variantSet.id for variantSet in self.getVariantSets(
                    pageSize=pageSize)]
            results.append(variantSetIds)
        for result in results[1:]:
            self.assertEqual(result, results[0])


class TestFileSystemBackend(TestAbstractBackend):
    """
    Tests proper initialization of the filesystem backend using indexed
    files in the tests/data directory.
    """
    def setUp(self):
        self._dataDir = os.path.join("tests", "data")
        self._variantsDir = os.path.join(self._dataDir, "variants")
        self._vcfs = {}
        self._variants = []
        self._referenceNames = set()
        self._chromFileMap = {}
        for relativePath in os.listdir(self._variantsDir):
            pathToFiles = os.path.join(self._variantsDir, relativePath)
            self._vcfs[relativePath] = []
            for vcfFile in glob.glob(os.path.join(
                    pathToFiles, "*.vcf.gz")):
                self._chromFileMap[relativePath] = {}
                self._vcfs[relativePath].append(vcfFile)
                vcf = pysam.VariantFile(filename=vcfFile)
                for chrom in vcf.index:
                    self._chromFileMap[relativePath][chrom] = vcf
        self._backend = backend.FileSystemBackend(self._dataDir)

    def testVariantSetIds(self):
        variantSets = [variantSet for variantSet in self.getVariantSets()]
        self.assertEqual(len(variantSets), len(self._vcfs))
        ids = set(variantSet.id for variantSet in variantSets)
        self.assertEqual(ids, set(self._vcfs.keys()))


# Everything below here is wormtable related, can be removed on resolving
# Issue #210 (Removing wormtable dependency)
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


class TestWormtableBackend(TestAbstractBackend):

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

    def getVariantSetIds(self):
        """
        Return all variantSetIds.
        """
        return [variantSet.id for variantSet in self.getVariantSets()]

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


class TestWormtableVariantSets(TestWormtableBackend):
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


class TestWormtableVariants(TestWormtableBackend):
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


class TestWormtableCallSets(TestWormtableBackend):
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

    def testSearchAllCallSets(self):
        for variantSet in self.getVariantSets():
            self.verifySearchCallSets(variantSet.id)

    def testSearchByCallSetIds(self):
        # TODO implement after schemas on name string search is determined
        pass

    def testUniqueCallSetIds(self):
        # verify that all callSetIds are unique
        for variantSetId in self.getVariantSetIds():
            callSetIds = set()
            for callSetId in self.getWormtableCallSetIds(variantSetId):
                self.assertTrue(callSetId not in callSetIds)
                callSetIds.add(callSetId)
