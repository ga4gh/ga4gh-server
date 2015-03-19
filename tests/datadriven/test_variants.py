"""
Data-driven tests for variants.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import vcf

import ga4gh.protocol as protocol
import ga4gh.datamodel.variants as variants
import tests.datadriven as datadriven
import tests.utils as utils


def testVariantSets():
    testDataDir = "tests/data/variants"
    for test in datadriven.makeTests(testDataDir, VariantSetTest):
        yield test


class VariantSetTest(datadriven.DataDrivenTest):
    """
    Data driven test class for variant sets. Builds an alternative model of
    a variant set, and verifies that it is consistent with the model
    built by the variants.VariantSet object.
    """

    def __init__(self, variantSetId, baseDir):
        super(VariantSetTest, self).__init__(variantSetId, baseDir)
        self._variantRecords = []
        self._referenceNames = set()
        # Read in all the VCF files in datadir and store each variant.
        for vcfFile in glob.glob(os.path.join(self._dataDir, "*.vcf.gz")):
            self._readVcf(vcfFile)

    def _readVcf(self, vcfFileName):
        """
        Reads all variants and metadata from the specified VCF file and
        store locally.
        """
        vcfReader = vcf.Reader(filename=vcfFileName)
        self.vcfSamples = vcfReader.samples
        metadata = vcfReader.metadata
        self._vcfVersion = metadata["fileformat"]
        self._infos = vcfReader.infos
        self._formats = vcfReader.formats
        for record in vcfReader:
            self._referenceNames.add(record.CHROM)
            self._variantRecords.append(record)

    def getDataModelClass(self):
        return variants.HtslibVariantSet

    def getProtocolClass(self):
        return protocol.GAVariantSet

    def _verifyInfoEqual(self, gaObjectInfo, pyvcfInfo):
        def _compareTwoFloats(a, b):
            afloat = float(a)
            bfloat = float(b)
            return (abs(afloat - bfloat) <= 0.001 * bfloat)

        def _compareGaValuePyvcfValue(gaValue, pyvcfValue):
            if type(pyvcfValue) is str:
                self.assertEqual(gaValue, pyvcfValue)
            elif type(pyvcfValue) is int or type(value) is bool:
                self.assertEqual(gaValue, str(pyvcfValue))
            elif type(pyvcfValue) is float:
                _compareTwoFloats(gaValue, pyvcfValue)
            elif pyvcfValue is None:
                self.assertEqual(gaValue, ".")
            else:
                print (gaValue, pyvcfValue)

        for key, value in pyvcfInfo.iteritems():
            if type(value) is list:
                self.assertEqual(len(gaObjectInfo[key]), len(value))
                for gaValue, pyvcfValue in (
                  zip(gaObjectInfo[key], pyvcfInfo[key])):
                    _compareGaValuePyvcfValue(gaValue, pyvcfValue)
            else:
                self.assertEqual(len(gaObjectInfo[key]), 1)

    def _verifyVariantsEqual(self, gaVariants, pyvcfVariants):
        """
        Verifies that the lists of GA4GH variants and pyvcf variants
        are equivalent.
        """
        def _verifyVariantCallEqual(gaCall, pyvcfCall, pyvcfInfo):
            genotype, phaseset = variants.convertVCFGenotype(
                    pyvcfCall.data.GT, pyvcfCall.phased)
            self.assertEqual(gaCall.callSetId, pyvcfCall.site.ID)
            self.assertEqual(gaCall.callSetName, pyvcfCall.sample)
            self.assertEqual(gaCall.genotype, genotype)
            # TODO: Need to check the phaseset!
            # gaCall.phaseset is currently not implemented?
            # self.assertEqual(gaCall.phaseset,phaseset)
            if len(gaCall.genotypeLikelihood) > 0:
                gaCallgenotypeLK = [(
                    '%.2f' % elem) for elem in gaCall.genotypeLikelihood]
                pyvcfCallgenotypeLK = [(
                    '%.2f' % elem) for elem in pyvcfCall.data.GL]
                self.assertEqual(gaCallgenotypeLK, pyvcfCallgenotypeLK)
            else:
                self.assertNotIn("GL", pyvcfCall.data)
            self._verifyInfoEqual(gaCall.info, pyvcfInfo)

        def _verifyVariantCalls(gaVariantCalls, pyvcfCallMap, pyvcfInfo):
            for gaCall in gaVariantCalls:
                self.assertTrue(protocol.GACall.validate(
                    gaCall.toJsonDict()))
                self.assertIn(gaCall.callSetName, pyvcfCallMap)
                pyvcfCall = pyvcfCallMap[gaCall.callSetName]
                _verifyVariantCallEqual(gaCall, pyvcfCall, pyvcfInfo)

        self.assertEqual(len(gaVariants), len(pyvcfVariants))
        for gaVariant, pyvcfVariant in zip(gaVariants, pyvcfVariants):
            pyvcfInfo = pyvcfVariant.INFO
            self.assertEqual(gaVariant.referenceName, pyvcfVariant.CHROM)
            self.assertEqual(gaVariant.referenceBases, pyvcfVariant.REF)
            # pyvcf uses 1-based indexing.
            self.assertEqual(gaVariant.start, pyvcfVariant.POS - 1)
            # When an END info tag is present it takes precedence
            end = pyvcfVariant.end
            if "END" in pyvcfInfo:
                end = pyvcfInfo["END"]
            self.assertEqual(gaVariant.end, end)
            alt = pyvcfVariant.ALT
            # PyVCF does something funny when no ALT allele is provided.
            # TODO we should clarify exactly what this means.
            if len(alt) == 1 and alt[0] is None:
                alt = []
            if pyvcfVariant.is_sv:
                self.assertEqual(len(alt), len(gaVariant.alternateBases))
                for alt1, alt2 in zip(alt, gaVariant.alternateBases):
                    self.assertEqual(str(alt1), str(alt2))
            else:
                self.assertEqual(gaVariant.alternateBases, alt)

            pyvcfCallMap = {}
            for call in pyvcfVariant:
                pyvcfCallMap[call.sample] = call
            _verifyVariantCalls(gaVariant.calls, pyvcfCallMap, pyvcfInfo)
            self._verifyInfoEqual(gaVariant.info, pyvcfInfo)

    def _verifyGaVariantsSample(self, gaVariants, sampleIds):
        for v in gaVariants:
            self.assertEqual(len(v.calls), len(sampleIds))
            for call in v.calls:
                self.assertIn(call.callSetName, sampleIds)

    def _verifyVariantsCallSetIds(self, searchVariants, searchsampleIds):
        """
        SearchVariants = None, this will get all variants
        If searchsampleIds = None, this will get all samples
        else if = [of something], this will iterate through the list
        """
        def filterLocalVariant(localVariant):
            updatedList = []
            for record in localVariants:
                if record.ID is not None and record.ID in searchVariants:
                    updatedList.append(record)
            return updatedList

        gaCallSetVariants = []
        for referenceName in self._referenceNames:
            end = 2**30  # TODO This is arbitrary, and pysam can choke. FIX!
            gaVariants = list(self._gaObject.getVariants(
                referenceName, 0, end, searchVariants, searchsampleIds))
            self._verifyGaVariantsSample(gaVariants, searchsampleIds)
            gaCallSetVariants += gaVariants
            localVariants = [
                v for v in self._variantRecords if (v.CHROM == referenceName)]
            if searchVariants is not None:
                localVariants = filterLocalVariant(localVariants)
            self._verifyVariantsEqual(gaVariants, localVariants)
        if searchVariants is None:
            self.assertEqual(len(gaCallSetVariants), len(self._variantRecords))

    def testSearchAllVariants(self):
        self._verifyVariantsCallSetIds(None, [])

    def testSearchCallSetIdsSystematic(self):
        for sampleIds in utils.powerset(self.vcfSamples, maxSets=10):
            self._verifyVariantsCallSetIds(None, list(sampleIds))

    def testVariantsValid(self):
        # NOTE: time consuming
        end = 2**30  # TODO This is arbitrary, and pysam can choke. FIX!
        for referenceName in self._referenceNames:
            iterator = self._gaObject.getVariants(
                referenceName, 0, end, None, None)
            for gaVariant in iterator:
                self.assertTrue(protocol.GAVariant.validate(
                    gaVariant.toJsonDict()))

    def testVariantSetMetadata(self):
        def convertPyvcfNumber(number):
            if number == -1:
                ret = "A"
            elif number == -2:
                ret = "G"
            elif number is None:
                ret = "."
            else:
                ret = str(number)
            return ret

        keyMap = {}
        for metadata in self._gaObject.getMetadata():
            keyMap[metadata.key] = metadata

        metadata = keyMap["version"]
        self.assertEqual(metadata.value, self._vcfVersion)

        gtCounter = 0
        for infoKey in self._infos.keys():
            key = "INFO.{}".format(infoKey)
            self.assertEqual(
                keyMap[key].type, self._infos[infoKey].type)
            self.assertEqual(keyMap[key].number, convertPyvcfNumber(
                self._infos[infoKey].num))

        for formatKey in self._formats.keys():
            if formatKey == "GT":
                gtCounter += 1
            else:
                key = "FORMAT.{}".format(formatKey)
                self.assertEqual(
                    keyMap[key].type, self._formats[formatKey].type)
                self.assertEqual(keyMap[key].number, convertPyvcfNumber(
                    self._formats[formatKey].num))

        testMetaLength = (
            1 + len(self._formats) + len(self._infos) - gtCounter)
        self.assertEqual(len(keyMap), testMetaLength)
