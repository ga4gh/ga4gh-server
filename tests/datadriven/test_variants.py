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
        metadata = vcfReader.metadata
        self._vcfVersion = metadata["fileformat"]
        self._infos = vcfReader.infos
        self._formats = vcfReader.formats
        self.vcfSamples = vcfReader.samples
        for record in vcfReader:
            self._referenceNames.add(record.CHROM)
            self._variantRecords.append(record)

    def getDataModelClass(self):
        return variants.HtslibVariantSet

    def getProtocolClass(self):
        return protocol.GAVariantSet

    def _floatsAgreeWithinTolerance(self, a, b, decimalPoints=7):
        afloat = float(a)
        bfloat = float(b)
        return (abs(afloat - bfloat) <= 10**(-decimalPoints) * abs(bfloat))

    def _compareTwoListFloats(self, a, b):
        for ai, bi in zip(a, b):
            if not self._floatsAgreeWithinTolerance(ai, bi):
                return False
        return True

    def _verifyInfoEqual(self, gaObjectInfo, pyvcfInfo):
        def _assertEquivalentGaVCFValues(gaValue, pyvcfValue):
            if isinstance(pyvcfValue, str):
                self.assertEqual(gaValue, pyvcfValue)
            elif isinstance(pyvcfValue, (int, bool)):
                self.assertEqual(gaValue, str(pyvcfValue))
            elif isinstance(pyvcfValue, float):
                self._floatsAgreeWithinTolerance(gaValue, pyvcfValue)
            elif pyvcfValue is None:
                self.assertEqual(gaValue, ".")
            else:
                raise Exception(key, (
                    " values are inconsistent",
                    "between ga4ghObject and pyvcf!"))

        for key, value in pyvcfInfo.iteritems():
            if isinstance(value, list):
                self.assertEqual(len(gaObjectInfo[key]), len(value))
                for gaValue, pyvcfValue in zip(
                  gaObjectInfo[key], pyvcfInfo[key]):
                    _assertEquivalentGaVCFValues(gaValue, pyvcfValue)
            else:
                self.assertEqual(len(gaObjectInfo[key]), 1)

    def _verifyVariantCallEqual(self, gaCall, pyvcfCall):
        genotype, phaseset = variants.convertVCFGenotype(
                pyvcfCall.data.GT, pyvcfCall.phased)
        self.assertEqual(gaCall.callSetId, pyvcfCall.site.ID)
        self.assertEqual(gaCall.callSetName, pyvcfCall.sample)
        self.assertEqual(gaCall.genotype, genotype)
        # TODO: Need to check the phaseset!
        # gaCall.phaseset is currently not implemented?
        # self.assertEqual(gaCall.phaseset,phaseset)
        if len(gaCall.genotypeLikelihood) > 0:
            self.assertTrue(self._compareTwoListFloats(
                gaCall.genotypeLikelihood, pyvcfCall.data.GL))
        else:
            self.assertNotIn("GL", pyvcfCall.data)
        for key, value in gaCall.info.items():
            if key != "GT" and key != "GL":
                if isinstance(value[0], (list, tuple)):
                    self._compareTwoListFloats(value[0], getattr(
                        pyvcfCall.data, key))
                elif isinstance(value[0], float):
                    self._compareTwoFloats(value[0], getattr(
                        pyvcfCall.data, key))

    def _verifyVariantsEqual(self, gaVariants, pyvcfVariants):
        """
        Verifies that the lists of GA4GH variants and pyvcf variants
        are equivalent.
        """
        def _verifyVariantCalls():
            for gaCall in gaVariant.calls:
                self.assertTrue(protocol.GACall.validate(
                    gaCall.toJsonDict()))
                self.assertIn(gaCall.callSetName, pyvcfCallMap)
                pyvcfCall = pyvcfCallMap[gaCall.callSetName]
                self._verifyVariantCallEqual(gaCall, pyvcfCall)

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
            self._verifyInfoEqual(gaVariant.info, pyvcfInfo)
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
            _verifyVariantCalls()

    def _verifyGaVariantsSample(self, gaVariants, sampleIds):
        for variant in gaVariants:
            self.assertEqual(len(variant.calls), len(sampleIds))
            for call in variant.calls:
                self.assertIn(call.callSetName, sampleIds)

    def _verifyVariantsCallSetIds(self, searchVariants, searchsampleIds):
        """
        Leaving searchVariants empty will get all variants,
        leaving searchSampleIds will get all samples.
        """
        gaCallSetVariants = []
        for referenceName in self._referenceNames:
            end = 2**30  # TODO This is arbitrary, and pysam can choke. FIX!
            gaVariants = list(self._gaObject.getVariants(
                referenceName, 0, end, searchVariants, searchsampleIds))
            self._verifyGaVariantsSample(gaVariants, searchsampleIds)
            gaCallSetVariants += gaVariants
            localVariants = filter(
                lambda v: v.CHROM == referenceName, self._variantRecords)
            self._verifyVariantsEqual(gaVariants, localVariants)
        if searchVariants is None:
            self.assertEqual(len(gaCallSetVariants), len(self._variantRecords))

    def testSearchAllVariants(self):
        self._verifyVariantsCallSetIds(None, [])

    def testSearchCallSetIdsSystematic(self):
        for sampleIds in utils.powerset(self.vcfSamples, maxSets=10):
            self._verifyVariantsCallSetIds(None, list(sampleIds))

    def testVariantsValid(self):
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
