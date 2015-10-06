"""
Data-driven tests for variants.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import hashlib

import vcf

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.datamodel.variants as variants
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
import tests.datadriven as datadriven
import tests.utils as utils


def testVariantSets():
    testDataDir = "tests/data/datasets/dataset1/variants"
    for test in datadriven.makeTests(testDataDir, VariantSetTest):
        yield test


class VariantSetTest(datadriven.DataDrivenTest):
    """
    Data driven test class for variant sets. Builds an alternative model of
    a variant set, and verifies that it is consistent with the model
    built by the variants.VariantSet object.
    """
    def __init__(self, variantSetId, baseDir):
        self._dataset = datasets.AbstractDataset("ds")
        super(VariantSetTest, self).__init__(variantSetId, baseDir)
        self._variantRecords = []
        self._referenceNames = set()
        # Read in all the VCF files in datadir and store each variant.
        for vcfFile in glob.glob(os.path.join(self._dataPath, "*.vcf.gz")):
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
            # When an END info tag is present it takes precedence
            if "END" in record.INFO:
                record.end = record.INFO["END"]
            self._variantRecords.append(record)

    def getDataModelInstance(self, localId, dataPath):
        return variants.HtslibVariantSet(
            self._dataset, localId, dataPath, None)

    def getProtocolClass(self):
        return protocol.VariantSet

    def _compareTwoListFloats(self, a, b):
        for ai, bi in zip(a, b):
            self.assertAlmostEqual(float(ai), float(bi), 5)

    def _verifyInfoEqual(self, gaObjectInfo, pyvcfInfo):
        def _assertEquivalentGaVCFValues(gaValue, pyvcfValue):
            if isinstance(pyvcfValue, str):
                self.assertEqual(gaValue, pyvcfValue)
            elif isinstance(pyvcfValue, (int, bool)):
                self.assertEqual(gaValue, str(pyvcfValue))
            elif isinstance(pyvcfValue, float):
                self.assertAlmostEqual(float(gaValue), float(pyvcfValue))
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
        # callSetId information is not available in pyvcf.model._Call
        self.assertEqual(gaCall.callSetName, pyvcfCall.sample)
        self.assertEqual(gaCall.genotype, genotype)
        phaseset = None
        if pyvcfCall.phased:
            phaseset = "*"
        self.assertEqual(gaCall.phaseset, phaseset)
        if len(gaCall.genotypeLikelihood) > 0:
            self._compareTwoListFloats(
                gaCall.genotypeLikelihood, pyvcfCall.data.GL)
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
                self.assertValid(protocol.Call, gaCall.toJsonDict())
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
            self.assertEqual(gaVariant.end, pyvcfVariant.end)
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

    def _verifyVariantsCallSetIds(self, searchsampleIds):
        """
        leaving searchSampleIds will get all samples.
        """
        gaCallSetVariants = []
        for referenceName in self._referenceNames:
            end = datamodel.PysamDatamodelMixin.vcfMax
            gaSearchIds = [
                str(datamodel.CallSetCompoundId(
                    self._gaObject.getCompoundId(), sampleId))
                for sampleId in searchsampleIds]
            gaVariants = list(self._gaObject.getVariants(
                referenceName, 0, end, gaSearchIds))
            self._verifyGaVariantsSample(gaVariants, searchsampleIds)
            gaCallSetVariants += gaVariants
            localVariants = filter(
                lambda v: v.CHROM == referenceName, self._variantRecords)
            self._verifyVariantsEqual(gaVariants, localVariants)
        self.assertEqual(len(gaCallSetVariants), len(self._variantRecords))

    def testSearchAllVariants(self):
        self._verifyVariantsCallSetIds(self.vcfSamples[:1])

    def testSearchCallSetIdsSystematic(self):
        for sampleIds in utils.powerset(self.vcfSamples, maxSets=10):
            self._verifyVariantsCallSetIds(list(sampleIds))

    def testVariantsValid(self):
        end = datamodel.PysamDatamodelMixin.vcfMax
        for referenceName in self._referenceNames:
            iterator = self._gaObject.getVariants(
                referenceName, 0, end)
            for gaVariant in iterator:
                self.assertValid(protocol.Variant, gaVariant.toJsonDict())

    def _getPyvcfVariants(
            self, referenceName, startPosition=0, endPosition=2**30):
        """
        variants with in interval [startPosition, endPosition)
        """
        localVariants = []
        for variant in self._variantRecords:
            case1 = (variant.start <= startPosition and
                     variant.end > startPosition)
            case2 = (variant.start > startPosition and
                     variant.start < endPosition)
            if (variant.CHROM == referenceName) and (case1 or case2):
                localVariants.append(variant)
        return localVariants

    def _assertEmptyVariant(self, referenceName, startPosition, endPosition):
        gaVariants = list(self._gaObject.getVariants(
            referenceName, startPosition, endPosition, []))
        self.assertEqual(len(gaVariants), 0)
        pyvcfVariants = self._getPyvcfVariants(
            referenceName, startPosition, endPosition)
        self.assertEqual(len(pyvcfVariants), 0)

    def _assertVariantsEqualInRange(
            self, referenceName, startPosition, endPosition):
        gaVariants = list(self._gaObject.getVariants(
            referenceName, startPosition, endPosition, []))
        pyvcfVariants = self._getPyvcfVariants(
            referenceName, startPosition, endPosition)
        self.assertGreaterEqual(len(gaVariants), 0)
        self._verifyVariantsEqual(gaVariants, pyvcfVariants)

    def testVariantInSegments(self):
        for referenceName in self._referenceNames:
            localVariants = self._getPyvcfVariants(referenceName)
            mini = localVariants[0].start
            # NOTE, the end of the last variant may not reflect the END
            maxi = max([v.end for v in localVariants])
            seglen = int(maxi-mini) // 3
            seg1 = mini + seglen
            seg2 = seg1 + seglen
            self._assertEmptyVariant(referenceName, -1, mini)
            self._assertEmptyVariant(referenceName, maxi, maxi+100)
            self._assertVariantsEqualInRange(referenceName, mini, seg1)
            self._assertVariantsEqualInRange(referenceName, seg1, seg2)
            self._assertVariantsEqualInRange(referenceName, seg2, maxi)

    def _gaVariantEqualsPyvcfVariant(self, gaVariant, pyvcfVariant):
        if gaVariant.start != pyvcfVariant.start:
            return False
        if gaVariant.end != pyvcfVariant.end:
            return False
        alt = pyvcfVariant.ALT
        if len(alt) == 1 and alt[0] is None:
            alt = []
        if not pyvcfVariant.is_sv and gaVariant.alternateBases != alt:
            return False
        self._verifyVariantsEqual([gaVariant], [pyvcfVariant])
        return True

    def _pyvcfVariantIsInGaVarants(
            self, pyvcfVariant, intervalStart, intervalEnd):
        isIn = False
        gaVariants = list(self._gaObject.getVariants(
            pyvcfVariant.CHROM, intervalStart, intervalEnd, []))
        for gaVariant in gaVariants:
            if self._gaVariantEqualsPyvcfVariant(gaVariant, pyvcfVariant):
                isIn = True
                break
        return isIn

    def testVariantFromToEveryVariant(self):
        for variant in self._variantRecords:
            variantStart = variant.start
            variantEnd = variant.end
            # Cases of search interval does not overlap with variant
            # interval on the left

            self.assertFalse(self._pyvcfVariantIsInGaVarants(
                variant, variantStart-2, variantStart-1))
            self.assertFalse(self._pyvcfVariantIsInGaVarants(
                variant, variantStart-1, variantStart))
            # interval on the right
            self.assertFalse(self._pyvcfVariantIsInGaVarants(
                variant, variantEnd, variantEnd+1))
            self.assertFalse(self._pyvcfVariantIsInGaVarants(
                variant, variantEnd+1, variantEnd+2))
            # case of search interval is within variant

            self.assertTrue(self._pyvcfVariantIsInGaVarants(
                variant, variantStart, variantEnd))
            if (variantEnd - variantStart) != 1:
                self.assertTrue(self._pyvcfVariantIsInGaVarants(
                    variant, variantStart, variantEnd-1))
                self.assertTrue(self._pyvcfVariantIsInGaVarants(
                    variant, variantStart+1, variantEnd))

            # case of search interval contains variant
            self.assertTrue(self._pyvcfVariantIsInGaVarants(
                variant, variantStart-1, variantEnd+1))
            # cases of search interval intersec with variant
            self.assertTrue(self._pyvcfVariantIsInGaVarants(
                variant, variantStart-1, variantStart+1))
            self.assertTrue(self._pyvcfVariantIsInGaVarants(
                variant, variantStart, variantStart+1))
            self.assertTrue(self._pyvcfVariantIsInGaVarants(
                variant, variantEnd-1, variantEnd))
            self.assertTrue(self._pyvcfVariantIsInGaVarants(
                variant, variantEnd-1, variantEnd+1))

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
        for prefix, content in [("FORMAT", self._formats),
                                ("INFO", self._infos)]:
            for contentKey in content.keys():
                key = "{0}.{1}".format(prefix, contentKey)
                if key == "FORMAT.GT":
                    gtCounter += 1
                else:
                    self.assertEqual(
                        keyMap[key].type, content[contentKey].type)
                    self.assertEqual(keyMap[key].number, convertPyvcfNumber(
                        content[contentKey].num))
                    self.assertEqual(
                        keyMap[key].description, content[contentKey].desc)
        testMetaLength = (
            1 + len(self._formats) + len(self._infos) - gtCounter)
        self.assertEqual(len(keyMap), testMetaLength)

    def testGetVariantsCallSets(self):
        variantSet = self._gaObject
        start = 0
        end = datamodel.PysamDatamodelMixin.vcfMax
        callSetIds = [callSet.getId() for callSet in variantSet.getCallSets()]
        someCallSetIds = callSetIds[0:3]
        for callSetId in callSetIds:
            returnedCallSet = variantSet.getCallSet(callSetId)
            self.assertEqual(callSetId, returnedCallSet.getId())
        for referenceName in self._referenceNames:
            # passing None as the callSetIds argument should be equivalent
            # to passing all of the possible callSetIds as an argument
            noneRecords = list(variantSet.getVariants(
                referenceName, start, end, None))
            allRecords = list(variantSet.getVariants(
                referenceName, start, end, callSetIds))
            self.assertEqual(len(noneRecords), len(allRecords))
            for noneRecord, allRecord in zip(noneRecords, allRecords):
                for noneCall, allCall in zip(
                        noneRecord.calls, allRecord.calls):
                    self.assertEqual(
                        noneCall.callSetName, allCall.callSetName)

            # passing an empty list as the callSetIds argument should
            # return no callsets for any variant
            emptyRecords = variantSet.getVariants(
                referenceName, start, end, [])
            for record in emptyRecords:
                self.assertEqual(len(record.calls), 0)

            # passing some callSetIds as the callSetIds argument should
            # return only those calls
            someRecords = list(variantSet.getVariants(
                referenceName, start, end, someCallSetIds))
            for record in someRecords:
                self.assertEqual(len(record.calls), len(someCallSetIds))
                for call, someId in zip(record.calls, someCallSetIds):
                    self.assertEqual(call.callSetId, someId)

    def testGetVariant(self):
        variantSet = self._gaObject
        for referenceName in self._referenceNames:
            refnameVariants = self._getPyvcfVariants(referenceName)
            for variant in refnameVariants:
                # positive test: get the expected variant
                md5 = self._hashVariant(variant)
                compoundId = datamodel.VariantCompoundId(
                    variantSet.getCompoundId(), referenceName,
                    variant.start, md5)
                gotVariant = variantSet.getVariant(compoundId)
                self.assertEqual(str(compoundId), gotVariant.id)

                # negative test: change start position to past variant
                wrongStart = variant.end
                compoundId = datamodel.VariantCompoundId(
                    variantSet.getCompoundId(), referenceName,
                    wrongStart, md5)
                try:
                    gotVariant = variantSet.getVariant(compoundId)
                    self.assertNotEqual(variant.start, gotVariant.start)
                except exceptions.ObjectNotFoundException:
                    pass

                # negative test: change reference name
                compoundId = datamodel.VariantCompoundId(
                    variantSet.getCompoundId(), "wrong reference name",
                    variant.start, md5)
                with self.assertRaises(exceptions.ObjectNotFoundException):
                    variantSet.getVariant(compoundId)

                # negative test: change hash
                compoundId = datamodel.VariantCompoundId(
                    variantSet.getCompoundId(), referenceName,
                    variant.start, "wrong hash")
                with self.assertRaises(exceptions.ObjectNotFoundException):
                    variantSet.getVariant(compoundId)

    def _hashVariant(self, record):
        if record.ALT[0] is None:
            alts = tuple()
        else:
            alts = tuple([str(substitution) for substitution in record.ALT])
        return hashlib.md5(record.REF + str(alts)).hexdigest()
