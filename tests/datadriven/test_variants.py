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

import ga4gh.server.datamodel as datamodel
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.datamodel.references as references
import ga4gh.server.datamodel.variants as variants
import ga4gh.server.exceptions as exceptions
import tests.datadriven as datadriven
import tests.paths as paths

import ga4gh.common.utils as utils
import ga4gh.schemas.protocol as protocol


def testVariantSets():
    testDataDir = os.path.join(
        paths.testDataDir, "datasets/dataset1/variants")
    for test in datadriven.makeTests(testDataDir, VariantSetTest):
        yield test


def convertVCFPhaseset(vcfPhaseset):
    """
    Parses the VCF phaseset string
    """
    if vcfPhaseset is not None and vcfPhaseset != ".":
        phaseset = vcfPhaseset
    else:
        phaseset = "*"
    return phaseset


def convertVCFGenotype(vcfGenotype):
    """
    Parses the VCF genotype
    """
    if vcfGenotype is not None:
        delim = "/"
        if "|" in vcfGenotype:
            delim = "|"
        if "." in vcfGenotype:
            genotype = [-1]
        else:
            genotype = map(int, vcfGenotype.split(delim))
    else:
        genotype = [-1]
    return genotype


class VariantSetTest(datadriven.DataDrivenTest):
    """
    Data driven test class for variant sets. Builds an alternative model of
    a variant set, and verifies that it is consistent with the model
    built by the variants.VariantSet object.
    """
    def __init__(self, variantSetId, baseDir):
        self._dataset = datasets.Dataset("ds")
        super(VariantSetTest, self).__init__(variantSetId, baseDir)
        self._variantRecords = []
        self._reference_names = set()
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
        self._filters = vcfReader.filters
        self.vcfSamples = vcfReader.samples
        for record in vcfReader:
            self._reference_names.add(record.CHROM)
            # When an END info tag is present it takes precedence
            if "END" in record.INFO:
                record.end = record.INFO["END"]
            self._variantRecords.append(record)

    def getDataModelInstance(self, localId, dataPath):
        variantSet = variants.HtslibVariantSet(self._dataset, localId)
        variantSet.populateFromDirectory(dataPath)
        referenceSet = references.AbstractReferenceSet("test")
        variantSet.setReferenceSet(referenceSet)
        return variantSet

    def getProtocolClass(self):
        return protocol.VariantSet

    def _compareTwoListFloats(self, a, b):
        for ai, bi in zip(a, b):
            self.assertAlmostEqual(float(ai), float(bi), 5)

    def _verifyInfoEqual(self, gaObjectInfo, pyvcfInfo):
        def _assertEquivalentGaVCFValues(gaValue, pyvcfValue):
            compareValue = protocol.getValueFromValue(gaValue)
            if isinstance(pyvcfValue, str):
                self.assertEqual(compareValue, pyvcfValue)
            elif isinstance(pyvcfValue, (int, bool)):
                self.assertEqual(compareValue, pyvcfValue)
            elif isinstance(pyvcfValue, float):
                self.assertAlmostEqual(float(compareValue), float(pyvcfValue))
            elif pyvcfValue is None:
                self.assertEqual(compareValue, ".")
            else:
                raise Exception(key, (
                    " values are inconsistent",
                    "between ga4ghObject and pyvcf!"))

        for key, value in pyvcfInfo.iteritems():
            if isinstance(value, list):
                self.assertEqual(len(gaObjectInfo[key].values), len(value))
                for gaValue, pyvcfValue in zip(
                        gaObjectInfo[key].values, pyvcfInfo[key]):
                    _assertEquivalentGaVCFValues(gaValue, pyvcfValue)
            else:
                self.assertEqual(len(gaObjectInfo[key].values), 1)

    def _verifyStructuralInfo(self, gaVariant, pyvcfInfo):
        pyvcfType = None
        pyvcfLen = None
        pyvcfPos = None
        pyvcfEnd = None
        for key, value in pyvcfInfo.iteritems():
            if key == 'SVTYPE':
                pyvcfType = value
            elif key == 'SVLEN':
                pyvcfLen = value[0]
            elif key == 'CIPOS':
                pyvcfPos = value
            elif key == 'CIEND':
                pyvcfEnd = value
        if gaVariant.variant_type is not None or pyvcfType is not None:
            self.assertEqual(gaVariant.variant_type, pyvcfType)
        if gaVariant.svlen != 0 or pyvcfLen is not None:
            self.assertEqual(gaVariant.svlen, pyvcfLen)
        if len(gaVariant.cipos) != 0 or pyvcfPos is not None:
            self.assertEqual(gaVariant.cipos, pyvcfPos)
        if len(gaVariant.ciend) != 0 or pyvcfEnd is not None:
            self.assertEqual(gaVariant.ciend, pyvcfEnd)

    def _verifyVariantCallEqual(self, gaCall, pyvcfCall):
        genotype = convertVCFGenotype(pyvcfCall.data.GT)
        # callSetId information is not available in pyvcf.model._Call
        self.assertEqual(gaCall.call_set_name, pyvcfCall.sample)
        self.assertEqual(
            [x.number_value for x in gaCall.genotype.values], genotype)
        if len(pyvcfCall.data.GT.split("/")) == 1:
            # corner case when there is only a single genotype pyvcf
            # and pysam disagree
            self.assertTrue(gaCall.phaseset)
        else:
            # Compare our value for phased with pyvcf
            phaseset = None
            if pyvcfCall.phased:
                phaseset = str(pyvcfCall.phased)
            self.assertEqual(gaCall.phaseset, phaseset)
        if len(gaCall.genotype_likelihood) > 0:
            self._compareTwoListFloats(
                gaCall.genotype_likelihood, pyvcfCall.data.GL)
        else:
            self.assertNotIn("GL", pyvcfCall.data)
        for key, value in gaCall.attributes.attr.items():
            if key != "GT" and key != "GL":
                if isinstance(value.values[0], (list, tuple)):
                    self._compareTwoListFloats(value.values[0], getattr(
                        pyvcfCall.data, key))
                elif isinstance(value.values[0], float):
                    self._compareTwoFloats(value.values[0], getattr(
                        pyvcfCall.data, key))

    def _verifyVariantsEqual(self, gaVariants, pyvcfVariants):
        """
        Verifies that the lists of GA4GH variants and pyvcf variants
        are equivalent.
        """
        def _verifyVariantCalls():
            for gaCall in gaVariant.calls:
                self.assertValid(protocol.Call, protocol.toJson(gaCall))
                self.assertIn(gaCall.call_set_name, pyvcfCallMap)
                pyvcfCall = pyvcfCallMap[gaCall.call_set_name]
                self._verifyVariantCallEqual(gaCall, pyvcfCall)

        self.assertEqual(len(gaVariants), len(pyvcfVariants))
        for gaVariant, pyvcfVariant in zip(gaVariants, pyvcfVariants):
            pyvcfInfo = pyvcfVariant.INFO
            self.assertEqual(gaVariant.reference_name, pyvcfVariant.CHROM)
            self.assertEqual(gaVariant.reference_bases, pyvcfVariant.REF)
            # pyvcf uses 1-based indexing.
            self.assertEqual(gaVariant.start, pyvcfVariant.POS - 1)
            self.assertEqual(gaVariant.end, pyvcfVariant.end)
            self._verifyInfoEqual(gaVariant.attributes.attr, pyvcfInfo)
            alt = pyvcfVariant.ALT
            # PyVCF does something funny when no ALT allele is provided.
            # TODO we should clarify exactly what this means.
            if len(alt) == 1 and alt[0] is None:
                alt = []
            if pyvcfVariant.is_sv:
                self.assertEqual(len(alt), len(gaVariant.alternate_bases))
                for alt1, alt2 in zip(alt, gaVariant.alternate_bases):
                    self.assertEqual(str(alt1), str(alt2))
            else:
                self.assertEqual(gaVariant.alternate_bases, alt)
            self._verifyStructuralInfo(gaVariant, pyvcfInfo)
            if not gaVariant.filters_applied:
                self.assertEqual(pyvcfVariant.FILTER, None)
            elif gaVariant.filters_passed:
                self.assertEqual(len(pyvcfVariant.FILTER), 0)
            else:
                self.assertEqual(
                        len(gaVariant.filters_failed),
                        len(pyvcfVariant.FILTER))

            pyvcfCallMap = {}
            for call in pyvcfVariant:
                pyvcfCallMap[call.sample] = call
            _verifyVariantCalls()

    def _verifyGaVariantsSample(self, gaVariants, sampleIds):
        for variant in gaVariants:
            self.assertEqual(len(variant.calls), len(sampleIds))
            for call in variant.calls:
                self.assertIn(call.call_set_name, sampleIds)

    def _verifyVariantsCallSetIds(self, searchsampleIds):
        """
        leaving searchSampleIds will get all samples.
        """
        gaCallSetVariants = []
        for reference_name in self._reference_names:
            end = datamodel.PysamDatamodelMixin.vcfMax
            gaSearchIds = [
                str(datamodel.CallSetCompoundId(
                    self._gaObject.getCompoundId(), sampleId))
                for sampleId in searchsampleIds]
            gaVariants = list(self._gaObject.getVariants(
                reference_name, 0, end, gaSearchIds))
            self._verifyGaVariantsSample(gaVariants, searchsampleIds)
            gaCallSetVariants += gaVariants
            localVariants = filter(
                lambda v: v.CHROM == reference_name, self._variantRecords)
            self._verifyVariantsEqual(gaVariants, localVariants)
        self.assertEqual(len(gaCallSetVariants), len(self._variantRecords))

    def testSearchAllVariants(self):
        self._verifyVariantsCallSetIds(self.vcfSamples[:1])

    def testSearchCallSetIdsSystematic(self):
        for sampleIds in utils.powerset(self.vcfSamples, maxSets=10):
            self._verifyVariantsCallSetIds(list(sampleIds))

    def testVariantsValid(self):
        end = datamodel.PysamDatamodelMixin.vcfMax
        for reference_name in self._reference_names:
            iterator = self._gaObject.getVariants(
                reference_name, 0, end)
            for gaVariant in iterator:
                self.assertValid(protocol.Variant, protocol.toJson(gaVariant))

    def _getPyvcfVariants(
            self, reference_name, startPosition=0, endPosition=2**30):
        """
        variants with in interval [startPosition, endPosition)
        """
        localVariants = []
        for variant in self._variantRecords:
            case1 = (variant.start <= startPosition and
                     variant.end > startPosition)
            case2 = (variant.start > startPosition and
                     variant.start < endPosition)
            if (variant.CHROM == reference_name) and (case1 or case2):
                localVariants.append(variant)
        return localVariants

    def _assertEmptyVariant(self, reference_name, startPosition, endPosition):
        gaVariants = list(self._gaObject.getVariants(
            reference_name, startPosition, endPosition, []))
        self.assertEqual(len(gaVariants), 0)
        pyvcfVariants = self._getPyvcfVariants(
            reference_name, startPosition, endPosition)
        self.assertEqual(len(pyvcfVariants), 0)

    def _assertVariantsEqualInRange(
            self, reference_name, startPosition, endPosition):
        gaVariants = list(self._gaObject.getVariants(
            reference_name, startPosition, endPosition, []))
        pyvcfVariants = self._getPyvcfVariants(
            reference_name, startPosition, endPosition)
        self.assertGreaterEqual(len(gaVariants), 0)
        self._verifyVariantsEqual(gaVariants, pyvcfVariants)

    def testVariantInSegments(self):
        for reference_name in self._reference_names:
            localVariants = self._getPyvcfVariants(reference_name)
            mini = localVariants[0].start
            # NOTE, the end of the last variant may not reflect the END
            maxi = max([v.end for v in localVariants])
            seglen = (maxi-mini) // 3
            seg1 = mini + seglen
            seg2 = seg1 + seglen
            self._assertEmptyVariant(reference_name, -1, mini)
            self._assertEmptyVariant(reference_name, maxi, maxi+100)
            self._assertVariantsEqualInRange(reference_name, mini, seg1)
            self._assertVariantsEqualInRange(reference_name, seg1, seg2)
            self._assertVariantsEqualInRange(reference_name, seg2, maxi)

    def _gaVariantEqualsPyvcfVariant(self, gaVariant, pyvcfVariant):
        if gaVariant.start != pyvcfVariant.start:
            return False
        if gaVariant.end != pyvcfVariant.end:
            return False
        alt = pyvcfVariant.ALT
        if len(alt) == 1 and alt[0] is None:
            alt = []
        if not pyvcfVariant.is_sv and gaVariant.alternate_bases != alt:
            return False
        self._verifyVariantsEqual([gaVariant], [pyvcfVariant])
        return True

    def _pyvcfVariantIsInGaVariants(
            self, pyvcfVariant, intervalStart, intervalEnd):
        isIn = False
        gaVariants = list(self._gaObject.getVariants(
            pyvcfVariant.CHROM, intervalStart, intervalEnd, []))
        for gaVariant in gaVariants:
            if self._gaVariantEqualsPyvcfVariant(gaVariant, pyvcfVariant):
                # FIXME theres an assertion that will cause this to
                # fail if false
                isIn = True
                break
        return isIn

    def testVariantFromToEveryVariant(self):
        for variant in self._variantRecords:
            variantStart = variant.start
            variantEnd = variant.end
            # Cases of search interval does not overlap with variant
            # interval on the left

            self.assertFalse(self._pyvcfVariantIsInGaVariants(
                variant, variantStart-2, variantStart-1))
            self.assertFalse(self._pyvcfVariantIsInGaVariants(
                variant, variantStart-1, variantStart))
            # interval on the right
            self.assertFalse(self._pyvcfVariantIsInGaVariants(
                variant, variantEnd, variantEnd+1))
            self.assertFalse(self._pyvcfVariantIsInGaVariants(
                variant, variantEnd+1, variantEnd+2))
            # case of search interval is within variant

            self.assertTrue(self._pyvcfVariantIsInGaVariants(
                variant, variantStart, variantEnd))
            if (variantEnd - variantStart) != 1:
                self.assertTrue(self._pyvcfVariantIsInGaVariants(
                    variant, variantStart, variantEnd-1))
                self.assertTrue(self._pyvcfVariantIsInGaVariants(
                    variant, variantStart+1, variantEnd))

            # case of search interval contains variant
            self.assertTrue(self._pyvcfVariantIsInGaVariants(
                variant, variantStart-1, variantEnd+1))
            # cases of search interval intersec with variant
            self.assertTrue(self._pyvcfVariantIsInGaVariants(
                variant, variantStart-1, variantStart+1))
            self.assertTrue(self._pyvcfVariantIsInGaVariants(
                variant, variantStart, variantStart+1))
            self.assertTrue(self._pyvcfVariantIsInGaVariants(
                variant, variantEnd-1, variantEnd))
            self.assertTrue(self._pyvcfVariantIsInGaVariants(
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
            # Check that ID is present and has the right 'structure'
            deobfuscatedId = datamodel.CompoundId.deobfuscate(metadata.id)
            splits = datamodel.CompoundId.split(deobfuscatedId)
            self.assertTrue(splits[-1].endswith(metadata.key))

        metadata = keyMap["version"]
        self.assertEqual(metadata.value, self._vcfVersion)

        gtCounter = 0
        hasPass = False
        for prefix, content in [("FORMAT", self._formats),
                                ("INFO", self._infos),
                                ("FILTER", self._filters)]:
            for contentKey in content.keys():
                key = "{0}.{1}".format(prefix, contentKey)
                if key == "FORMAT.GT":
                    gtCounter += 1
                elif key == "FILTER.PASS":
                    hasPass = True
                else:
                    self.assertEqual(
                        keyMap[key].description, content[contentKey].desc)
                    if prefix != "FILTER":
                        self.assertEqual(
                            keyMap[key].type, content[contentKey].type)
                        self.assertEqual(
                                keyMap[key].number,
                                convertPyvcfNumber(content[contentKey].num))
        testMetaLength = (
            1 + len(self._formats) + len(self._infos) + len(self._filters)
            - gtCounter)
        if not hasPass:  # meta-data always has FILTER.PASS
            testMetaLength += 1
        self.assertEqual(len(keyMap), testMetaLength)

    def testGetVariantsCallSets(self):
        variantSet = self._gaObject
        start = 0
        end = datamodel.PysamDatamodelMixin.vcfMax
        call_set_ids = [cs.getId() for cs in variantSet.getCallSets()]
        somecall_set_ids = call_set_ids[0:3]
        for call_set_id in call_set_ids:
            returnedCallSet = variantSet.getCallSet(call_set_id)
            self.assertEqual(call_set_id, returnedCallSet.getId())
        for reference_name in self._reference_names:
            # passing None as the callSetIds argument should be equivalent
            # to passing all of the possible callSetIds as an argument
            noneRecords = list(variantSet.getVariants(
                reference_name, start, end, None))
            allRecords = list(variantSet.getVariants(
                reference_name, start, end, call_set_ids))
            self.assertEqual(len(noneRecords), len(allRecords))
            for noneRecord, allRecord in zip(noneRecords, allRecords):
                for noneCall, allCall in zip(
                        noneRecord.calls, allRecord.calls):
                    self.assertEqual(
                        noneCall.call_set_name, allCall.call_set_name)

            # passing an empty list as the call_set_ids argument should
            # return no callsets for any variant
            emptyRecords = list(variantSet.getVariants(
                reference_name, start, end, []))
            for record in emptyRecords:
                self.assertEqual(len(record.calls), 0)

            # passing some call_set_ids as the call_set_ids argument should
            # return only those calls
            someRecords = list(variantSet.getVariants(
                reference_name, start, end, somecall_set_ids))
            for record in someRecords:
                self.assertEqual(len(record.calls), len(somecall_set_ids))
                for call, someId in zip(record.calls, somecall_set_ids):
                    self.assertEqual(call.call_set_id, someId)

    def testGetVariant(self):
        variantSet = self._gaObject
        for reference_name in self._reference_names:
            refnameVariants = self._getPyvcfVariants(reference_name)
            for variant in refnameVariants:
                # positive test: get the expected variant
                md5 = self._hashVariant(variant)
                compoundId = datamodel.VariantCompoundId(
                    variantSet.getCompoundId(), reference_name,
                    str(variant.start), md5)
                gotVariant = variantSet.getVariant(compoundId)
                self.assertEqual(str(compoundId), gotVariant.id)

                # negative test: change start position to past variant
                wrongStart = variant.end
                compoundId = datamodel.VariantCompoundId(
                    variantSet.getCompoundId(), reference_name,
                    str(wrongStart), md5)
                try:
                    gotVariant = variantSet.getVariant(compoundId)
                    self.assertNotEqual(variant.start, gotVariant.start)
                except exceptions.ObjectNotFoundException:
                    pass

                # negative test: change reference name
                compoundId = datamodel.VariantCompoundId(
                    variantSet.getCompoundId(), "wrong reference name",
                    str(variant.start), md5)
                with self.assertRaises(exceptions.ObjectNotFoundException):
                    variantSet.getVariant(compoundId)

                # negative test: change hash
                compoundId = datamodel.VariantCompoundId(
                    variantSet.getCompoundId(), reference_name,
                    str(variant.start), "wrong hash")
                with self.assertRaises(exceptions.ObjectNotFoundException):
                    variantSet.getVariant(compoundId)

    def _hashVariant(self, record):
        if record.ALT[0] is None:
            alts = tuple()
        else:
            alts = tuple([unicode(sub) for sub in record.ALT])
        hash_str = record.REF + str(alts)
        return hashlib.md5(hash_str).hexdigest()
