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
        self._variants = []
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
        for record in vcfReader:
            self._variants.append(record)
            self._referenceNames.add(record.CHROM)

    def getDataModelClass(self):
        return variants.HtslibVariantSet

    def getProtocolClass(self):
        return protocol.GAVariantSet

    def verifyVariantsEqual(self, gaVariants, pyvcfVariants):
        """
        Verifies that the lists of GA4GH variants and pyvcf variants
        are equivalent.
        """
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
            self.assertEqual(gaVariant.alternateBases, alt)
            # TODO check INFO fields

    def testSearchAllVariants(self):
        allVariants = []
        for referenceName in self._referenceNames:
            end = 2**30  # TODO This is arbitrary, and pysam can choke. FIX!
            variants = list(self._gaObject.getVariants(
                referenceName, 0, end, None, None))
            allVariants += variants
            localVariants = [
                v for v in self._variants if v.CHROM == referenceName]
            self.verifyVariantsEqual(variants, localVariants)
        self.assertEqual(len(allVariants), len(self._variants))

    def testVariantSetMetadata(self):
        # TODO get the metadata and test it.
        pass
