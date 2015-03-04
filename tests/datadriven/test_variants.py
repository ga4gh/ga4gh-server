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
        return variants.TabixVariantSet

    def getProtocolClass(self):
        return protocol.GAVariantSet

    def verifyVariantsEqual(self, gaVariants, pyvcfVariants):
        """
        Verifies that the lists of GA4GH variants and pyvcf variants
        are equivalent.
        """
        assert len(gaVariants) == len(pyvcfVariants)
        for gaVariant, pyvcfVariant in zip(gaVariants, pyvcfVariants):
            assert gaVariant.referenceName == pyvcfVariant.CHROM
            assert gaVariant.referenceBases == pyvcfVariant.REF
            # TODO: more tests!

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
        assert len(allVariants) == len(self._variants)
