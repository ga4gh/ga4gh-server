"""
Data-driven tests for variants.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob

import vcf

import ga4gh.server.datamodel as datamodel
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.datamodel.variants as variants
import ga4gh.server.datamodel.references as references
import ga4gh.server.datamodel.ontologies as ontologies
import tests.datadriven as datadriven
import tests.paths as paths

import ga4gh.schemas.protocol as protocol


def testVariantAnnotationSets():
    testDataDir = "tests/data/datasets/dataset1/variants"
    for test in datadriven.makeTests(testDataDir, VariantAnnotationSetTest):
        yield test


class VariantAnnotationSetTest(datadriven.DataDrivenTest):
    """
    Data driven test class for variant annotation sets. Builds an
    alternative model of a variant annotation set using pyvcf, and
    verifies that it is consistent with the model built by the
    variants.VariantAnnotationSet object.
    """
    def __init__(self, variantAnnotationSetId, baseDir):
        super(VariantAnnotationSetTest, self).__init__(
            variantAnnotationSetId, baseDir)
        self._variantRecords = []
        self._referenceNames = set()
        # Only read in VCF files that are annotated.
        for vcfFile in glob.glob(os.path.join(self._dataPath, "*.vcf.gz")):
            if self._isAnnotated():
                self._readVcf(vcfFile)
        self._isCsq = self._hasConsequenceField()

    def _isAnnotated(self):
        """
        Determines whether the variant set under test is annotated.

        :return: Boolean
        """
        pyvcfreader = vcf.Reader(
            filename=glob.glob(
                os.path.join(self._dataPath, "*.vcf.gz"))[0])
        items = [x for x in pyvcfreader.infos]
        return ('ANN' in items) or ('CSQ' in items)

    def _hasConsequenceField(self):
        pyvcfreader = vcf.Reader(
            filename=glob.glob(
                os.path.join(self._dataPath, "*.vcf.gz"))[0])
        items = [x for x in pyvcfreader.infos]
        return 'CSQ' in items

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
        if "VEP" in metadata:
            self._isVEP = True
        else:
            self._isVEP = False
        for record in vcfReader:
            self._referenceNames.add(record.CHROM)
            # When an END info tag is present it takes precedence
            if "END" in record.INFO:
                record.end = record.INFO["END"]
            self._variantRecords.append(record)

    def splitAnnField(self, annfield):
        """
        Splits the `ANN=` record of a VCF into a dictionary of
        key, value. This is brittle but covers the cases of the

        :param annfield:
        :return: dictionary of k, v pairs along from the ann string
        """
        # Have yet to see if other files follow this pattern
        if self._isVEP and not self._isCsq:
            fields = ["alt", "effects", "impact", "symbol", "geneName",
                      "featureType", "featureId", "trBiotype", "exon",
                      "intron", "hgvsC", "hgvsP",
                      "cdnaPos", "cdsPos", "protPos",
                      "aminos", "codons", "existingVar", "distance",
                      "strand", "symbolSource", "hgncId", "hgvsOffset"]
        else:
            fields = ["alt", "effects", "impact", "geneName",
                      "geneId", "featureType",
                      "featureId", "trBiotype", "rank", "hgvsC",
                      "hgvsP", "cdnaPos", "cdsPos", "protPos",
                      "distance", "errsWarns"]
        values = annfield.split('|')
        values[1] = values[1].split('&')  # split effects
        return dict(zip(fields, values))

    def getDataModelInstance(self, localId, dataPath):
        dataset = datasets.Dataset("ds")
        variantSet = variants.HtslibVariantSet(dataset, localId)
        variantSet.populateFromDirectory(dataPath)
        referenceSet = references.AbstractReferenceSet("rs")
        variantSet.setReferenceSet(referenceSet)
        if variantSet.isAnnotated():
            sequenceOntology = ontologies.Ontology(paths.ontologyName)
            sequenceOntology.populateFromFile(paths.ontologyPath)
            annotationSet = variantSet.getVariantAnnotationSets()[0]
            annotationSet.setOntology(sequenceOntology)
            return annotationSet
        else:
            return variantSet

    def getProtocolClass(self):
        if self._isAnnotated():
            return protocol.VariantAnnotationSet
        else:
            return protocol.VariantSet

    def _compareTwoListFloats(self, a, b):
        for ai, bi in zip(a, b):
            self.assertAlmostEqual(float(ai), float(bi), 5)

    def _verifyVariantAnnotationsEqual(
            self, gaVariantAnnotations, pyvcfVariants):
        """
        Verifies that the lists of GA4GH variants and pyvcf variants
        are equivalent.
        """
        self.assertEqual(len(gaVariantAnnotations), len(pyvcfVariants))
        for pair, pyvcfVariant in zip(gaVariantAnnotations, pyvcfVariants):
            gaVariant, gaVariantAnnotation = pair
            # TODO parse and compare with analysis.info
            # pyvcfInfo = pyvcfVariant.INFO
            # pyvcf uses 1-based indexing.
            self.assertEqual(gaVariant.id, gaVariantAnnotation.variant_id)
            self.assertEqual(gaVariant.start, pyvcfVariant.POS - 1)
            self.assertEqual(gaVariant.end, pyvcfVariant.end)
            # Annotated VCFs contain an ANN field in the record info
            if 'ANN' in pyvcfVariant.INFO:
                pyvcfAnn = pyvcfVariant.INFO['ANN']
                i = 0
                for pyvcfEffect, gaEff in \
                        zip(pyvcfAnn, gaVariantAnnotation.transcript_effects):
                    effectDict = self.splitAnnField(pyvcfEffect)
                    self.assertEqual(
                        map(lambda e: e.term, gaEff.effects),
                        effectDict['effects'])
                    self.assertEqual(
                        gaEff.attributes.attr['impact'].values[0].string_value,
                        effectDict['impact'])
                    self.assertEqual(
                        gaEff.alternate_bases, effectDict['alt'])
                    self.assertEqual(
                        gaEff.feature_id, effectDict['featureId'])
                    self.assertEqual(
                        gaEff.hgvs_annotation.transcript,
                        effectDict['hgvsC'])
                    self.assertEqual(
                        gaEff.hgvs_annotation.protein, effectDict['hgvsP'])
                    if 'HGVS.g' in pyvcfVariant.INFO:
                        # Not all VCF have this field
                        index = i % len(pyvcfVariant.INFO['HGVS.g'])
                        self.assertEqual(
                            gaEff.hgvs_annotation.genomic,
                            pyvcfVariant.INFO['HGVS.g'][index])
                    i += 1
            elif 'CSQ' in pyvcfVariant.INFO:
                pyvcfAnn = pyvcfVariant.INFO['CSQ']
                transcriptEffects = [
                    self._splitCsqEffects(ann) for ann in pyvcfAnn]
                for treff, gaEff in zip(
                        transcriptEffects,
                        gaVariantAnnotation.transcript_effects):
                    self.assertEqual(
                        map(lambda e: e.term, gaEff.effects),
                        treff['effects'])
                    self.assertEqual(gaEff.alternate_bases, treff['alt'])
                    self.assertEqual(gaEff.feature_id, treff['featureId'])
            self.assertGreater(len(gaVariantAnnotation.transcript_effects), 0)

    def _splitCsqEffects(self, annStr):
        (alt, gene, featureId, featureType, effects, cdnaPos,
            cdsPos, protPos, aminos, codons, existingVar,
            distance, strand, sift, polyPhen, motifName,
            motifPos, highInfPos,
            motifScoreChange) = annStr.split('|')
        effect = {}
        effect['featureId'] = featureId
        effect['alt'] = alt
        effect['effects'] = effects.split("&")
        return effect

    def testVariantsValid(self):
        end = datamodel.PysamDatamodelMixin.vcfMax
        for referenceName in self._referenceNames:
            iterator = self._gaObject.getVariantAnnotations(
                referenceName, 0, end)
            for gaVariant, gaVariantAnnotation in iterator:
                self.assertValid(protocol.VariantAnnotation,
                                 protocol.toJson(gaVariantAnnotation))

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
        gaVariantAnnotations = list(self._gaObject.getVariantAnnotations(
            referenceName, startPosition, endPosition))
        self.assertEqual(len(gaVariantAnnotations), 0)
        pyvcfVariants = self._getPyvcfVariants(
            referenceName, startPosition, endPosition)
        self.assertEqual(len(pyvcfVariants), 0)

    def _assertVariantsEqualInRange(
            self, referenceName, startPosition, endPosition):
        gaVariantAnnotations = list(
            self._gaObject.getVariantAnnotations(
                referenceName, startPosition, endPosition))
        pyvcfVariants = self._getPyvcfVariants(
            referenceName, startPosition, endPosition)
        self.assertGreaterEqual(len(gaVariantAnnotations), 0)
        self._verifyVariantAnnotationsEqual(
            gaVariantAnnotations, pyvcfVariants)

    def testVariantAnnotationInSegments(self):
        for referenceName in self._referenceNames:
            localVariants = self._getPyvcfVariants(referenceName)
            mini = localVariants[0].start
            # NOTE, the end of the last variant may not reflect the END
            maxi = max([v.end for v in localVariants])
            seglen = (maxi-mini) // 3
            seg1 = mini + seglen
            seg2 = seg1 + seglen
            self._assertEmptyVariant(referenceName, -1, mini)
            self._assertEmptyVariant(referenceName, maxi, maxi+100)
            self._assertVariantsEqualInRange(referenceName, mini, seg1)
            self._assertVariantsEqualInRange(referenceName, seg1, seg2)
            self._assertVariantsEqualInRange(referenceName, seg2, maxi)

    def _gaVariantAnnotationEqualsPyvcfVariantAnnotation(
            self, gaVariant, gaAnnotation, pyvcfVariant):
        if gaVariant.start != pyvcfVariant.start:
            return False
        if gaVariant.end != pyvcfVariant.end:
            return False
        return True

    def _pyvcfVariantAnnotationIsInGaVariantAnnotations(
            self, pyvcfVariant, intervalStart, intervalEnd):
        isIn = False
        pairs = list(self._gaObject.getVariantAnnotations(
            pyvcfVariant.CHROM, intervalStart, intervalEnd))
        for gaVariant, gaAnnotation in pairs:
            if self._gaVariantAnnotationEqualsPyvcfVariantAnnotation(
                    gaVariant, gaAnnotation, pyvcfVariant):
                isIn = True
                break
        return isIn

    def testVariantAnnotationFromToEveryVariantAnnotation(self):
        for variant in self._variantRecords:
            variantStart = variant.start
            variantEnd = variant.end
            # Cases of search interval does not overlap with variant
            # interval on the left

            self.assertFalse(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantStart-2, variantStart-1))
            self.assertFalse(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantStart-1, variantStart))
            # interval on the right
            self.assertFalse(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantEnd, variantEnd+1))
            self.assertFalse(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantEnd+1, variantEnd+2))
            # case of search interval is within variant

            self.assertTrue(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantStart, variantEnd))
            if (variantEnd - variantStart) != 1:
                self.assertTrue(
                    self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                        variant, variantStart, variantEnd-1))
                self.assertTrue(
                    self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                        variant, variantStart+1, variantEnd))

            # case of search interval contains variant
            self.assertTrue(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantStart-1, variantEnd+1))
            # cases of search interval intersec with variant
            self.assertTrue(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantStart-1, variantStart+1))
            self.assertTrue(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantStart, variantStart+1))
            self.assertTrue(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantEnd-1, variantEnd))
            self.assertTrue(
                self._pyvcfVariantAnnotationIsInGaVariantAnnotations(
                    variant, variantEnd-1, variantEnd+1))
