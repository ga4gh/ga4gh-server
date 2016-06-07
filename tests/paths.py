"""
Centralizes hardcoded paths, names, etc. used in tests
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os


testDir = 'tests'
testDataDir = os.path.join(testDir, 'data')
testDataRepo = os.path.join(testDataDir, 'repo.db')

# references
referenceSetName = 'chr17'
faPath = os.path.join(testDataDir, 'referenceSets/Default.fa.gz')
faPath2 = os.path.join(testDataDir, 'referenceSets/example_1.fa.gz')
faPath3 = os.path.join(testDataDir, 'referenceSets/example_2.fa.gz')
ncbi37FaPath = os.path.join(testDataDir, 'referenceSets/NCBI37.fa.gz')

# variants
variantSetName = '1kgPhase1'
vcfDirPath = os.path.join(
    testDataDir, 'datasets/dataset1/variants/1kgPhase1')
vcfDirPath2 = os.path.join(
    testDataDir, 'datasets/dataset1/variants/1kgPhase3')
vcfPath1 = os.path.join(vcfDirPath, 'chr1.vcf.gz')
vcfPath2 = os.path.join(
    vcfDirPath, 'chr2.vcf.gz')
vcfIndexPath1 = os.path.join(
    vcfDirPath, 'chr1.vcf.gz.tbi')
vcfIndexPath2 = os.path.join(
    vcfDirPath, 'chr2.vcf.gz.tbi')
annotatedVcfPath = os.path.join(
    testDataDir, 'datasets/dataset1/variants/1kg.3.annotations')

# Ontologies
ontologyName = "so-xp-simple"
ontologyPath = os.path.join(testDataDir, 'ontologies/so-xp-simple.obo')

# reads
readGroupSetName = 'chr17.1-250'
bamDir = os.path.join(
    testDataDir, 'datasets/dataset1/reads')
bamPath = os.path.join(
    bamDir, 'chr17.1-250.bam')
bamPath2 = os.path.join(
    bamDir, 'wgEncodeUwRepliSeqBg02esG1bAlnRep1_sample.bam')
bamIndexPath = os.path.join(
    bamDir, 'chr17.1-250.bam.bai')
bamIndexPath2 = os.path.join(
    bamDir, 'wgEncodeUwRepliSeqBg02esG1bAlnRep1_sample.bam.bai')

# sequence annotations
featureSetName = 'gencodeV21Set1'
featuresPath = os.path.join(
    testDataDir, 'datasets/dataset1/sequenceAnnotations/gencodeV21Set1.db')
featuresPath2 = os.path.join(
    testDataDir, 'datasets/dataset1/sequenceAnnotations/specialCasesTest.db')

# misc.
landingMessageHtml = os.path.join(testDataDir, "test.html")
