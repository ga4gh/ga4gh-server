"""
Centralizes hardcoded paths, names, etc. used in tests
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os


packageName = 'ga4gh'


def getProjectRootFilePath():
    # assumes we're in a directory one level below the project root
    return os.path.dirname(os.path.dirname(__file__))


def getGa4ghFilePath():
    return os.path.join(getProjectRootFilePath(), packageName)


testDir = 'tests'
testDataDir = os.path.join(testDir, 'data')
testDataRepo = os.path.join(testDataDir, 'registry.db')

# datasets
datasetName = "dataset1"
datasetsDir = os.path.join(testDataDir, "datasets")
datasetDir = os.path.join(datasetsDir, "dataset1")

# references
referenceSetName = 'chr17'
referenceSetsDir = os.path.join(testDataDir, 'referenceSets')
faPath = os.path.join(referenceSetsDir, 'Default.fa.gz')
faPath2 = os.path.join(referenceSetsDir, 'example_1.fa.gz')
faPath3 = os.path.join(referenceSetsDir, 'example_2.fa.gz')
ncbi37FaPath = os.path.join(referenceSetsDir, 'NCBI37.fa.gz')

# variants
variantSetName = '1kgPhase1'
variantsDir = os.path.join(datasetDir, 'variants')
vcfDirPath = os.path.join(variantsDir, '1kgPhase1')
vcfDirPath2 = os.path.join(variantsDir, '1kgPhase3')
vcfPath1 = os.path.join(vcfDirPath, 'chr1.vcf.gz')
vcfPath2 = os.path.join(vcfDirPath, 'chr2.vcf.gz')
vcfIndexPath1 = os.path.join(vcfDirPath, 'chr1.vcf.gz.tbi')
vcfIndexPath2 = os.path.join(vcfDirPath, 'chr2.vcf.gz.tbi')
annotatedVcfPath = os.path.join(variantsDir, '1kg.3.annotations')

# Ontologies
ontologyName = "so-xp-simple"
ontologiesDir = os.path.join(testDataDir, 'ontologies')
ontologyPath = os.path.join(ontologiesDir, 'so-xp-simple.obo')

# reads
readGroupSetName = 'chr17.1-250'
bamDir = os.path.join(datasetDir, 'reads')
bamPath = os.path.join(bamDir, 'chr17.1-250.bam')
bamPath2 = os.path.join(
    bamDir, 'wgEncodeUwRepliSeqBg02esG1bAlnRep1_sample.bam')
bamIndexPath = os.path.join(bamDir, 'chr17.1-250.bam.bai')
bamIndexPath2 = os.path.join(
    bamDir, 'wgEncodeUwRepliSeqBg02esG1bAlnRep1_sample.bam.bai')

# sequence annotations
featureSetName = 'gencodeV21Set1'
featuresDir = os.path.join(datasetDir, 'sequenceAnnotations')
featuresPath = os.path.join(featuresDir, 'gencodeV21Set1.db')
featuresPath2 = os.path.join(featuresDir, 'specialCasesTest.db')

# continuous
continuousSetName = 'bigwig_1'
continuousDir = os.path.join(datasetDir, 'continuous')
continuousPath = os.path.join(continuousDir, 'bigwig_1.bw')

# g2p
phenotypesDir = os.path.join(datasetDir, 'phenotypes')
phenotypeAssociationSetPath = os.path.join(phenotypesDir, 'cgd')

# rna
rnaQuantDir = os.path.join(datasetDir, "rnaQuant")
rnaQuantificationSetName = "ENCFF305LZB"
rnaQuantificationSetDbPath = os.path.join(rnaQuantDir, "ENCFF305LZB.db")

# misc.
landingMessageHtml = os.path.join(testDataDir, "test.html")

peerUrl = "http://1kgenomes.ga4gh.org"
