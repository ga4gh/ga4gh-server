"""
Builds the test data registry DB.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import glob
import os.path

import ga4gh.common.utils as utils


def run(*args):
    cmd = "python repo_dev.py {}".format(" ".join(args))
    print("running:", cmd)
    utils.runCommand(cmd)


def buildTestData(
        dataDirectory='tests/data', relativePaths=False, force=False):
    prefix = dataDirectory
    repoFile = os.path.join(prefix, "registry.db")
    if os.path.exists(repoFile):
        if force:
            print("deleting repo at '{}'".format(repoFile))
            os.unlink(repoFile)
        else:
            print("'{}' already exists".format(repoFile))
            return
    print("building repo at '{}'".format(repoFile))
    sequenceOntologyName = "so-xp-simple"
    useRelativePath = '-r' if relativePaths else ''
    run("init", "-f", repoFile)

    run("add-peer", repoFile, useRelativePath, "http://example.ga4gh.org")

    pattern = os.path.join(prefix, "referenceSets", "*.fa.gz")
    for dataFile in glob.glob(pattern):
        run("add-referenceset", repoFile, useRelativePath, dataFile,
            "--species ", '\'{"termId": "NCBI:9606", '
                          '"term": "Homo sapiens"}\'')

    pattern = os.path.join(prefix, "ontologies", "*.obo")
    for dataFile in glob.glob(pattern):
        run("add-ontology", repoFile, useRelativePath, dataFile)

    datasetName = "dataset1"
    run("add-dataset", repoFile, datasetName)

    pattern = os.path.join(prefix, "datasets/dataset1/reads", "*.bam")
    for dataFile in glob.glob(pattern):
        run("add-readgroupset", repoFile, datasetName, useRelativePath,
            dataFile)

    pattern = os.path.join(prefix, "datasets/dataset1/variants", "*")
    for j, dataFile in enumerate(glob.glob(pattern)):
        name = "vs_{}".format(j)
        run(
            "add-variantset", repoFile, datasetName, useRelativePath,
            dataFile, "-R NCBI37", "-n ", name, "-aO", sequenceOntologyName)

    pattern = os.path.join(
        prefix, "datasets/dataset1/sequenceAnnotations", "*.db")
    for j, dataFile in enumerate(glob.glob(pattern)):
        run(
            "add-featureset", repoFile, datasetName, useRelativePath,
            dataFile, "-R NCBI37", "-O", sequenceOntologyName,
            "-C ga4gh.datamodel.sequence_annotations.Gff3DbFeatureSet")

    pattern = os.path.join(
        prefix, "datasets/dataset1/continuous", "*.bw")
    for dataFile in glob.glob(pattern):
        run("add-continuousset", repoFile, datasetName, useRelativePath,
            dataFile, "-R NCBI37",
            "-C ga4gh.datamodel.continuous.FileContinuousSet")

    pattern = os.path.join(prefix, "datasets/dataset1/phenotypes", "*")
    for dataFile in glob.glob(pattern):
        # coordinate featureset name and g2p name
        name = dataFile.split("/")[-1]
        run(
            "add-phenotypeassociationset", repoFile,
            datasetName, dataFile, "-n {}".format(name))
        run(
            "add-featureset", repoFile, datasetName, useRelativePath,
            dataFile, "-R NCBI37",  "-O", sequenceOntologyName,
            "-C ga4gh.datamodel.genotype_phenotype_featureset."
            "PhenotypeAssociationFeatureSet")

    pattern = os.path.join(
        prefix, "datasets/dataset1/rnaQuant", "*.db")
    for j, dataFile in enumerate(glob.glob(pattern)):
        name = "rnaseq_{}".format(j)
        run(
            "add-rnaquantificationset", repoFile, datasetName, dataFile,
            "-R NCBI37", "-n ", name)


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--data-directory", default='tests/data', action="store_true",
        help="root directory for test data")
    parser.add_argument(
        "-r", "--relativePaths", default=False, action="store_true",
        help="store relative paths in database")
    parser.add_argument(
        "-f", "--force", default=False, action="store_true",
        help="delete previous database and build a new one")
    args = parser.parse_args()
    return args


@utils.Timed()
def main():
    args = parseArgs()
    buildTestData(args.data_directory, args.relativePaths, args.force)


if __name__ == "__main__":
    main()
