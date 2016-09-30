"""
Builds the test data registry DB.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import glob
import os
import shutil

import utils


dataDirectoryDefault = 'tests/data'
relativePathsDefault = False
forceDefault = False
targetDirDefault = 'tests/data/repo'


def run(*args):
    cmd = "python repo_dev.py {}".format(" ".join(args))
    print("running:", cmd)
    utils.runCommand(cmd)


def buildTestData(
        dataDirectory=dataDirectoryDefault,
        relativePaths=relativePathsDefault,
        force=forceDefault,
        targetDir=targetDirDefault):
    prefix = dataDirectory
    repoDir = targetDir
    repoFile = os.path.join(repoDir, "registry.db")
    if os.path.exists(repoDir) and os.path.isdir(repoDir):
        if force:
            print("deleting repo at '{}'".format(repoDir))
        else:
            print("'{}' already exists".format(repoDir))
            return
    print("building repo at '{}'".format(repoDir))
    sequenceOntologyName = "so-xp-simple"
    useRelativePath = '-r' if relativePaths else ''
    run("init", "-f", repoDir)

    pattern = os.path.join(prefix, "referenceSets", "*.fa.gz")
    for dataFile in glob.glob(pattern):
        run("add-referenceset", repoDir, useRelativePath, dataFile)

    pattern = os.path.join(prefix, "ontologies", "*.obo")
    for dataFile in glob.glob(pattern):
        run("add-ontology", repoDir, useRelativePath, dataFile)

    datasetName = "dataset1"
    run("add-dataset", repoDir, datasetName)

    pattern = os.path.join(prefix, "datasets/dataset1/reads", "*.bam")
    for dataFile in glob.glob(pattern):
        run("add-readgroupset", repoDir, datasetName, useRelativePath,
            dataFile)

    pattern = os.path.join(prefix, "datasets/dataset1/variants", "*")
    for j, dataFile in enumerate(glob.glob(pattern)):
        name = "vs_{}".format(j)
        run(
            "add-variantset", repoDir, datasetName, useRelativePath,
            dataFile, "-R NCBI37", "-n ", name, "-aO", sequenceOntologyName)

    pattern = os.path.join(
        prefix, "datasets/dataset1/sequenceAnnotations", "*.db")
    for j, dataFile in enumerate(glob.glob(pattern)):
        run(
            "add-featureset", repoDir, datasetName, useRelativePath,
            dataFile, "-R NCBI37", "-O", sequenceOntologyName,
            "-C ga4gh.datamodel.sequence_annotations.Gff3DbFeatureSet")

    pattern = os.path.join(prefix, "datasets/dataset1/phenotypes", "*")
    for dataFile in glob.glob(pattern):
        # coordinate featureset name and g2p name
        name = dataFile.split("/")[-1]
        run(
            "add-phenotypeassociationset", repoDir,
            datasetName, dataFile, "-n {}".format(name))
        run(
            "add-featureset", repoDir, datasetName, useRelativePath,
            dataFile, "-R NCBI37",  "-O", sequenceOntologyName,
            "-C ga4gh.datamodel.genotype_phenotype_featureset."
            "PhenotypeAssociationFeatureSet")

    pattern = os.path.join(
        prefix, "datasets/dataset1/rnaQuant", "*.db")
    for j, dataFile in enumerate(glob.glob(pattern)):
        name = "rnaseq_{}".format(j)
        run(
            "add-rnaquantificationset", repoDir, datasetName, dataFile,
            "-R NCBI37", "-n ", name)


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--data-directory", default=dataDirectoryDefault,
        help="root directory for test data")
    parser.add_argument(
        "-r", "--relativePaths",
        default=relativePathsDefault, action="store_true",
        help="store relative paths in database")
    parser.add_argument(
        "-f", "--force", default=forceDefault, action="store_true",
        help="delete previous database and build a new one")
    parser.add_argument(
        "-t", "--target-dir", default=targetDirDefault,
        help="where the repo should be initialized")
    args = parser.parse_args()
    return args


@utils.Timed()
def main():
    args = parseArgs()
    buildTestData(
        args.data_directory, args.relativePaths, args.force,
        args.target_dir)


if __name__ == "__main__":
    main()
