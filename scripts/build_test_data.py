"""
Builds the test data repo DB.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import glob
import os.path

import utils


def run(*args):
    cmd = "python repo_dev.py {}".format(" ".join(args))
    print("running:", cmd)
    utils.runCommand(cmd)


def main(args):
    prefix = args.data_directory
    repoFile = os.path.join(prefix, "repo.db")
    sequenceOntologyName = "so-xp-simple"
    useRelativePath = '-r' if args.relativePaths else ''
    run("init", "-f", repoFile)

    pattern = os.path.join(prefix, "referenceSets", "*.fa.gz")
    for dataFile in glob.glob(pattern):
        run("add-referenceset", repoFile, useRelativePath, dataFile)

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
            dataFile, "-R NCBI37", "-O", sequenceOntologyName)


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--data-directory", default='tests/data', action="store_true",
        help="root directory for test data")
    parser.add_argument(
        "-r", "--relativePaths", default=False, action="store_true",
        help="store relative paths in database")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parseArgs()
    main(args)
