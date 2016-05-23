"""
Builds the test data repo DB.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import glob
import os.path

import utils


def run(*args):
    cmd = "python repo_dev.py {}".format(" ".join(args))
    print("running:", cmd)
    utils.runCommand(cmd)


def main():
    prefix = "tests/data"
    repoFile = os.path.join(prefix, "repo.db")
    sequenceOntologyName = "so-xp-simple"
    run("init", "-f", repoFile)

    pattern = os.path.join(prefix, "referenceSets", "*.fa.gz")
    for dataFile in glob.glob(pattern):
        run("add-referenceset", repoFile, dataFile)

    pattern = os.path.join(prefix, "ontologies", "*.obo")
    for dataFile in glob.glob(pattern):
        run("add-ontology", repoFile, dataFile)

    datasetName = "dataset1"
    run("add-dataset", repoFile, datasetName)

    pattern = os.path.join(prefix, "datasets/dataset1/reads", "*.bam")
    for dataFile in glob.glob(pattern):
        run("add-readgroupset", repoFile, datasetName, dataFile)

    pattern = os.path.join(prefix, "datasets/dataset1/variants", "*")
    for j, dataFile in enumerate(glob.glob(pattern)):
        name = "vs_{}".format(j)
        run(
            "add-variantset", repoFile, datasetName, dataFile, "-R NCBI37",
            "-n ", name, "-aO", sequenceOntologyName)

    pattern = os.path.join(
        prefix, "datasets/dataset1/sequenceAnnotations", "*.db")
    for j, dataFile in enumerate(glob.glob(pattern)):
        run(
            "add-featureset", repoFile, datasetName, dataFile, "-R NCBI37",
            "-O", sequenceOntologyName)


if __name__ == "__main__":
    main()
