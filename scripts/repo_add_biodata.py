from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os


import ga4gh.repo_manager as repo_manager


def main():
    rp = repo_manager.RepoManager("ga4gh-example-data")
    for filePath in os.listdir("scripts/biodata/individuals"):
        rp.addIndividual("brca1", os.path.join(
            "scripts/biodata/individuals", filePath), "move")
    for filePath in os.listdir("scripts/biodata/biosamples"):
        rp.addBioSample("brca1", os.path.join(
            "scripts/biodata/biosamples", filePath), "move")

if __name__ == "__main__":
    main()
