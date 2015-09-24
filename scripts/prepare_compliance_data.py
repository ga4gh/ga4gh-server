"""
A script that takes the compliance dataset (the released version
of which is at https://github.com/ga4gh/compliance/tree/master/test-data)
and turns it into a directory bundle of binary and JSON files suitable
for use by the reference server.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import shutil
import argparse
import os
import sys
import glob
import pysam
import utils


class ComplianceDataMunger(object):
    def __init__(self, args):
        self.inputDirectory = args.inputDirectory
        self.outputDirectory = args.outputDirectory

        # get all the reference files (they'll be the ones with .fa extension)
        self.referenceFiles = map(
            os.path.basename, glob.glob(
                os.path.join(self.inputDirectory, "*.fa")))

        self.refsetsDirectory = os.path.join(
            self.outputDirectory, "referenceSets")
        self.hg37Directory = os.path.join(self.refsetsDirectory, "hg37")

        # datasets
        self.datasetsDirectory = os.path.join(self.outputDirectory, "datasets")
        self.readFiles = map(
            os.path.basename, glob.glob(
                os.path.join(self.inputDirectory, "*.sam")))
        self.variantFiles = map(
            os.path.basename, glob.glob(
                os.path.join(self.inputDirectory, "*.vcf")))

        self.datasets = [d for d in set([p.split('_')[0] for p in
                                         self.readFiles + self.variantFiles])]
        self.datasetReads = dict()
        self.datasetVariants = dict()
        for ds in self.datasets:
            self.datasetReads[ds] = [r for r in
                                     self.readFiles if r.startswith(ds)]

        # Variants themselves are split into groups,
        # based on second part of the _ split:
        for ds in self.datasets:
            self.datasetVariants[ds] = dict()
            # only those variants inside this dataset
            dsvlist = [v for v in self.variantFiles if v.startswith(ds)]
            # create nested dictionary based on group belonging
            for dsv in dsvlist:
                dsvGroup = dsv.split('_')[1]
                self.datasetVariants[ds][dsvGroup] = \
                    self.datasetVariants[ds].get(dsvGroup, []) + [dsv]

        self.datasetDirs = [os.path.join(self.outputDirectory, ds)
                            for ds in self.datasets]

    def run(self):
        if not os.path.exists(self.outputDirectory):
            os.makedirs(self.outputDirectory)

        # Clean out, make and re-populate references directory
        # For now, assume a single, statically-named referenceSet
        print("Converting references...", file=sys.stderr)
        shutil.rmtree(self.refsetsDirectory, ignore_errors=True)
        os.makedirs(self.refsetsDirectory)
        shutil.copy(
            os.path.join(self.inputDirectory, "referenceset_hg37.json"),
            os.path.join(self.refsetsDirectory, "hg37.json"))

        os.makedirs(self.hg37Directory)
        for refFile in self.referenceFiles:
            refBase = os.path.splitext(refFile)[0]
            destFastaFilename = os.path.join(
                self.hg37Directory, refBase) + ".fa"
            shutil.copy(os.path.join(self.inputDirectory, refBase) + ".fa",
                        destFastaFilename)
            pysam.tabix_compress(destFastaFilename, destFastaFilename + ".gz")
            refFasta = pysam.FastaFile(destFastaFilename + ".gz")
            refFasta.close()
            os.remove(destFastaFilename)
            shutil.copy(
                os.path.join(self.inputDirectory, refBase) + ".json",
                os.path.join(self.hg37Directory, refBase) + ".json")

        # Clean out, make and repopulate dataset directories
        shutil.rmtree(self.datasetsDirectory, ignore_errors=True)
        os.makedirs(self.datasetsDirectory)

        for ds in self.datasets:
            dsdir = os.path.join(self.datasetsDirectory, ds)
            os.makedirs(dsdir)

            # Reads
            print("Converting reads...", file=sys.stderr)
            dsReadsdir = os.path.join(dsdir, "reads")
            os.makedirs(dsReadsdir)
            for readFile in self.datasetReads[ds]:
                destFile = os.path.join(
                    dsReadsdir,
                    readFile.split('_')[1].split('.')[0]) + ".bam"
                readSrc = pysam.AlignmentFile(
                    os.path.join(self.inputDirectory, readFile), "r")
                readDest = pysam.AlignmentFile(destFile, "wb",
                                               header=readSrc.header)
                destFilePath = readDest.filename

                for readData in readSrc:
                    readDest.write(readData)
                readDest.close()
                readSrc.close()
                pysam.index(destFilePath)

            # Variants
            print("Converting variants...", file=sys.stderr)
            dsVariantsdir = os.path.join(dsdir, "variants")
            os.makedirs(dsVariantsdir)
            for vgroup in self.datasetVariants[ds].keys():
                vgroupdir = os.path.join(dsVariantsdir, vgroup)
                os.makedirs(vgroupdir)
                for variantFile in self.datasetVariants[ds][vgroup]:
                    destFile = os.path.join(
                        vgroupdir, variantFile.split('_')[2])
                    shutil.copy(
                        os.path.join(
                            self.inputDirectory, variantFile), destFile)
                    # Pysam's tabix_index automatically compresses the file
                    # in place, creates a tabix index.
                    pysam.tabix_index(destFile, preset="vcf")

        print("done converting compliance data.", file=sys.stderr)


@utils.Timed()
def main():
    parser = argparse.ArgumentParser(
        description="Script to generate data bundle from a locally stored "
        "(and possibly locally edited) version of the compliance dataset.")
    parser.add_argument(
        "--outputDirectory", "-o", default="ga4gh-compliance-data",
        help="The directory to output the server-ready data bundle to.")
    parser.add_argument(
        "--inputDirectory", "-i",
        help="Path to local directory containing compliance dataset.",
        default='.')
    parser.add_argument('--verbose', '-v', action='count', default=0)
    args = parser.parse_args()
    cdm = ComplianceDataMunger(args)
    cdm.run()

if __name__ == "__main__":
    main()
