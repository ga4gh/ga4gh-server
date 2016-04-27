"""
A script that takes the compliance dataset (the released version
of which is at https://github.com/ga4gh/compliance/tree/master/test-data)
and turns it into a directory bundle of binary and JSON files suitable
for use by the reference server.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import shutil
import glob
import pysam
import utils
import sys
import generate_gff3_db
import tempfile
import zipfile


class ComplianceDataMunger(object):
    def __init__(self, inputDirectory, outputDirectory):
        """
        Converts human readable dataset from compliance repository,
        and translates it into a reference-server readable filesystem
        with binary files.
        :param inputDirectory: location of
            the human readable compliance dataset
        :param outputDirectory: location of
            the file hierarchy suitable for deploying on the reference server
        """
        self.inputDirectory = inputDirectory
        self.outputDirectory = outputDirectory
        self.tempdir = None

        # If no input directory is specified download from GitHub
        if inputDirectory is None:
            utils.log("Downloading test data...")
            self.tempdir = tempfile.mkdtemp()
            assert(os.path.exists(self.tempdir))
            url = "https://github.com/ga4gh/compliance/archive/master.zip"
            filePath = os.path.join(self.tempdir, 'compliance-master.zip')
            downloader = utils.HttpFileDownloader(url, filePath)
            downloader.download()
            utils.log("Extracting test data...")
            with zipfile.ZipFile(filePath, "r") as z:
                z.extractall(self.tempdir)
            self.inputDirectory = os.path.join(
                self.tempdir, 'compliance-master', 'test-data')

        # get all the reference files (they'll be the ones with .fa extension)
        self.referenceFiles = map(
            os.path.basename, glob.glob(
                os.path.join(self.inputDirectory, "*.fa")))

        self.refsetsDirectory = os.path.join(
            self.outputDirectory, "referenceSets")
        self.hg37Directory = os.path.join(self.refsetsDirectory, "hg37")

        # datasets
        self.datasetsDirectory = os.path.join(self.outputDirectory, "datasets")

        readFiles = map(os.path.basename, glob.glob(
            os.path.join(self.inputDirectory, "*.sam")))
        variantFiles = map(os.path.basename, glob.glob(
            os.path.join(self.inputDirectory, "*.vcf")))
        sequenceAnnotationFiles = map(os.path.basename, glob.glob(
            os.path.join(self.inputDirectory, "*.gff3")))

        self.datasets = [d for d in set(
            [p.split('_')[0] for p in (readFiles +
                                       variantFiles +
                                       sequenceAnnotationFiles)])]
        self.datasetDirs = [os.path.join(
            self.outputDirectory, ds) for ds in self.datasets]

        # Create maps of arrays of files based on dataset.
        self.datasetReads = dict()
        self.datasetVariants = dict()
        self.datasetSequenceAnnotations = dict()

        for ds in self.datasets:
            self.datasetReads[ds] = [r for r in readFiles if r.startswith(ds)]
            self.datasetSequenceAnnotations[ds] = [sa for sa in (
                sequenceAnnotationFiles) if sa.startswith(ds)]

            # Variants themselves are split into groups,
            # based on second part of the _ split:
            self.datasetVariants[ds] = dict()
            # only those variants inside this dataset
            dsvlist = [v for v in variantFiles if v.startswith(ds)]
            # create nested dictionary based on group belonging
            for dsv in dsvlist:
                dsvGroup = dsv.split('_')[1]
                self.datasetVariants[ds][dsvGroup] = \
                    self.datasetVariants[ds].get(dsvGroup, []) + [dsv]

    def run(self):
        if not os.path.exists(self.outputDirectory):
            os.makedirs(self.outputDirectory)

        # Export sequence ontologies
        print("Exporting ontologies...", file=sys.stderr)
        ontologiesDir = os.path.join(self.outputDirectory, "ontologymaps")
        sequenceOntologyDir = os.path.join(ontologiesDir, "sequence_ontology")
        os.makedirs(sequenceOntologyDir)
        shutil.copy(os.path.join(self.inputDirectory, "sequence_ontology.txt"),
                    os.path.join(sequenceOntologyDir, "sequence_ontology.txt"))

        # Clean out, make and re-populate references directory
        # For now, assume a single, statically-named referenceSet
        utils.log("Converting references...")
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
            utils.log("Converting reads...")
            dsReadsdir = os.path.join(dsdir, "reads")
            os.makedirs(dsReadsdir)
            for readFile in self.datasetReads[ds]:
                destFile = os.path.join(
                    dsReadsdir,
                    readFile.split('_')[1].split('.')[0]) + ".bam"
                readSrc = pysam.AlignmentFile(
                    os.path.join(self.inputDirectory, readFile), "r")
                readDest = pysam.AlignmentFile(
                    destFile, "wb", header=readSrc.header)
                destFilePath = readDest.filename

                for readData in readSrc:
                    readDest.write(readData)
                readDest.close()
                readSrc.close()
                pysam.index(destFilePath)

            # Variants
            utils.log("Converting variants...")
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

            # Sequence Annotations
            print("Converting sequence annotations...", file=sys.stderr)
            dsSeqAnndir = os.path.join(dsdir, "sequenceAnnotations")
            os.makedirs(dsSeqAnndir)
            for seqAnnFile in self.datasetSequenceAnnotations[ds]:
                seqAnnDest = os.path.join(
                    dsSeqAnndir,
                    seqAnnFile.split('_')[1].split('.')[0]) + ".db"
                seqAnnSrc = os.path.join(self.inputDirectory, seqAnnFile)
                dbgen = generate_gff3_db.Gff32Db(seqAnnSrc, seqAnnDest)
                dbgen.run()

        print("done converting compliance data.", file=sys.stderr)

    def cleanup(self):
        if self.tempdir is not None:
            shutil.rmtree(self.tempdir)
        utils.log("Done converting compliance data.")
        utils.log("Result in '{}'".format(self.outputDirectory))


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
        help="Path to local directory containing compliance dataset. "
        "If no directory is provided this script will attempt to "
        "download the compliance test-data from github",
        default=None)
    parser.add_argument('--verbose', '-v', action='count', default=0)
    args = parser.parse_args()
    cdm = ComplianceDataMunger(args.inputDirectory, args.outputDirectory)
    try:
        cdm.run()
    finally:
        cdm.cleanup()

if __name__ == "__main__":
    main()
