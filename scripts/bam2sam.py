"""
Convert a BAM file to a small SAM file
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

import pysam

import utils


def parseArgs():
    parser = argparse.ArgumentParser(
        description="BAM to SAM conversion tool")
    parser.add_argument(
        "inputFile", help="the name of the BAM file to read")
    parser.add_argument(
        "--outputFile", "-o", default='out.sam',
        help="the name of the SAM file to write")
    parser.add_argument(
        "--numLines", "-n", default=10,
        help="the number of lines to write")
    args = parser.parse_args()
    return args


def bam2sam(args):
    bam = pysam.AlignmentFile(args.inputFile, "rb")
    sam = pysam.AlignmentFile(args.outputFile, "wh", header=bam.header)
    for _ in xrange(args.numLines):
        alignedSegment = bam.next()
        sam.write(alignedSegment)
    bam.close()
    sam.close()


@utils.Timed()
def main():
    args = parseArgs()
    bam2sam(args)


if __name__ == '__main__':
    main()
