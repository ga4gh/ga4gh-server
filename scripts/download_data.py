"""
Constructs a data source for the ga4gh server by downloading data from
authoritative remote servers.

Note: checkpoint algorithm: htslib often crashes, and we use it a lot in
this script.  To prevent against having to run the entire script again,
which could take a long time, we checkpoint progress made by leaving a file
named X.checkpoint in certain directories, where X could be a chromosome,
etc. that has been successfully processed to completion.  If the script is
started again after a failed run, it sees the checkpoint files, skips
processing of that object, and resumes at the point of the crash.  At the
end of the script, we delete all of the checkpoint files.  (We can not
simply use the presence of the expected processed file as evidence of
successfully processing that object since that file, while extant, may be
corrupted / half-finished.  To the extent that we can, it is also nice
to be able to not assume anything about how the object is processed, only
that some processing happens and it either succeeds or fails; this way we
keep a separation of concerns between checkpointing and processing.)
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import gzip
import json
import hashlib
import os
import requests
import tempfile
import urllib2

import pysam

import utils


def getCheckpointFileName(objName):
    """
    Return the name of a checkpoint file for a given object
    """
    return "{}.checkpoint".format(objName)


def isCheckpointFileExtant(objName):
    """
    Return true if a checkpoint file for an object exists, false otherwise
    """
    checkpointFileName = getCheckpointFileName(objName)
    return os.path.exists(checkpointFileName)


def createCheckpointFile(objName):
    """
    Create a checkpoint file for the object
    """
    open(getCheckpointFileName(objName), 'a').close()


def logCheckpoint(objName):
    """
    Indicate the object is being skipped because it has been checkpointed
    """
    utils.log("Checkpoint found for {} ... skipping ...".format(objName))


def mkdirAndChdirList(dirNameList):
    """
    For each entry in dirNameList, make a directory if it does not exist and
    cd into it
    """
    for directory in dirNameList:
        mkdirAndChdir(directory)


def mkdirAndChdir(dirName):
    """
    Make a directory if it does not exist and cd into it
    """
    if not os.path.exists(dirName):
        os.mkdir(dirName)
    os.chdir(dirName)


def getReferenceChecksum(fastaFile):
    """
    Returns the md5 checksum for the reference sequence in the specified
    FASTA file. This is the MD5 of the upper case sequence letters.
    """
    inputFile = pysam.FastaFile(fastaFile)
    bases = inputFile.fetch(inputFile.references[0])
    inputFile.close()
    return hashlib.md5(bases.upper()).hexdigest()


def escapeDir(levels):
    """
    Walk back through the number of parent directories equal to levels
    """
    for _ in range(levels):
        os.chdir('..')


def dumpDictToFileAsJson(data, filename):
    """
    Writes a dictionary to a file as a json dump
    """
    with open(filename, 'w') as dumpFile:
        json.dump(data, dumpFile, indent=4)


class ChromMinMax(object):
    """
    A container class for storing the min and max position seen
    for every chromosome
    """
    defaultMinPos = 2**30
    defaultMaxPos = 0

    class MinMax(object):
        def __init__(self):
            self.minPos = ChromMinMax.defaultMinPos
            self.maxPos = ChromMinMax.defaultMaxPos

    def __init__(self):
        self.chromMap = {}

    def addPos(self, chrom, position):
        if chrom not in self.chromMap:
            self.chromMap[chrom] = self.MinMax()
        minMax = self.chromMap[chrom]
        if position > minMax.maxPos:
            minMax.maxPos = position
        if position < minMax.minPos:
            minMax.minPos = position

    def getMinPos(self, chrom):
        minMax = self.chromMap[chrom]
        return minMax.minPos

    def getMaxPos(self, chrom):
        minMax = self.chromMap[chrom]
        return minMax.maxPos


def _fetchSequence(ac, startIndex=None, endIndex=None):
    """Fetch sequences from NCBI using the eutils interface.

    An interbase interval may be optionally provided with startIndex and
    endIndex. NCBI eutils will return just the requested subsequence, which
    might greatly reduce payload sizes (especially with chromosome-scale
    sequences). When wrapped is True, return list of sequence lines rather
    than concatenated sequence.

    >>> len(_fetchSequence('NP_056374.2'))
    1596

    Pass the desired interval rather than using Python's [] slice
    operator.

    >>> _fetchSequence('NP_056374.2',0,10)
    'MESRETLSSS'

    >>> _fetchSequence('NP_056374.2')[0:10]
    'MESRETLSSS'

    """
    urlFmt = ("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
              "db=nucleotide&id={ac}&rettype=fasta")
    if startIndex is None or endIndex is None:
        url = urlFmt.format(ac=ac)
    else:
        urlFmt += "&seq_start={start}&seq_stop={stop}"
        url = urlFmt.format(ac=ac, start=startIndex + 1, stop=endIndex)
    resp = requests.get(url)
    resp.raise_for_status()
    seqlines = resp.content.splitlines()[1:]
    print("{ac}[{s},{e}) => {n} lines ({u})".format(
        ac=ac, s=startIndex, e=endIndex, n=len(seqlines), u=url))
    # return response as list of lines, already line wrapped
    return seqlines


class AbstractFileDownloader(object):
    """
    Base class for individual site genome file downloaders
    """
    def __init__(self, args):
        self.excludeReferenceMin = args.exclude_reference_min
        self.maxVariants = args.num_variants
        self.maxReads = args.num_reads
        self.samples = args.samples.split(',')
        self.numChromosomes = args.num_chromosomes
        self.chromosomes = [str(j + 1) for j in range(self.numChromosomes)]
        self.dirName = args.dir_name
        self.datasetName = '1kg-p3-subset'
        self.variantSetName = 'mvncall'
        self.referenceSetName = 'GRCh37-subset'
        self.chromMinMax = ChromMinMax()
        self.accessions = {
            '1': 'CM000663.1',
            '2': 'CM000664.1',
            '3': 'CM000665.1',
            '4': 'CM000666.1',
            '5': 'CM000667.1',
            '6': 'CM000668.1',
            '7': 'CM000669.1',
            '8': 'CM000670.1',
            '9': 'CM000671.1',
            '10': 'CM000672.1',
            '11': 'CM000673.1',
            '12': 'CM000674.1',
            '13': 'CM000675.1',
            '14': 'CM000676.1',
            '15': 'CM000677.1',
            '16': 'CM000678.1',
            '17': 'CM000679.1',
            '18': 'CM000680.1',
            '19': 'CM000681.1',
            '20': 'CM000682.1',
            '21': 'CM000683.1',
            '22': 'CM000684.1',
        }
        self.studyMap = {
            'HG00096': 'GBR',
            'HG00533': 'CHS',
            'HG00534': 'CHS',
        }

    def getVcfBaseUrl(self):
        return os.path.join(self.getBaseUrl(), 'ftp/release/20130502/')

    def getBamBaseUrl(self):
        return os.path.join(self.getBaseUrl(), 'ftp/phase3/data/')

    def _updatePositions(self, fileName):
        localVariantFile = pysam.VariantFile(fileName)
        localIterator = localVariantFile.fetch()
        for record in localIterator:
            self.chromMinMax.addPos(record.chrom, record.start)
        localIterator = None
        localVariantFile.close()
        utils.log('chrom: {}, maxPos: {}, minPos: {}'.format(
            record.chrom, self.chromMinMax.getMaxPos(record.chrom),
            self.chromMinMax.getMinPos(record.chrom)))

    def _writeVcfTempFile(self, localTempFileName, data):
        with tempfile.NamedTemporaryFile() as binaryFile:
            binaryFile.write(data)
            binaryFile.flush()
            gzipFile = gzip.open(binaryFile.name, "r")
            outputFile = open(localTempFileName, "w")
            lineCount = 0
            for line in gzipFile:
                outputFile.write(line)
                if not line.startswith("#"):
                    lineCount += 1
                if lineCount >= self.maxVariants:
                    break
            assert lineCount == self.maxVariants
            outputFile.close()
            gzipFile.close()

    def _downloadVcf(self, chromosome):
        sourceFileName = (
            "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a"
            ".20130502.genotypes.vcf.gz").format(chromosome)
        url = os.path.join(self.getVcfBaseUrl(), sourceFileName)
        utils.log("Downloading '{}'".format(url))
        response = urllib2.urlopen(url)
        megabyte = 1024 * 1024
        data = response.read(megabyte)
        localFileName = "{}.vcf".format(chromosome)
        localCompressedFileName = "{}.gz".format(localFileName)
        localTempFileName = localFileName + '.unsampled'
        utils.log("Writing '{}'".format(localTempFileName))
        self._writeVcfTempFile(localTempFileName, data)
        utils.log("Sampling '{}'".format(localTempFileName))
        utils.runCommand(
            'bcftools view --force-samples -s {} {} -o {}'.format(
                args.samples, localTempFileName, localFileName))
        os.remove(localTempFileName)
        utils.log("Compressing '{}'".format(localFileName))
        utils.runCommand('bgzip -f {}'.format(localFileName))
        utils.log("Indexing '{}'".format(localCompressedFileName))
        utils.runCommand('tabix {}'.format(localCompressedFileName))
        self._updatePositions(localCompressedFileName)

    def downloadVcfs(self):
        dirList = [
            self.dirName, "datasets", self.datasetName, 'variants',
            self.variantSetName]
        mkdirAndChdirList(dirList)
        for chromosome in self.chromosomes:
            if isCheckpointFileExtant(chromosome):
                logCheckpoint(chromosome)
                localCompressedFileName = '{}.vcf.gz'.format(chromosome)
                self._updatePositions(localCompressedFileName)
            else:
                self._downloadVcf(chromosome)
                createCheckpointFile(chromosome)
        escapeDir(5)

    def createBamHeader(self, baseHeader):
        """
        Creates a new bam header based on the specified header from the
        parent BAM file.
        """
        header = dict(baseHeader)
        newSequences = []
        for index, referenceInfo in enumerate(header['SQ']):
            if index < self.numChromosomes:
                referenceName = referenceInfo['SN']
                # The sequence dictionary in the BAM file has to match up
                # with the sequence ids in the data, so we must be sure
                # that these still match up.
                assert referenceName == self.chromosomes[index]
                newReferenceInfo = {
                    'AS': self.referenceSetName,
                    'SN': referenceName,
                    'LN': 0,  # FIXME
                    'UR': 'http://example.com',
                    'M5': 'dbb6e8ece0b5de29da56601613007c2a',  # FIXME
                    'SP': 'Human'
                }
                newSequences.append(newReferenceInfo)
        header['SQ'] = newSequences
        return header

    def _downloadBam(self, sample):
        samplePath = '{}/alignment/'.format(sample)
        study = self.studyMap[sample]
        sourceFileName = (
            '{}.mapped.ILLUMINA.bwa.{}.'
            'low_coverage.20120522.bam'.format(sample, study))
        destFileName = "{}.bam".format(sample)
        baseUrl = self.getBamBaseUrl()
        sampleUrl = os.path.join(baseUrl, samplePath, sourceFileName)
        utils.log("Downloading index for '{}'".format(sampleUrl))
        remoteFile = pysam.AlignmentFile(sampleUrl)
        header = self.createBamHeader(remoteFile.header)
        utils.log("Writing '{}'".format(destFileName))
        localFile = pysam.AlignmentFile(
            destFileName, 'wb', header=header)
        for chromosome in self.chromosomes:
            utils.log("chromosome {}".format(chromosome))
            iterator = remoteFile.fetch(
                chromosome.encode('utf-8'),
                start=self.chromMinMax.getMinPos(chromosome),
                end=self.chromMinMax.getMaxPos(chromosome))
            for index, record in enumerate(iterator):
                # We only write records where we have the references
                # for the next mate. TODO we should take the positions
                # of these reads into account later when calculating
                # our reference bounds.
                if record.next_reference_id < self.numChromosomes:
                    if index >= self.maxReads:
                        break
                    localFile.write(record)
            utils.log("{} records written".format(index))
        remoteFile.close()
        localFile.close()
        baiFileName = sourceFileName + '.bai'
        os.remove(baiFileName)
        utils.log("Indexing '{}'".format(destFileName))
        pysam.index(destFileName.encode('utf-8'))

    def downloadBams(self):
        dirList = [self.dirName, "datasets", self.datasetName, 'reads']
        mkdirAndChdirList(dirList)
        for sample in self.samples:
            if isCheckpointFileExtant(sample):
                logCheckpoint(sample)
            else:
                self._downloadBam(sample)
                createCheckpointFile(sample)
        escapeDir(4)

    def _downloadFasta(self, chromosome):
        accession = self.accessions[chromosome]
        fileName = '{}.fa'.format(chromosome)
        minPos = 0
        if self.excludeReferenceMin:
            minPos = self.chromMinMax.getMinPos(chromosome)
        maxPos = self.chromMinMax.getMaxPos(chromosome)
        with open(fileName, "w") as outFasta:
            print(">{}".format(chromosome), file=outFasta)
            sequence = _fetchSequence(accession, minPos, maxPos)
            for line in sequence:
                print(line, file=outFasta)
        utils.log("Compressing {}".format(fileName))
        utils.runCommand("bgzip -f {}".format(fileName))
        compressedFileName = fileName + '.gz'
        utils.log("Indexing {}".format(compressedFileName))
        utils.runCommand("samtools faidx {}".format(compressedFileName))
        # Assemble the metadata.
        metadata = {
            "md5checksum": getReferenceChecksum(compressedFileName),
            "sourceUri": None,
            "ncbiTaxonId": 9606,
            "isDerived": False,
            "sourceDivergence": None,
            "sourceAccessions": [accession + ".subset"],
        }
        metadataFilename = "{}.json".format(chromosome)
        dumpDictToFileAsJson(metadata, metadataFilename)

    def downloadFastas(self):
        dirList = [self.dirName, 'referenceSets']
        mkdirAndChdirList(dirList)
        # Assemble reference set metadata
        referenceSetMetadata = {
            "assemblyId": 'TODO',
            "description": 'TODO',
            "isDerived": False,
            "ncbiTaxonId": 9606,
            "sourceAccessions": [],
            "sourceUri": 'TODO',
        }
        referenceSetMetadataFilename = "{}.json".format(
            self.referenceSetName)
        dumpDictToFileAsJson(
            referenceSetMetadata, referenceSetMetadataFilename)

        # Download chromosomes
        mkdirAndChdirList([self.referenceSetName])
        for chromosome in self.chromosomes:
            if isCheckpointFileExtant(chromosome):
                logCheckpoint(chromosome)
            else:
                self._downloadFasta(chromosome)
                createCheckpointFile(chromosome)
        escapeDir(3)

    def cleanup(self):
        utils.log('Removing checkpoint files')
        for root, dirs, files in os.walk(self.dirName):
            for currentFile in files:
                if currentFile.endswith('.checkpoint'):
                    os.remove(os.path.join(root, currentFile))


class NcbiFileDownloader(AbstractFileDownloader):
    """
    Downloads files from NCBI
    """
    def __init__(self, args):
        super(NcbiFileDownloader, self).__init__(args)

    def getBaseUrl(self):
        return 'ftp://ftp-trace.ncbi.nih.gov/1000genomes'


class EbiFileDownloader(AbstractFileDownloader):
    """
    Downloads files from EBI
    """
    def __init__(self, args):
        super(EbiFileDownloader, self).__init__(args)

    def getBaseUrl(self):
        return 'ftp://ftp.1000genomes.ebi.ac.uk/vol1'


sources = {
    "ncbi": NcbiFileDownloader,
    "ebi": EbiFileDownloader,
}


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source", default="ncbi", choices=sources.keys(),
        help="the source to download from")
    parser.add_argument(
        "--dir-name", default="ga4gh-downloaded-data",
        help="the name of the directory that the data is downloaded to")
    parser.add_argument(
        "--samples", default='HG00096,HG00533,HG00534',
        help="a comma-seperated list of samples to download")
    parser.add_argument(
        "--num-reads", "-r", default=1000, type=int,
        help="the number of reads to download per reference")
    parser.add_argument(
        "--num-variants", "-V", default=1000, type=int,
        help="the maximum number of variants to download per VCF file.")
    parser.add_argument(
        "--num-chromosomes", default="3", type=int,
        help=(
            "the number of chromosomes whose corresponding reads should "
            "be downloaded"))
    parser.add_argument(
        "--exclude-reference-min", default=False, action="store_true",
        help="Exclude bases in the reference before the minimum position")
    args = parser.parse_args()
    return args


@utils.Timed()
def main(args):
    downloaderClass = sources[args.source]
    downloader = downloaderClass(args)
    downloader.downloadVcfs()
    downloader.downloadBams()
    downloader.downloadFastas()
    downloader.cleanup()


if __name__ == '__main__':
    args = parseArgs()
    main(args)
