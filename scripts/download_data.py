"""
Constructs a data source for the ga4gh server by downloading data from
authoritative remote servers.
"""
# TODO
# - need some kind of checkpoint functionality to resume process
#   where it left off since getting a clean run is so rare...
# - references support
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import glob
import gzip
import os
import tempfile
import urllib
import urllib2

import pysam

import utils


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


def cleanDir():
    """
    Removes genomic files in current directory
    """
    cwd = os.getcwd()
    utils.log("Cleaning out directory '{}'".format(cwd))
    globs = [
        "*.tbi", "*.vcf", "*.vcf.gz", "*.bam", "*.bam.bai", "*.fa.gz",
        "*.fa", "*.fa.gz.fai", "*.fa.gz.gzi"]
    for fileGlob in globs:
        fileNames = glob.glob(fileGlob)
        for fileName in fileNames:
            os.remove(fileName)


def escapeDir(levels=4):
    # back to orig dir
    for _ in range(levels):
        os.chdir('..')


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


class AbstractFileDownloader(object):
    """
    Base class for individual site genome file downloaders
    """
    def __init__(self, args):
        self.args = args
        self.datasetName = 'dataset1'
        self.variantSetName = self.args.source
        self.chromMinMax = ChromMinMax()
        self.chromosomes = self.args.chromosomes.split(',')
        self.accessions = {
            '1': 'CM000663.2',
            '2': 'CM000664.2',
            '3': 'CM000665.2',
            '4': 'CM000666.2',
            '5': 'CM000667.2',
            '6': 'CM000668.2',
            '7': 'CM000669.2',
            '8': 'CM000670.2',
            '9': 'CM000671.2',
            '10': 'CM000672.2',
            '11': 'CM000673.2',
            '12': 'CM000674.2',
            '13': 'CM000675.2',
            '14': 'CM000676.2',
            '15': 'CM000677.2',
            '16': 'CM000678.2',
            '17': 'CM000679.2',
            '18': 'CM000680.2',
            '19': 'CM000681.2',
            '20': 'CM000682.2',
            '21': 'CM000683.2',
            '22': 'CM000684.2',
            'X': 'CM000685.2',
            'Y': 'CM000686.2',
        }

    def _getVcfFilenames(self):
        baseFileName = (
            "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a"
            ".20130502.genotypes.vcf.gz")
        chrNames = [
            str(i) for i in range(1, 23)
            if str(i) in self.chromosomes]
        fileNames = [baseFileName.format(chrName) for chrName in chrNames]
        # the X and Y files use a different filename prefix
        if 'X' in self.chromosomes:
            fileNames.append(
                'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.'
                '20130502.genotypes.vcf.gz')
        if 'Y' in self.chromosomes:
            fileNames.append(
                'ALL.chrY.phase3_integrated_v1a.'
                '20130502.genotypes.vcf.gz')
        return fileNames

    def getVcfBaseUrl(self):
        return os.path.join(self.getBaseUrl(), 'ftp/release/20130502/')

    def getBamBaseUrl(self):
        return os.path.join(self.getBaseUrl(), 'ftp/phase3/data/')

    def _prepareDir(self):
        dirList = [
            self.args.dir_name, self.datasetName, 'variants',
            self.variantSetName]
        mkdirAndChdirList(dirList)
        cleanDir()

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

    def _processFileName(self, fileName):
        url = os.path.join(self.getVcfBaseUrl(), fileName)
        utils.log("Downloading '{}'".format(url))
        response = urllib2.urlopen(url)
        megabyte = 1024 * 1024
        data = response.read(megabyte)
        lineCountQuota = 1000
        localFileName, _ = os.path.splitext(fileName)
        utils.log("Writing '{}'".format(localFileName))
        with tempfile.NamedTemporaryFile() as binaryFile:
            binaryFile.write(data)
            binaryFile.flush()
            gzipFile = gzip.open(binaryFile.name, "r")
            outputFile = open(localFileName, "w")
            lineCount = 0
            for line in gzipFile:
                outputFile.write(line)
                if not line.startswith("#"):
                    lineCount += 1
                if lineCount >= lineCountQuota:
                    break
            assert lineCount == lineCountQuota
            outputFile.close()
            gzipFile.close()
        utils.log("Compressing '{}'".format(localFileName))
        utils.runCommand('bgzip -f {}'.format(localFileName))
        utils.log("Indexing '{}'".format(fileName))
        utils.runCommand('tabix {}'.format(fileName))
        self._updatePositions(fileName)

    def downloadVcfs(self):
        self._prepareDir()
        fileNames = self._getVcfFilenames()
        for fileName in fileNames:
            self._processFileName(fileName)
        escapeDir()

    def downloadBams(self):
        dirList = [
            self.args.dir_name, self.datasetName, 'reads', 'low-coverage']
        mkdirAndChdirList(dirList)
        cleanDir()
        studyMap = {
            'HG00096': 'GBR',
            'HG00533': 'CHS',
            'HG00534': 'CHS',
        }
        baseUrl = self.getBamBaseUrl()
        samples = self.args.samples.split(',')
        for sample in samples:
            samplePath = '{}/alignment/'.format(sample)
            study = studyMap[sample]
            fileName = (
                '{}.mapped.ILLUMINA.bwa.{}.'
                'low_coverage.20120522.bam'.format(sample, study))
            sampleUrl = os.path.join(baseUrl, samplePath, fileName)
            utils.log("Downloading index for '{}'".format(sampleUrl))
            remoteFile = pysam.AlignmentFile(sampleUrl)
            header = remoteFile.header
            utils.log("Writing '{}'".format(fileName))
            localFile = pysam.AlignmentFile(
                fileName, 'wb', header=header)
            for chromosome in self.chromosomes:
                utils.log("chromosome {}".format(chromosome))
                iterator = remoteFile.fetch(
                    chromosome.encode('utf-8'),
                    start=self.chromMinMax.getMinPos(chromosome),
                    end=self.chromMinMax.getMaxPos(chromosome))
                for index, record in enumerate(iterator):
                    if index >= args.num_reads:
                        break
                    localFile.write(record)
                utils.log("{} records written".format(index))
            iterator = None
            remoteFile.close()
            localFile.close()
            baiFileName = fileName + '.bai'
            os.remove(baiFileName)
            utils.log("Indexing '{}'".format(fileName))
            pysam.index(fileName.encode('utf-8'))
        escapeDir()

    def downloadFastas(self):
        dirList = [
            self.args.dir_name, 'references', 'ebi']
        mkdirAndChdirList(dirList)
        cleanDir()
        baseUrl = 'http://www.ebi.ac.uk/ena/data/view/'
        for chromosome in self.chromosomes:
            accession = self.accessions[chromosome]
            path = os.path.join(baseUrl, accession)
            maxPos = self.chromMinMax.getMaxPos(chromosome)
            args = urllib.urlencode({
                'display': 'fasta',
                'range': '{}-{}'.format(0, maxPos)})
            url = '{}%26{}'.format(path, args)
            fileName = '{}.fa'.format(chromosome)
            downloader = utils.HttpFileDownloader(url, fileName)
            downloader.download()
            utils.log("Compressing {}".format(fileName))
            utils.runCommand("bgzip {}".format(fileName))
            compressedFileName = fileName + '.gz'
            utils.log("Indexing {}".format(compressedFileName))
            utils.runCommand("samtools faidx {}".format(compressedFileName))
        escapeDir(3)


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
        "--num-reads", default=1000,
        help="the number of reads to download per reference")
    parser.add_argument(
        "--chromosomes", default="1,2,3,X,Y",
        help="the chromosomes whose corresponding reads should be downloaded")
    args = parser.parse_args()
    return args


@utils.Timed()
def main(args):
    downloaderClass = sources[args.source]
    downloader = downloaderClass(args)
    downloader.downloadVcfs()
    downloader.downloadBams()
    downloader.downloadFastas()


if __name__ == '__main__':
    args = parseArgs()
    main(args)
