"""
Constructs a data source for the ga4gh server by downloading data from
authoritative remote servers.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import gzip
import os
import requests
import shutil
import subprocess
import tempfile
import urllib2

import pysam

import ga4gh.common.utils as utils
import glue

glue.ga4ghImportGlue()

# We need to turn off QA because of the import glue
import ga4gh.server.datarepo as datarepo  # NOQA
import ga4gh.server.datamodel.references as references  # NOQA
import ga4gh.server.datamodel.datasets as datasets  # NOQA
import ga4gh.server.datamodel.variants as variants  # NOQA
import ga4gh.server.datamodel.reads as reads  # NOQA


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
    """Fetch sequences from NCBI using the eself interface.

    An interbase interval may be optionally provided with startIndex and
    endIndex. NCBI eself will return just the requested subsequence, which
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
    urlFmt = (
        "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=nucleotide&id={ac}&rettype=fasta&retmode=text")
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
        self.tempDir = tempfile.mkdtemp(prefix="ga4gh-download")
        self.numChromosomes = args.num_chromosomes
        self.chromosomes = [str(j + 1) for j in range(self.numChromosomes)]
        self.dirName = args.destination
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
        self.vcfFilePaths = []
        self.bamFilePaths = []
        self.fastaFilePath = None
        if os.path.exists(self.dirName) and args.force:
            shutil.rmtree(self.dirName)
        os.mkdir(self.dirName)
        self.repoPath = os.path.join(self.dirName, "registry.db")

    def log(self, message):
        print(message)

    def runCommand(self, command):
        """
        Runs the specified command.
        """
        subprocess.check_call(command, shell=True)

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
        self.log('chrom: {}, maxPos: {}, minPos: {}'.format(
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
        self.log("Downloading '{}'".format(url))
        response = urllib2.urlopen(url)
        megabyte = 1024 * 1024
        data = response.read(megabyte)
        localFileName = os.path.join(
            self.dirName, "chr{}.vcf".format(chromosome))
        localCompressedFileName = "{}.gz".format(localFileName)
        localTempFileName = localFileName + '.unsampled'
        self.log("Writing '{}'".format(localTempFileName))
        self._writeVcfTempFile(localTempFileName, data)
        self.log("Sampling '{}'".format(localTempFileName))
        self.runCommand(
            'bcftools view --force-samples -s {} {} -o {}'.format(
                args.samples, localTempFileName, localFileName))
        os.remove(localTempFileName)
        self.log("Compressing '{}'".format(localFileName))
        self.runCommand('bgzip -f {}'.format(localFileName))
        self.log("Indexing '{}'".format(localCompressedFileName))
        self.runCommand('tabix {}'.format(localCompressedFileName))
        self._updatePositions(localCompressedFileName)
        self.vcfFilePaths.append(
            (localCompressedFileName, localCompressedFileName + ".tbi"))

    def downloadVcfs(self):
        for chromosome in self.chromosomes:
            self._downloadVcf(chromosome)

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

    def _downloadIndex(self, indexUrl, localIndexFile):
        self.log("Downloading index from {} to {}".format(
            indexUrl, localIndexFile))
        response = urllib2.urlopen(indexUrl)
        with open(localIndexFile, "w") as destFile:
            destFile.write(response.read())

    def _downloadBam(self, sample):
        samplePath = '{}/alignment/'.format(sample)
        study = self.studyMap[sample]
        sourceFileName = (
            '{}.mapped.ILLUMINA.bwa.{}.'
            'low_coverage.20120522.bam'.format(sample, study))
        destFileName = os.path.join(
            self.dirName, "{}.bam".format(sample))
        baseUrl = self.getBamBaseUrl()
        sampleUrl = os.path.join(baseUrl, samplePath, sourceFileName)
        indexUrl = sampleUrl + ".bai"
        localIndexFile = os.path.join(self.tempDir, sourceFileName + ".bai")
        self._downloadIndex(indexUrl, localIndexFile)
        remoteFile = pysam.AlignmentFile(
            sampleUrl, filepath_index=localIndexFile)
        header = self.createBamHeader(remoteFile.header)
        self.log("Writing '{}'".format(destFileName))
        localFile = pysam.AlignmentFile(
            destFileName, 'wb', header=header)
        for chromosome in self.chromosomes:
            self.log("chromosome {}".format(chromosome))
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
            self.log("{} records written".format(index))
        remoteFile.close()
        localFile.close()
        self.log("Indexing '{}'".format(destFileName))
        pysam.index(destFileName.encode('utf-8'))
        self.bamFilePaths.append(
            (destFileName, destFileName + ".bai"))

    def downloadBams(self):
        for sample in self.samples:
            self._downloadBam(sample)

    def _downloadFasta(self, chromosomes):
        fileName = os.path.join(self.dirName, "GRCh37-subset.fa")
        with open(fileName, "w") as outFasta:
            for chromosome in chromosomes:
                accession = self.accessions[chromosome]
                minPos = 0
                if self.excludeReferenceMin:
                    minPos = self.chromMinMax.getMinPos(chromosome)
                maxPos = self.chromMinMax.getMaxPos(chromosome)
                print(minPos, maxPos)
                print(">{}".format(chromosome), file=outFasta)
                sequence = _fetchSequence(accession, minPos, maxPos)
                for line in sequence:
                    print(line, file=outFasta)
        self.log("Compressing {}".format(fileName))
        self.runCommand("bgzip -f {}".format(fileName))
        compressedFileName = fileName + '.gz'
        self.log("Indexing {}".format(compressedFileName))
        self.runCommand("samtools faidx {}".format(compressedFileName))
        self.fastaFilePath = compressedFileName

    def downloadReference(self):
        self._downloadFasta(self.chromosomes)

    def createRepo(self):
        """
        Creates the repository for all the data we've just downloaded.
        """
        repo = datarepo.SqlDataRepository(self.repoPath)
        repo.open("w")
        repo.initialise()

        referenceSet = references.HtslibReferenceSet("GRCh37-subset")
        referenceSet.populateFromFile(self.fastaFilePath)
        referenceSet.setDescription("Subset of GRCh37 used for demonstration")
        referenceSet.setSpeciesFromJson(
                '{"id": "9606",'
                + '"term": "Homo sapiens", "source_name": "NCBI"}')
        for reference in referenceSet.getReferences():
            reference.setSpeciesFromJson(
                '{"id": "9606",'
                + '"term": "Homo sapiens", "source_name": "NCBI"}')
            reference.setSourceAccessions(
                self.accessions[reference.getName()] + ".subset")
        repo.insertReferenceSet(referenceSet)

        dataset = datasets.Dataset("1kg-p3-subset")
        dataset.setDescription("Sample data from 1000 Genomes phase 3")
        repo.insertDataset(dataset)

        variantSet = variants.HtslibVariantSet(dataset, "mvncall")
        variantSet.setReferenceSet(referenceSet)
        dataUrls = [vcfFile for vcfFile, _ in self.vcfFilePaths]
        indexFiles = [indexFile for _, indexFile in self.vcfFilePaths]
        variantSet.populateFromFile(dataUrls, indexFiles)
        variantSet.checkConsistency()
        repo.insertVariantSet(variantSet)

        for sample, (bamFile, indexFile) in zip(
                self.samples, self.bamFilePaths):
            readGroupSet = reads.HtslibReadGroupSet(dataset, sample)
            readGroupSet.populateFromFile(bamFile, indexFile)
            readGroupSet.setReferenceSet(referenceSet)
            repo.insertReadGroupSet(readGroupSet)

        repo.commit()
        repo.close()
        self.log("Finished creating the repository; summary:\n")
        repo.open("r")
        repo.printSummary()

    def cleanup(self):
        self.log('Removing temporary files')
        shutil.rmtree(self.tempDir)


class NcbiFileDownloader(AbstractFileDownloader):
    """
    Downloads files from NCBI
    """
    def getBaseUrl(self):
        return 'ftp://ftp-trace.ncbi.nih.gov/1000genomes'


class EbiFileDownloader(AbstractFileDownloader):
    """
    Downloads files from EBI
    """
    def getBaseUrl(self):
        return 'ftp://ftp.1000genomes.ebi.ac.uk/vol1'


sources = {
    "ncbi": NcbiFileDownloader,
    "ebi": EbiFileDownloader,
}


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "destination",
        help="the name of the directory that the data is downloaded to")
    parser.add_argument(
        "-f", "--force", default=False, action="store_true",
        help="Overwrite an existing directory with the same name")
    parser.add_argument(
        "--source", default="ncbi", choices=sources.keys(),
        help="the source to download from")
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
    requiredExecutables = ['bcftools', 'tabix', 'bgzip', 'samtools']
    utils.requireExecutables(requiredExecutables)
    downloaderClass = sources[args.source]
    downloader = downloaderClass(args)
    try:
        downloader.downloadVcfs()
        downloader.downloadReference()
        downloader.downloadBams()
        downloader.createRepo()
    finally:
        downloader.cleanup()


if __name__ == '__main__':
    args = parseArgs()
    main(args)
