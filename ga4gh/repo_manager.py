"""
Manages a GA4GH data repository
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import hashlib
import json
import os
import shlex
import shutil
import subprocess

import pysam

import ga4gh.exceptions as exceptions
import ga4gh.datarepo as datarepo
import ga4gh.datamodel.datasets as datasets


def getReferenceChecksum(fastaFile):
    """
    Returns the md5 checksum for the reference sequence in the specified
    FASTA file. This is the MD5 of the upper case sequence letters.
    """
    inputFile = pysam.FastaFile(fastaFile)
    bases = inputFile.fetch(inputFile.references[0])
    inputFile.close()
    return hashlib.md5(bases.upper()).hexdigest()


def filenameWithoutExtension(filepath, extension):
    """
    Return a filename without the extension suffix

    (os.path.splitext(filename)[0] messes up
    on filenames with more than one period)
    """
    filename = os.path.basename(filepath)
    index = filename.index(extension)
    return filename[:index]


def runCommandSplits(splits, silent=False):
    """
    Run a shell command given the command's parsed command line
    """
    if silent:
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(splits, stdout=devnull, stderr=devnull)
    else:
        subprocess.check_call(splits)


def runCommand(command, silent=False):
    """
    Run a shell command
    """
    splits = shlex.split(command)
    runCommandSplits(splits, silent=silent)


class RepoManager(object):
    """
    Performs operations on a GA4GH data repository
    """
    datasetsDirName = datarepo.FileSystemDataRepository.datasetsDirName
    ontologiesDirName = datarepo.FileSystemDataRepository.ontologiesDirName
    referenceSetsDirName = \
        datarepo.FileSystemDataRepository.referenceSetsDirName
    readsDirName = datasets.FileSystemDataset.readsDirName
    variantsDirName = datasets.FileSystemDataset.variantsDirName
    fastaExtension = '.fa.gz'
    fastaIndexExtensionFai = '.fa.gz.fai'
    fastaIndexExtensionGzi = '.fa.gz.gzi'
    jsonExtension = '.json'
    vcfExtension = '.vcf.gz'
    vcfIndexExtension = '.vcf.gz.tbi'
    bamExtension = '.bam'
    bamIndexExtension = '.bam.bai'

    def __init__(self, repoPath):
        self._repoPath = repoPath
        self._topStructure = [
            self.datasetsDirName,
            self.ontologiesDirName,
            self.referenceSetsDirName]
        self._datasetStructure = [
            self.readsDirName, self.variantsDirName]

    def _assertFileExists(
            self, filePath, text='File', inRepo=False, emitName=None):
        if not os.path.exists(filePath):
            if emitName is None:
                emitName = filePath
            message = "{} '{}' does not exist".format(text, emitName)
            if inRepo:
                self._raiseException(message)
            else:
                raise exceptions.RepoManagerException(message)

    def _assertPathEmpty(
            self, path, text='Path', inRepo=False, emitName=None):
        if os.path.exists(path):
            if emitName is None:
                emitName = path
            message = "{} '{}' already exists".format(text, emitName)
            if inRepo:
                self._raiseException(message)
            else:
                raise exceptions.RepoManagerException(message)

    def _assertDirectory(
            self, dirPath, text='File', inRepo=False, emitName=None):
        if not os.path.isdir(dirPath):
            if emitName is None:
                emitName = dirPath
            message = "{} '{}' is not a directory".format(text, emitName)
            if inRepo:
                self._raiseException(message)
            else:
                raise exceptions.RepoManagerException(message)

    def _emit(self, message):
        self._emitIndent(message, 0)

    def _emitIndent(self, message, indentLevel=1):
        print('\t' * indentLevel + message)

    def _repoEmit(self, message):
        header = "Repo at '{}'".format(self._repoPath)
        self._emit(header)
        self._emit(message)

    def _raiseException(self, message):
        exceptionMessage = "Exception for repo at '{}'\n{}".format(
            self._repoPath, message)
        raise exceptions.RepoManagerException(exceptionMessage)

    def _getDatasetPath(self, datasetName):
        datasetPath = os.path.join(
            self._repoPath, self.datasetsDirName, datasetName)
        return datasetPath

    def _getReferenceSetPath(self, referenceSetName):
        referenceSetPath = os.path.join(
            self._repoPath, self.referenceSetsDirName, referenceSetName)
        return referenceSetPath

    def _getReferenceSetJsonPath(self, referenceSetName):
        jsonPath = os.path.join(
            self._repoPath, self.referenceSetsDirName,
            referenceSetName + self.jsonExtension)
        return jsonPath

    def _getReadGroupSetPath(self, datasetName, readGroupSetName):
        readGroupSetPath = os.path.join(
            self._repoPath, self.datasetsDirName, datasetName,
            self.readsDirName, readGroupSetName) + self.bamExtension
        return readGroupSetPath

    def _getReadGroupSetIndexPath(self, datasetName, readGroupSetName):
        indexPath = os.path.join(
            self._repoPath, self.datasetsDirName, datasetName,
            self.readsDirName, readGroupSetName)
        return indexPath + self.bamIndexExtension

    def _getVariantSetPath(self, datasetName, variantSetName):
        variantSetPath = os.path.join(
            self._repoPath, self.datasetsDirName, datasetName,
            self.variantsDirName, variantSetName)
        return variantSetPath

    def _moveFile(self, src, dst, mode):
        if mode == 'move':
            os.rename(src, dst)
        elif mode == 'copy':
            if os.path.isdir(src):
                shutil.copytree(src, dst)
            else:
                shutil.copy(src, dst)
        elif mode == 'link':
            os.symlink(os.path.abspath(src), os.path.abspath(dst))
            os.stat(dst)  # check that the link is not broken
        else:
            self._raiseException(
                "Unrecognized move file mode '{}'".format(mode))

    def _createJsonFile(self, jsonFilePath, metadata):
        with open(jsonFilePath, 'w') as jsonFile:
            json.dump(metadata, jsonFile, indent=4)

    def _removePath(self, path):
        if os.path.isdir(path):
            if os.path.islink(path):
                os.unlink(path)
            else:
                shutil.rmtree(path)
        else:
            os.unlink(path)

    def _checkFile(self, filePath, fileExt):
        self._assertFileExists(filePath)
        if not filePath.endswith(fileExt):
            raise exceptions.RepoManagerException(
                "File '{}' does not have a '{}' extension".format(
                    filePath, fileExt))

    def _check(self):
        if not os.path.exists(self._repoPath):
            self._raiseException("No repo at path")
        if not os.path.isdir(self._repoPath):
            self._raiseException("Another file at path")
        topLevelDirs = os.listdir(self._repoPath)
        for topDir in self._topStructure:
            if topDir not in topLevelDirs:
                self._raiseException(
                    "top-level directory does not contain required "
                    "directory '{}'".format(topDir))
        # TODO more checks here ...
        # or perhaps just init a FileSystemDataRepository

    def _checkDataset(self, datasetName):
        datasetDir = self._getDatasetPath(datasetName)
        self._assertFileExists(datasetDir, 'Dataset', True, datasetName)
        self._assertDirectory(datasetDir, 'Dataset', True, datasetName)
        for directory in self._datasetStructure:
            datasetSubDir = os.path.join(datasetDir, directory)
            self._assertFileExists(datasetSubDir, 'Dataset', True)
            self._assertDirectory(datasetSubDir, 'Dataset', True)

    def _checkFolder(self, folderPath, fileExt):
        self._assertFileExists(folderPath, 'Directory')
        self._assertDirectory(folderPath)
        folderFiles = os.listdir(folderPath)
        if len(folderFiles) == 0:
            raise exceptions.RepoManagerException(
                "Directory '{}' has no contents; "
                "needs files with extension '{}'".format(
                    folderPath, fileExt))
        vcfPresent = any(
            [folderFile.endswith(fileExt) for folderFile in folderFiles])
        if not vcfPresent:
            raise exceptions.RepoManagerException(
                "Directory '{}' does not contain a file "
                "with a '{}' extension".format(
                    folderPath, fileExt))

    ############
    # Commands #
    ############

    def init(self):
        """
        Initialize a repository
        """
        if os.path.exists(self._repoPath):
            self._assertDirectory(self._repoPath)
            raise exceptions.RepoManagerException(
                "Directory already exists at '{}'".format(self._repoPath))
        else:
            os.mkdir(self._repoPath)
            for directory in self._topStructure:
                newDir = os.path.join(self._repoPath, directory)
                os.mkdir(newDir)
            self._repoEmit("Created")

    def destroy(self):
        """
        Removes the repository
        """
        self._check()
        self._removePath(self._repoPath)
        self._repoEmit("Destroyed")

    def check(self):
        """
        Check the repository for well-formedness
        """
        self._check()
        self._repoEmit("Well-formed".format(self._repoPath))

    def addDataset(self, datasetName):
        """
        Add a dataset to the repository
        """
        self._check()
        datasetPath = self._getDatasetPath(datasetName)
        self._assertPathEmpty(datasetPath, 'Dataset', datasetName)
        os.mkdir(datasetPath)
        for directory in self._datasetStructure:
            newDir = os.path.join(datasetPath, directory)
            os.mkdir(newDir)
        self._repoEmit("Dataset '{}' added".format(datasetName))

    def removeDataset(self, datasetName):
        """
        Remove a dataset from the repository
        """
        self._check()
        self._checkDataset(datasetName)
        datasetPath = self._getDatasetPath(datasetName)
        self._removePath(datasetPath)
        self._repoEmit("Dataset '{}' removed".format(datasetName))

    def addReferenceSet(self, filePath, moveMode, metadata):
        """
        Add a reference set to the repo
        """
        # move the fasta file
        self._check()
        self._checkFile(filePath, self.fastaExtension)
        fileName = filenameWithoutExtension(filePath, self.fastaExtension)
        destPath = os.path.join(
            self._repoPath, self.referenceSetsDirName, fileName)
        self._assertPathEmpty(destPath, inRepo=True)
        os.mkdir(destPath)
        fileDestPath = os.path.join(destPath, os.path.basename(filePath))
        self._moveFile(filePath, fileDestPath, moveMode)

        # move the index files if they exist, otherwise do indexing
        indexPathFai = os.path.join(
            os.path.split(
                filePath)[0], fileName + self.fastaIndexExtensionFai)
        indexPathGzi = os.path.join(
            os.path.split(
                filePath)[0], fileName + self.fastaIndexExtensionGzi)
        indexedMessage = ""
        if os.path.exists(indexPathFai) and os.path.exists(indexPathGzi):
            self._moveFile(
                indexPathFai,
                os.path.join(destPath, os.path.basename(indexPathFai)),
                moveMode)
            self._moveFile(
                indexPathGzi,
                os.path.join(destPath, os.path.basename(indexPathGzi)),
                moveMode)
        else:
            runCommand("samtools faidx {}".format(fileDestPath))
            indexedMessage = " (and indexed)"

        # create reference set json file
        referenceSetJsonFileName = fileName + self.jsonExtension
        referenceSetJsonFilePath = os.path.join(
            self._repoPath, self.referenceSetsDirName,
            referenceSetJsonFileName)
        referenceSetMetadata = {
            "assemblyId": 'TODO',
            "description": metadata.get('description', 'TODO'),
            "isDerived": False,
            "ncbiTaxonId": 9606,
            "sourceAccessions": [],
            "sourceUri": 'TODO',
        }
        self._createJsonFile(referenceSetJsonFilePath, referenceSetMetadata)

        # create reference json file
        referenceJsonFilePath = os.path.join(
            destPath, fileName + self.jsonExtension)
        md5checksum = getReferenceChecksum(fileDestPath)
        referenceMetadata = {
            "md5checksum": md5checksum,
            "sourceUri": None,
            "ncbiTaxonId": 9606,
            "isDerived": False,
            "sourceDivergence": None,
            "sourceAccessions": [],
        }
        self._createJsonFile(referenceJsonFilePath, referenceMetadata)

        # finish
        self._repoEmit("ReferenceSet '{}' added{}".format(
            fileName, indexedMessage))

    def removeReferenceSet(self, referenceSetName):
        self._check()
        referenceSetPath = self._getReferenceSetPath(referenceSetName)
        self._assertDirectory(referenceSetPath)
        self._removePath(referenceSetPath)
        jsonPath = self._getReferenceSetJsonPath(referenceSetName)
        self._removePath(jsonPath)
        self._repoEmit("ReferenceSet '{}' removed".format(
            referenceSetName))

    def addReadGroupSet(self, datasetName, filePath, moveMode):
        """
        Add a read group set to the repo
        """
        # move the bam file
        self._check()
        self._checkDataset(datasetName)
        self._checkFile(filePath, self.bamExtension)
        fileName = os.path.basename(filePath)
        readGroupSetName = filenameWithoutExtension(
            fileName, self.bamExtension)
        destPath = os.path.join(
            self._repoPath, self.datasetsDirName, datasetName,
            self.readsDirName, fileName)
        self._assertPathEmpty(destPath, inRepo=True)
        self._moveFile(filePath, destPath, moveMode)

        # move the index file if it exists, otherwise do indexing
        indexPath = os.path.join(
            os.path.split(filePath)[0],
            readGroupSetName + self.bamIndexExtension)
        indexMessage = ""
        if os.path.exists(indexPath):
            dstDir = os.path.split(destPath)[0]
            self._moveFile(
                indexPath,
                os.path.join(dstDir, os.path.basename(indexPath)),
                moveMode)
        else:
            pysam.index(destPath.encode('utf-8'))
            indexMessage = " (and indexed)"

        # finish
        self._repoEmit("ReadGroupSet '{}' added to dataset '{}'{}".format(
            fileName, datasetName, indexMessage))

    def removeReadGroupSet(self, datasetName, readGroupSetName):
        """
        Remove a read group set from the repo
        """
        self._check()
        self._checkDataset(datasetName)
        readGroupSetPath = self._getReadGroupSetPath(
            datasetName, readGroupSetName)
        self._assertFileExists(readGroupSetPath, inRepo=True)
        self._removePath(readGroupSetPath)
        indexPath = self._getReadGroupSetIndexPath(
            datasetName, readGroupSetName)
        self._assertFileExists(indexPath, inRepo=True)
        self._removePath(indexPath)
        self._repoEmit("ReadGroupSet '{}/{}' removed".format(
            datasetName, readGroupSetName))

    def addVariantSet(self, datasetName, filePath, moveMode):
        """
        Add a variant set to the repo
        """
        # move the vcf file
        self._check()
        self._checkDataset(datasetName)
        self._checkFolder(filePath, self.vcfExtension)
        dirName = os.path.basename(filePath)
        destPath = os.path.join(
            self._repoPath, self.datasetsDirName, datasetName,
            self.variantsDirName, dirName)
        self._assertPathEmpty(destPath, inRepo=True)
        self._moveFile(filePath, destPath, moveMode)

        # do indexing if the indexes do not exist
        indexCount = 0
        dirFiles = os.listdir(destPath)
        for dirFile in dirFiles:
            if dirFile.endswith(self.vcfExtension):
                vcfName = filenameWithoutExtension(
                    dirFile, self.vcfExtension)
                indexName = vcfName + self.vcfIndexExtension
                indexPath = os.path.join(destPath, indexName)
                if not os.path.exists(indexPath):
                    vcfPath = os.path.join(destPath, dirFile)
                    runCommand('tabix {}'.format(vcfPath))
                    indexCount += 1

        # finish
        indexMessage = ""
        if indexCount > 0:
            indexMessage = " ({} vcfs indexed)".format(indexCount)
        self._repoEmit("VariantSet '{}' added to dataset '{}'{}".format(
            dirName, datasetName, indexMessage))

    def removeVariantSet(self, datasetName, variantSetName):
        """
        Remove a variant set from the repo
        """
        self._check()
        self._checkDataset(datasetName)
        variantSetPath = self._getVariantSetPath(
            datasetName, variantSetName)
        self._assertDirectory(variantSetPath)
        self._removePath(variantSetPath)
        self._repoEmit("Variant set '{}/{}' removed".format(
            datasetName, variantSetName))

    def list(self):
        """
        List the contents of the repo
        """
        self._check()
        self._repoEmit("Listing")

        dataRepo = datarepo.FileSystemDataRepository(
            self._repoPath, doConsistencyCheck=False)
        self._emit(self.referenceSetsDirName)
        for referenceSet in dataRepo.getReferenceSets():
            self._emitIndent(referenceSet.getLocalId())
            for reference in referenceSet.getReferences():
                self._emitIndent(reference.getLocalId(), 2)
        self._emit(self.datasetsDirName)
        for dataset in dataRepo.getDatasets():
            self._emitIndent(dataset.getLocalId())
            self._emitIndent(self.readsDirName, 2)
            for readGroupSet in dataset.getReadGroupSets():
                self._emitIndent(readGroupSet.getLocalId(), 3)
            self._emitIndent(self.variantsDirName, 2)
            for variantSet in dataset.getVariantSets():
                self._emitIndent(variantSet.getLocalId(), 3)
                for chromFile in sorted(variantSet._chromFileMap.keys()):
                    self._emitIndent(chromFile, 4)
