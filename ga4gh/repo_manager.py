"""
Manages a GA4GH data repository
"""
import os


class RepoManagerException(Exception):
    pass


class RepoManager(object):

    datasetsDirName = "datasets"
    referenceSetsDirName = "referenceSets"
    readsDirName = "reads"
    variantsDirName = "variants"

    def __init__(self, repoPath):
        self.repoPath = repoPath
        self.topStructure = [
            self.datasetsDirName, self.referenceSetsDirName]
        self.datasetStructure = [
            self.readsDirName, self.variantsDirName]

    def _emit(self, message):
        print(message)

    def _repoEmit(self, message):
        header = "Repo at '{}'".format(self.repoPath)
        self._emit(header)
        self._emit(message)

    def _raiseException(self, message):
        exceptionMessage = "Exception for repo at '{}'\n{}".format(
            self.repoPath, message)
        raise RepoManagerException(exceptionMessage)

    def init(self):
        if os.path.exists(self.repoPath):
            if os.path.isdir(self.repoPath):
                raise RepoManagerException(
                    "Directory already exists at '{}'".format(self.repoPath))
            else:
                raise RepoManagerException(
                    "File at '{}' is not a directory".format(self.repoPath))
        else:
            os.mkdir(self.repoPath)
            for directory in self.topStructure:
                newDir = os.path.join(self.repoPath, directory)
                os.mkdir(newDir)
            self._repoEmit("Created")

    def check(self):
        self._check()
        self._repoEmit("Well-formed".format(self.repoPath))

    def _check(self):
        if not os.path.exists(self.repoPath):
            self._raiseException("No repo at path")
        if not os.path.isdir(self.repoPath):
            self._raiseException("Another file at path")
        topLevelDirs = os.listdir(self.repoPath)
        if len(topLevelDirs) != len(self.topStructure):
            self._raiseException(
                "top-level directory contents should be '{}' "
                "but are actually '{}'".format(
                    self.topStructure, topLevelDirs))
        for directory in topLevelDirs:
            if directory not in self.topStructure:
                self._raiseException(
                    "top-level directory does not contain required "
                    "directory '{}'".format(directory))
        # TODO more checks here

    def addDataset(self, datasetName):
        self._check()
        datasetDir = os.path.join(
            self.repoPath, self.datasetsDirName, datasetName)
        os.mkdir(datasetDir)
        for directory in self.datasetStructure:
            newDir = os.path.join(datasetDir, directory)
            os.mkdir(newDir)
        self._repoEmit("Dataset '{}' added".format(datasetName))

    def addReferenceSet(self, filePath):
        self._check()
        # TODO

    def _checkDataset(self, datasetName):
        datasetDir = os.path.join(
            self.repoPath, self.datasetsDirName, datasetName)
        if not os.path.exists(datasetDir):
            self._raiseException(
                "Dataset '{}' does not exist".format(datasetName))
        if not os.path.isdir(datasetDir):
            self._raiseException(
                "Dataset '{}' is not a directory".format(datasetName))

    def _checkFile(self, datasetName, filePath, fileExt):
        if not os.path.exists(filePath):
            raise RepoManagerException(
                "File '{}' does not exist".format(filePath))
        fileName = os.path.basename(filePath)
        extension = os.path.splitext(fileName)[1]
        if extension != '.' + fileExt:
            raise RepoManagerException(
                "File '{}' does not have a '{}' extension".format(
                    filePath, fileExt))

    def addReadGroupSet(self, datasetName, filePath):
        self._check()
        self._checkDataset(datasetName)
        self._checkFile(datasetName, filePath, 'bam')
        fileName = os.path.basename(filePath)
        destPath = os.path.join(
            self.repoPath, self.datasetsDirName, datasetName,
            self.readsDirName, fileName)
        if os.path.exists(destPath):
            self._raiseException(
                "Path '{}' already exists".format(destPath))
        os.rename(filePath, destPath)  # TODO add option for copy / link
        # TODO index bam
        self._repoEmit("ReadGroupSet '{}' added to dataset '{}'".format(
            fileName, datasetName))

    def addVariantSet(self, datasetName, filePath):
        self._check()
        # TODO

    def list(self):
        self._check()
        # TODO
