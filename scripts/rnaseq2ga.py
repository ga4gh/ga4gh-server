"""
Script to parse the output file produced by cufflinks.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sqlite3
import argparse

import utils

utils.ga4ghImportGlue()

import ga4gh.datarepo as datarepo  # NOQA


def data_repo(path):
    dataRepository = datarepo.SqlDataRepository(path)
    dataRepository.open(datarepo.MODE_READ)
    return dataRepository


class RNASqliteStore(object):
    """
    Defines a sqlite store for RNA data as well as methods for loading the
    tables.
    """
    def __init__(self, sqliteFileName, repoPath):
        # TODO: check to see if the rnaQuantId is in the db and exit if it is
        # since this is a generator and not an updater
        sqlFilePath = sqliteFileName
        self._dbConn = sqlite3.connect(sqlFilePath)
        self._cursor = self._dbConn.cursor()
        self.createTables(self._cursor)
        self._dbConn.commit()

        self._batchSize = 100
        self._rnaValueList = []
        self._expressionValueList = []

    def createTables(self, cursor):
        # annotationIds is a comma separated list
        cursor.execute('''CREATE TABLE RnaQuantification (
                       id text,
                       feature_set_ids text,
                       description text,
                       name text,
                       read_group_ids text,
                       programs text)''')
        cursor.execute('''CREATE TABLE Expression (
                       id text,
                       rna_quantification_id text,
                       name text,
                       feature_id text,
                       expression real,
                       is_normalized boolean,
                       raw_read_count real,
                       score real,
                       units integer,
                       conf_low real,
                       conf_hi real)''')

    def addRNAQuantification(self, datafields):
        """
        Adds an RNAQuantification to the db.  Datafields is a tuple in the
        order:
        id, feature_set_ids, description, name, read_group_ids, programs
        """
        self._rnaValueList.append(datafields)
        if len(self._rnaValueList) >= self._batchSize:
            self.batchAddRNAQuantification()

    def batchAddRNAQuantification(self):
        if len(self._rnaValueList) > 0:
            sql = "INSERT INTO RnaQuantification VALUES (?,?,?,?,?,?)"
            self._cursor.executemany(sql, self._rnaValueList)
            self._dbConn.commit()
            self._rnaValueList = []

    def addExpression(self, datafields):
        """
        Adds an Expression to the db.  Datafields is a tuple in the order:
        id, rna_quantification_id, name, feature_id, expression,
        is_normalized, raw_read_count, score, units, conf_low, conf_hi
        """
        self._expressionValueList.append(datafields)
        if len(self._expressionValueList) >= self._batchSize:
            self.batchAddExpression()

    def batchAddExpression(self):
        if len(self._expressionValueList) > 0:
            sql = "INSERT INTO Expression VALUES (?,?,?,?,?,?,?,?,?,?,?)"
            self._cursor.executemany(sql, self._expressionValueList)
            self._dbConn.commit()
            self._expressionValueList = []


class AbstractWriter(object):
    """
    Base class to use for the rna quantification writers
    """
    def __init__(self, rnaDB, repoPath, featureType="gene"):
        self._db = rnaDB
        self._isNormalized = None
        self._units = 0  # EXPRESSION_UNIT_UNSPECIFIED
        self._expressionLevelCol = None
        self._idCol = None
        self._nameCol = None
        self._featureCol = None
        self._countCol = None
        self._confColLow = None
        self._confColHi = None
        self._dataRepo = data_repo(repoPath)
        self._featureType = featureType

    def writeExpression(self, rnaQuantificationId, quantfile):
        """
        Reads the quantification results file and adds entries to the
        specified database.
        """
        isNormalized = self._isNormalized
        units = self._units
        quantfile.readline()  # skip header
        for expression in quantfile:
            fields = expression.strip().split("\t")
            expressionLevel = fields[self._expressionLevelCol]
            expressionId = fields[self._idCol]
            name = fields[self._nameCol]
            rawCount = 0.0
            if self._countCol is not None:
                rawCount = fields[self._countCol]
            confidenceLow = 0.0
            confidenceHi = 0.0
            score = 0.0
            if (self._confColLow is not None and self._confColHi is not None):
                confidenceLow = float(fields[self._confColLow])
                confidenceHi = float(fields[self._confColHi])
                score = (confidenceLow + confidenceHi)/2

            featureName = fields[self._featureCol]
            dataset = self._dataRepo.getDatasets()[0]
            featureSet = dataset.getFeatureSets()[0]
            featureId = ""
            for feature in featureSet.getFeatures(name=featureName):
                featureId = feature[0].id
            datafields = (expressionId, rnaQuantificationId, name, featureId,
                          expressionLevel, isNormalized,
                          rawCount, score, units, confidenceLow, confidenceHi)
            self._db.addExpression(datafields)
        self._db.batchAddExpression()


class CufflinksWriter(AbstractWriter):
    """
    Class to parse and write expression data from an input file generated by
    Cufflinks.

    cufflinks header:
        tracking_id    class_code    nearest_ref_id    gene_id
        gene_short_name    tss_id    locus    length    coverage    FPKM
        FPKM_conf_lo    FPKM_conf_hi    FPKM_status
    """
    def __init__(self, rnaDB, repoPath, featureType):
        super(CufflinksWriter, self).__init__(rnaDB, repoPath, featureType)
        self._isNormalized = True
        self._units = 1  # FPKM
        self._expressionLevelCol = 9
        self._idCol = 0
        self._nameCol = 4
        self._featureCol = 3
        self._confColLow = 10
        self._confColHi = 11


class RsemWriter(AbstractWriter):
    """
    Class to parse and write expression data from an input file generated by
    RSEM.

    RSEM header (gene quantification):
    gene_id    transcript_id(s)    length    effective_length    expected_count
    TPM    FPKM    posterior_mean_count
    posterior_standard_deviation_of_count    pme_TPM    pme_FPKM
    TPM_ci_lower_bound    TPM_ci_upper_bound    FPKM_ci_lower_bound
    FPKM_ci_upper_bound

    RSEM header (transcript quantification):
    transcript_id    gene_id    length    effective_length    expected_count
    TPM    FPKM    IsoPct    posterior_mean_count
    posterior_standard_deviation_of_count    pme_TPM    pme_FPKM
    IsoPct_from_pme_TPM    TPM_ci_lower_bound    TPM_ci_upper_bound
    FPKM_ci_lower_bound    FPKM_ci_upper_bound
    """
    def __init__(self, rnaDB, repoPath, featureType):
        super(RsemWriter, self).__init__(rnaDB, repoPath, featureType)
        self._isNormalized = True
        self._units = 2  # TPM
        self._expressionLevelCol = 5
        self._idCol = 0
        self._nameCol = self._idCol
        self._countCol = 4
        if self._featureType is "transcript":
            self._featureCol = 1
            self._confColLow = 13
            self._confColHi = 14
        else:
            self._featureCol = 0
            self._confColLow = 11
            self._confColHi = 12


class KallistoWriter(AbstractWriter):
    """
    Class to parse and write expression data from an input file generated by
    RSEM.

    kallisto header:
        target_id    length    eff_length    est_counts    tpm
    """
    def __init__(self, rnaDB):
        super(KallistoWriter, self).__init__(rnaDB)
        self._isNormalized = True
        self._units = 2  # TPM
        self._expressionLevelCol = 4
        self._idCol = 0
        self._nameCol = 0
        self._featureCol = 0
        self._countCol = 3


def writeRnaseqTable(rnaDB, analysisIds, description, annotationId,
                     readGroupId="", programs=""):
    if readGroupId is None:
        readGroupId = ""
    for analysisId in analysisIds:
        datafields = (analysisId, annotationId, description, analysisId,
                      readGroupId, programs)
        rnaDB.addRNAQuantification(datafields)
    rnaDB.batchAddRNAQuantification()


def writeExpressionTable(writer, data):
    for rnaQuantId, quantfile in data:
        writer.writeExpression(rnaQuantId, quantfile)


def rnaseq2ga(dataFolder, sqlFilename, repoPath="ga4gh-example-data/repo.db",
              controlFile="rna_control_file.tsv", featureType="gene"):
    """
    Reads RNA Quantification data in one of several formats and stores the data
    in a sqlite database for use by the GA4GH reference server.

    Quantifications are specified in a tab delimited control file with columns:
    rna_quant_name    filename        type    feature_set_name
    read_group_set_name     description    programs

    Supports the following quantification output type:
    Cufflinks, kallisto, RSEM
    """
    controlFilePath = os.path.join(dataFolder, controlFile)
    dataRepo = data_repo(repoPath)
    rnaDB = RNASqliteStore(sqlFilename, repoPath)
    with open(controlFilePath, "r") as rnaDatasetsFile:
        rnaDatasetsFile.readline()  # skip header
        for line in rnaDatasetsFile:
            fields = line.split("\t")
            rnaType = fields[2]
            dataset = dataRepo.getDatasets()[0]
            annotationId = fields[3].strip()
            featureSet = dataset.getFeatureSetByName(annotationId)
            if rnaType == "cufflinks":
                writer = CufflinksWriter(rnaDB, repoPath)
            elif rnaType == "kallisto":
                writer = KallistoWriter(rnaDB, repoPath)
            elif rnaType == "rsem":
                writer = RsemWriter(rnaDB, repoPath, featureType)
            else:
                raise Exception("Unknown RNA file type: {}".format(rnaType))
            rnaQuantId = fields[0].strip()
            quantFilename = os.path.join(dataFolder, fields[1].strip())
            readGroupSetName = fields[4].strip()
            readGroupSet = dataset.getReadGroupSetByName(readGroupSetName)
            readGroupIds = ",".join(
                [x.getId() for x in readGroupSet.getReadGroups()])
            description = fields[5].strip()
            programs = fields[6].strip()
            writeRnaseqTable(rnaDB, [rnaQuantId], description,
                             featureSet.getId(),
                             readGroupId=readGroupIds, programs=programs)
            quantFile = open(quantFilename, "r")
            writeExpressionTable(writer, [(rnaQuantId, quantFile)])


def parseArgs():
    parser = argparse.ArgumentParser(
        description="Script to generate SQLite database corresponding to "
        "input RNA Quantification experiment files.")
    parser.add_argument(
        "--outputFile", "-o", default="rnaseq.db",
        help="The file to output the server-ready database to.")
    parser.add_argument(
        "--inputDir", "-i",
        help="Path to input directory containing RNA quant files.",
        default='.')
    parser.add_argument(
        "--controlFile", "-c",
        help="Name of control file (.tsv format) in the inputDir",
        default="rna_control_file.tsv")
    parser.add_argument('--verbose', '-v', action='count', default=0)
    args = parser.parse_args()
    return args


@utils.Timed()
def main():
    args = parseArgs()
    rnaseq2ga(args.inputDir, args.outputFile, args.controlFile)


if __name__ == '__main__':
    main()
