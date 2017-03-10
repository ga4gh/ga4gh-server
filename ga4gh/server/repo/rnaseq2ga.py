from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sqlite3
import csv

import ga4gh.server.exceptions as exceptions


SUPPORTED_RNA_INPUT_FORMATS = ["cufflinks", "kallisto", "rsem"]


class RnaSqliteStore(object):
    """
    Defines a sqlite store for RNA data as well as methods for loading the
    tables.
    """
    def __init__(self, sqliteFileName):
        self._dbConn = sqlite3.connect(sqliteFileName)
        self._cursor = self._dbConn.cursor()
        self._batchSize = 2000
        self._rnaValueList = []
        self._expressionValueList = []

    def createTables(self):
        # annotationIds is a comma separated list
        self._cursor.execute('''CREATE TABLE RnaQuantification (
                       id TEXT NOT NULL PRIMARY KEY,
                       feature_set_ids TEXT,
                       description TEXT,
                       name TEXT,
                       read_group_ids TEXT,
                       programs TEXT,
                       biosample_id TEXT)''')
        self._cursor.execute('''CREATE TABLE Expression (
                       id INTEGER,
                       rna_quantification_id TEXT,
                       name TEXT,
                       expression REAL,
                       is_normalized BOOLEAN,
                       raw_read_count REAL,
                       score REAL,
                       units INTEGER,
                       conf_low REAL,
                       conf_hi REAL,
                       PRIMARY KEY (id, rna_quantification_Id))''')
        self._dbConn.commit()

    def addRNAQuantification(self, datafields):
        """
        Adds an RNAQuantification to the db.  Datafields is a tuple in the
        order:
        id, feature_set_ids, description, name,
        read_group_ids, programs, biosample_id
        """
        self._rnaValueList.append(datafields)
        if len(self._rnaValueList) >= self._batchSize:
            self.batchaddRNAQuantification()

    def batchaddRNAQuantification(self):
        if len(self._rnaValueList) > 0:
            sql = "INSERT INTO RnaQuantification VALUES (?,?,?,?,?,?,?)"
            self._cursor.executemany(sql, self._rnaValueList)
            self._dbConn.commit()
            self._rnaValueList = []

    def addExpression(self, datafields):
        """
        Adds an Expression to the db.  Datafields is a tuple in the order:
        id, rna_quantification_id, name, expression,
        is_normalized, raw_read_count, score, units, conf_low, conf_hi
        """
        self._expressionValueList.append(datafields)
        if len(self._expressionValueList) >= self._batchSize:
            self.batchAddExpression()

    def batchAddExpression(self):
        if len(self._expressionValueList) > 0:
            sql = "INSERT INTO Expression VALUES (?,?,?,?,?,?,?,?,?,?)"
            self._cursor.executemany(sql, self._expressionValueList)
            self._dbConn.commit()
            self._expressionValueList = []

    def createIndices(self):
        """
        Index columns that are queried. The expression index can
        take a long time.
        """

        sql = '''CREATE INDEX name_index
                 ON Expression (name)'''
        self._cursor.execute(sql)
        self._dbConn.commit()

        sql = '''CREATE INDEX expression_index
                 ON Expression (expression)'''
        self._cursor.execute(sql)
        self._dbConn.commit()


class AbstractWriter(object):
    """
    Base class to use for the rna quantification writers
    """
    def __init__(self, rnaDB, featureType="gene", dataset=None):
        self._db = rnaDB
        self._isNormalized = None
        self._units = 0  # EXPRESSION_UNIT_UNSPECIFIED
        self._expressionLevelCol = None
        self._idCol = None
        self._nameCol = None
        self._countCol = None
        self._confColLow = None
        self._confColHi = None
        self._dataRepo = None
        self._dataset = dataset
        self._featureType = featureType

    def setUnits(self, units):
        if units == "fpkm":
            self._units = 1
        elif units == "tpm":
            self._units = 2

    def setColNum(self, header, name, defaultNum=None):
        colNum = defaultNum
        try:
            colNum = header.index(name)
        except:
            if defaultNum is None:
                raise exceptions.RepoManagerException(
                       "Missing {} column in expression table.".format(name))
        return colNum

    def writeExpression(self, rnaQuantificationId, quantfilename):
        """
        Reads the quantification results file and adds entries to the
        specified database.
        """
        isNormalized = self._isNormalized
        units = self._units
        with open(quantfilename, "r") as quantFile:
            quantificationReader = csv.reader(quantFile, delimiter=b"\t")
            header = next(quantificationReader)
            expressionLevelColNum = self.setColNum(
                                        header, self._expressionLevelCol)
            nameColNum = self.setColNum(header, self._nameCol)
            countColNum = self.setColNum(header, self._countCol, -1)
            confColLowNum = self.setColNum(header, self._confColLow, -1)
            confColHiNum = self.setColNum(header, self._confColHi, -1)
            expressionId = 0
            for expression in quantificationReader:
                expressionLevel = expression[expressionLevelColNum]
                name = expression[nameColNum]
                rawCount = 0.0
                if countColNum != -1:
                    rawCount = expression[countColNum]
                confidenceLow = 0.0
                confidenceHi = 0.0
                score = 0.0
                if confColLowNum != -1 and confColHiNum != -1:
                    confidenceLow = float(expression[confColLowNum])
                    confidenceHi = float(expression[confColHiNum])
                    score = (confidenceLow + confidenceHi)/2

                datafields = (expressionId, rnaQuantificationId, name,
                              expressionLevel, isNormalized, rawCount, score,
                              units, confidenceLow, confidenceHi)
                self._db.addExpression(datafields)
                expressionId += 1
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
    def __init__(self, rnaDB, featureType, units="fpkm", dataset=None):
        super(CufflinksWriter, self).__init__(
            rnaDB, featureType, dataset=dataset)
        self._isNormalized = True
        self._expressionLevelCol = "FPKM"
        self._idCol = "tracking_id"
        self._nameCol = "tracking_id"
        self._confColLow = "FPKM_conf_lo"
        self._confColHi = "FPKM_conf_hi"
        self.setUnits(units)


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
    def __init__(self, rnaDB, featureType, units="tpm", dataset=None):
        super(RsemWriter, self).__init__(
            rnaDB, featureType=featureType, dataset=dataset)
        self._isNormalized = True
        self._expressionLevelCol = "TPM"
        self._confColLow = "TPM_ci_lower_bound"
        self._confColHi = "TPM_ci_upper_bound"
        self._countCol = "expected_count"
        if self._featureType is "transcript":
            self._idCol = "transcript_id"
        else:
            self._idCol = "gene_id"
        self._nameCol = self._idCol
        self.setUnits(units)


class KallistoWriter(AbstractWriter):
    """
    Class to parse and write expression data from an input file generated by
    RSEM.

    kallisto header:
        target_id    length    eff_length    est_counts    tpm
    """
    def __init__(self, rnaDB, featureType, units="tpm", dataset=None):
        super(KallistoWriter, self).__init__(
            rnaDB, featureType, dataset=dataset)
        self._isNormalized = True
        self._expressionLevelCol = "tpm"
        self._idCol = "target_id"
        self._nameCol = "target_id"
        self._countCol = "est_counts"
        self.setUnits(units)


def writeRnaseqTable(rnaDB, analysisIds, description, annotationId,
                     readGroupId="", programs="", biosampleId=""):
    if readGroupId is None:
        readGroupId = ""
    for analysisId in analysisIds:
        datafields = (analysisId, annotationId, description, analysisId,
                      readGroupId, programs, biosampleId)
        rnaDB.addRNAQuantification(datafields)
    rnaDB.batchaddRNAQuantification()


def writeExpressionTable(writer, data):
    for rnaQuantId, quantFilename in data:
        writer.writeExpression(rnaQuantId, quantFilename)


def rnaseq2ga(quantificationFilename, sqlFilename, localName, rnaType,
              dataset=None, featureType="gene",
              description="", programs="", featureSetNames="",
              readGroupSetNames="", biosampleId=""):
    """
    Reads RNA Quantification data in one of several formats and stores the data
    in a sqlite database for use by the GA4GH reference server.

    Supports the following quantification output types:
    Cufflinks, kallisto, RSEM.
    """
    readGroupSetName = ""
    if readGroupSetNames:
        readGroupSetName = readGroupSetNames.strip().split(",")[0]
    featureSetIds = ""
    readGroupIds = ""
    if dataset:
        featureSetIdList = []
        if featureSetNames:
            for annotationName in featureSetNames.split(","):
                featureSet = dataset.getFeatureSetByName(annotationName)
                featureSetIdList.append(featureSet.getId())
            featureSetIds = ",".join(featureSetIdList)
        # TODO: multiple readGroupSets
        if readGroupSetName:
            readGroupSet = dataset.getReadGroupSetByName(readGroupSetName)
            readGroupIds = ",".join(
                [x.getId() for x in readGroupSet.getReadGroups()])
    if rnaType not in SUPPORTED_RNA_INPUT_FORMATS:
        raise exceptions.UnsupportedFormatException(rnaType)
    rnaDB = RnaSqliteStore(sqlFilename)
    if rnaType == "cufflinks":
        writer = CufflinksWriter(rnaDB, featureType, dataset=dataset)
    elif rnaType == "kallisto":
        writer = KallistoWriter(rnaDB, featureType, dataset=dataset)
    elif rnaType == "rsem":
        writer = RsemWriter(rnaDB, featureType, dataset=dataset)
    writeRnaseqTable(rnaDB, [localName], description, featureSetIds,
                     readGroupId=readGroupIds, programs=programs,
                     biosampleId=biosampleId)
    writeExpressionTable(writer, [(localName, quantificationFilename)])
