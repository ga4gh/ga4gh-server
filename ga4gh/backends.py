"""
Server classes for the GA4GH reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import datetime

import wormtable as wt
import pysam

import ga4gh.protocol as protocol


class WormtableDataset(object):
    """
    Class representing a single dataset backed by a wormtable directory.
    We assume that VCF data has been converted to wormtable format using
    vcf2wt, which maps VCF data to wormtable columns in the following
    way:
    o The zeroth column (row_id) in a row is the primary key;
    o The next 7 columns correspond to the fixed fields in VCF, and are
      called CHROM, POS, etc.
    o There is then a single column for each INFO field defined in the
      VCF header. The column names are prefixed with "INFO.", and so
      we might have INFO.AF, INFO.THETA etc.
    o There is then a set of columns for each sample, one for each field
      defined in the FORMAT rows of the header. These are prefixed by the
      sample ID and suffixed with the field name, and so we have fields
      like HG00096.GT, NA20828.GL etc.
    The ``wtadmin show`` command is useful to see the columns in a table.

    We assume that indexes on the combination of the CHROM and POS columns
    and the CHROM and ID columns have been built. This is to support
    efficient retrieval of rows based on both coordinates and names
    within a chromosome.
    """
    # Positions of the fixed columns within a row
    CHROM_COL = 1
    POS_COL = 2
    ID_COL = 3
    REF_COL = 4
    ALT_COL = 5
    QUAL_COL = 6
    FILTER_COL = 7
    # These must be bytes literals for integration with wormtable.
    GENOTYPE_LIKELIHOOD_NAME = b"GL"
    GENOTYPE_NAME = b"GT"

    def __init__(self, variantSetId, wtDir):
        """
        Allocates a new WormtableDataset with the specified variantSetId
        based on the specified wormtable directory.
        """
        self._variantSetId = variantSetId
        self._wtDir = wtDir
        self._table = wt.open_table(wtDir)
        self._chromPosIndex = self._table.open_index("CHROM+POS")
        self._chromIdIndex = self._table.open_index("CHROM+ID")
        self._sampleCols = {}
        self._sampleNames = []
        self._infoCols = []
        self._firstSamplePosition = -1
        t = int(os.path.getctime(wtDir) * 1000)
        # ctime is in seconds, and we want milliseconds since the epoch
        self._creationTime = t
        self._updatedTime = t
        cols = self._table.columns()[self.FILTER_COL + 1:]
        # We build lookup tables for the INFO and sample columns so they can
        # be easily found during conversion. For the sample columns we make
        # a dictionary mapping the sample name to a the list of (name, col)
        # tuples for that sample.
        for c in cols:
            colName = c.get_name()
            if colName.startswith("INFO"):
                s = colName.split(".")[1]
                self._infoCols.append((s, c))
            else:
                if self._firstSamplePosition == -1:
                    # This must be a sample specific column
                    self._firstSamplePosition = c.get_position()
                sampleName, infoName = colName.split(".")
                if sampleName not in self._sampleCols:
                    self._sampleCols[sampleName] = []
                    self._sampleNames.append(sampleName)
                self._sampleCols[sampleName].append((infoName, c))

    def convertInfoField(self, v):
        """
        Converts the specified value into an appropriate format for a protocol
        info field. Info fields are arrays of strings.
        """
        if isinstance(v, tuple):
            ret = [str(x) for x in v]
        else:
            ret = [str(v)]
        return ret

    def convertGenotype(self, call, genotype):
        """
        Updates the specified call to reflect the value encoded in the
        specified VCF genotype string.
        """
        if genotype is not None:
            # Is this the correct interpretation of | and /?
            delim = "/"
            if "|" in genotype:
                delim = "|"
                # TODO what is the phaseset value supposed to be?
                call.phaseset = "True"
            try:
                call.genotype = map(int, genotype.split(delim))
            except ValueError:
                # TODO what is the correct interpretation of .|.?
                call.genotype = [-1]

    def convertVariant(self, row, sampleRowPositions):
        """
        Converts the specified wormtable row into a GAVariant object including
        the specified set of callSetIds.
        """
        v = protocol.GAVariant()
        v.created = self._creationTime
        v.updated = self._updatedTime
        v.variantSetId = self._variantSetId
        v.referenceName = row[self.CHROM_COL]
        v.start = row[self.POS_COL]
        names = row[self.ID_COL]
        if names is not None:
            v.names = names.split(";")
        v.referenceBases = row[self.REF_COL]
        v.end = v.start + len(v.referenceBases)
        v.id = "{0}.{1}".format(v.variantSetId, row[0])
        alt = row[self.ALT_COL]
        if alt is not None:
            v.alternateBases = alt.split(",")
        v.info = {}
        for infoField, col in self._infoCols:
            pos = col.get_position()
            if row[pos] is not None:
                v.info[infoField] = self.convertInfoField(row[pos])
        v.calls = []
        # All of the remaining values in the row correspond to Calls.
        for csid, rowPositions in sampleRowPositions.items():
            call = protocol.GACall()
            # Why do we have both of these? Wouldn't the ID be sufficient
            # to call back to searchCallSets to get more info?
            call.callSetId = csid
            call.callSetName = csid
            for (info, col), j in zip(self._sampleCols[csid], rowPositions):
                if info == self.GENOTYPE_LIKELIHOOD_NAME:
                    call.genotypeLikelihood = row[j]
                elif info == self.GENOTYPE_NAME:
                    self.convertGenotype(call, row[j])
                else:
                    if row[j] is not None:
                        # Missing values are not included in the info array
                        call.info[info] = self.convertInfoField(row[j])
            v.calls.append(call)
        return v

    def searchVariants(self, request):
        """
        Serves the specified GASearchVariantsRequest and returns a
        GASearchVariantsResponse. If the number of variants to be returned is
        greater than request.pageSize then the nextPageToken is set to a
        non-null value. Subsequent request objects should provide this value in
        the pageToken attribute to obtain the next page of results.
        """
        # This is temporary until we do this properly based on the output
        # buffer size
        if request.pageSize is None:
            request.pageSize = 100
        response = protocol.GASearchVariantsResponse()
        # First we get the set of columns for wormtable to retrieve. We don't
        # want to waste time reading in and decoding data for columns that
        # are not needed for the output when we specify callSetIds
        if request.callSetIds is None:
            readCols = self._table.columns()
            # These must be in the correct order, so we cannot use the keys of
            # self._sampleCols
            callSetIds = self._sampleNames
        else:
            readCols = self._table.columns()[:self._firstSamplePosition]
            callSetIds = request.callSetIds
            for csid in request.callSetIds:
                cols = [c for name, c in self._sampleCols[csid]]
                readCols.extend(cols)
        # Now we get the row positions for the sample columns
        sampleRowPositions = {}
        j = self._firstSamplePosition
        for csid in callSetIds:
            l = []
            for c in self._sampleCols[csid]:
                l.append(j)
                j += 1
            sampleRowPositions[csid] = l
        v = []
        # TODO encoding issues here? Wormtable only deals with bytes.
        chrom = request.referenceName.encode()
        if request.variantName is None:
            start = (chrom, request.start)
            end = (chrom, request.end)
            if request.pageToken is not None:
                # TODO we must sanitise the input here! This should be an int.
                # We must also check chrom to make sure it's not nasty.
                start = (chrom, request.pageToken)
            # A wormtable cursor on an index returns an iterator over the rows
            # in the table in the order defined by the index, starting with
            # the first two >= to the index key start. Each row returned is a
            # list, where each entry is the decoded value for the
            # corresponding column in readCols.
            cursor = self._chromPosIndex.cursor(readCols, start, end)
            # Normally, we would use a for loop to iterate over the rows in
            # a cursor, but here we use the next() function so that we can
            # end the iteration early when we have generated enough results.
            r = next(cursor, None)
            while r is not None and len(v) < request.pageSize:
                v.append(self.convertVariant(r, sampleRowPositions))
                r = next(cursor, None)
            if r is not None:
                response.nextPageToken = r[self.POS_COL]
        else:
            # TODO more encoding issues.
            name = request.variantName.encode()
            start = (chrom, name)
            cursor = self._chromIdIndex.cursor(readCols, start)
            r = next(cursor, None)
            # for now, we assume that there are 0 or 1 results.
            # TODO are there issues with having several names for a variant?
            if r is not None:
                # The result must still be within the range and must match
                # the specified name exactly. The cursor is positioned at
                # the first row >= the specified key.
                if request.start <= r[self.POS_COL] < request.end \
                        and r[self.ID_COL] == name:
                    v.append(self.convertVariant(r, sampleRowPositions))
        response.variants = v
        return response

    def getMetadata(self):
        """
        Returns a list of GAVariantSetMetadata objects for this variant set.
        """
        def f(infoField, col):
            md = protocol.GAVariantSetMetadata()
            md.key = infoField
            md.value = ""
            md.id = ""
            # What are the encodings here? With the lack of any precise
            # definitions I'm just using wormtable values for now.
            md.type = col.get_type_name()
            md.number = str(col.get_num_elements())
            md.description = col.get_description()
            return md
        ret = []
        for infoField, col in self._infoCols:
            ret.append(f(infoField, col))
        if len(self._sampleCols) > 0:
            # TODO this is pretty nasty, making a list just to take the head.
            sampleName = list(self._sampleCols.keys())[0]
            for infoField, col in self._sampleCols[sampleName]:
                if infoField != self.GENOTYPE_NAME:
                    ret.append(f(infoField, col))
        return ret


class TabixDataset(object):
    """
    Class representing a single dataset backed by a tabix directory.
    """
    def __init__(self, variantSetId, vcfPath):
        self._variantSetId = variantSetId
        self._created = protocol.convertDatetime(datetime.datetime.now())

        # Read in Tabix-index VCF files from vcfPath is a path.
        # self.tb chrom -> vcf dict
        self._tabixMap = {}

        for fn in glob.glob(os.path.join(vcfPath, "*.vcf.gz")):
            tabixfile = pysam.Tabixfile(fn)
            for chrom in tabixfile.contigs:
                if chrom in self._tabixMap:
                    raise Exception("cannot have overlapping VCF files.")
                self._tabixMap[chrom] = tabixfile

    def convertVariant(self, record):
        v = protocol.GAVariant()
        record = record.split('\t')
        position = int(record[1])

        v.id = "{0}:{1}:{2}".format(self._variantSetId, record[0], position)
        v.variantSetId = self._variantSetId
        v.referenceName = record[0]
        v.names = []  # What's a good model to generate these?
        # TODO How should we populate these from VCF?
        v.created = self._created
        v.updated = self._created
        v.start = position
        v.end = position + 1  # TODO support non SNP variaints
        v.referenceBases = record[3]
        v.alternateBases = [record[4].split(",")]
        v.calls = []
        for j in range(9, len(record)):
            c = protocol.GACall()
            c.genotype = record[j].split(":")[0].split("[\/\|]")
            # "(\d+)([\|\/])(\d+)"
            # Need to make up genotypeLikelihood for dev vcf.
            # make something more robust later
            c.genotypeLikelihood = [-100, -100, -100]
            v.calls.append(c)
        return v

    def searchVariants(self, request):
        """
        Serves the specified GASearchVariantsRequest and returns a
        GASearchVariantsResponse. If the number of variants to be returned is
        greater than request.maxResults then the nextPageToken is set to a
        non-null value. Subsequent request objects should provide this value in
        the pageToken attribute to obtain the next page of results.
        """
        # This is temporary until we do this properly based on the output
        # buffer size
        if request.pageSize is None:
            request.pageSize = 100
        response = protocol.GASearchVariantsResponse()
        response.variantSetId = self._variantSetId

        v = []
        j = int(request.start)
        if request.pageToken is not None:
            j = int(request.pageToken)
        if request.end is None:
            request.end = request.start + 1

        cursor = self._tabixMap[request.referenceName]\
            .fetch(str(request.referenceName), j-1, int(request.end))
        # Weird with the -1 but seems required. 0 based i guess?
        while len(v) < request.pageSize:
            try:
                record = cursor.next()
            except StopIteration:
                break
            v.append(self.convertVariant(record))
        try:
            record = cursor.next()
        except StopIteration:
            pass
        else:
            response.nextPageToken = self.convertVariant(record).start
        response.variants = v

        return response

    def getMetadata(self):
        # TODO: Implement this
        ret = []
        return ret


class Backend(object):
    """
    Superclass of GA4GH protocol backends.
    """
    def __init__(self, dataDir, Dataset):
        self._dataDir = dataDir
        self._variantSets = {}
        # All files in datadir are assumed to correspond to Datasets.
        for vsid in os.listdir(self._dataDir):
            f = os.path.join(self._dataDir, vsid)
            self._variantSets[vsid] = Dataset(vsid, f)
        self._variantSetIds = sorted(self._variantSets.keys())

    def searchVariants(self, request):
        assert len(request.variantSetIds) > 0
        variantSetIndex = 0  # x
        if request.pageToken is not None:
            # parse the pageToken and change what will be passed in
            variantSetIndex, pageToken = request.pageToken.split(':')
            variantSetIndex = int(variantSetIndex)
            if pageToken == '':
                pageToken = None  # None is represented as '' in the pageToken
            else:
                pageToken = int(pageToken)
            request.pageToken = pageToken
        ds = self._variantSets[request.variantSetIds[variantSetIndex]]
        response = ds.searchVariants(request)
        # Add the index of the variant set of the next
        # page to the nextPageToken
        if response.nextPageToken is None:
            if variantSetIndex < len(request.variantSetIds) - 1:
                # if not, there are no more results, so leave the token as None
                response.nextPageToken = '{0}:'.format(variantSetIndex + 1)
                # this token will give results from the beginning of the
                # next variantSet
        else:
            response.nextPageToken = "{0}:{1}".format(variantSetIndex,
                                                      response.nextPageToken)
        return response

    def searchVariantSets(self, request):
        """
        Returns a GASearchVariantSets response for the specified
        GASearchVariantSetsRequest object.
        """
        # This is temporary until we do this properly based on the output
        # buffer size
        if request.pageSize is None:
            request.pageSize = 100
        response = protocol.GASearchVariantSetsResponse()
        j = 0
        if request.pageToken is not None:
            j = int(request.pageToken)
        v = []
        while j < len(self._variantSetIds) and len(v) < request.pageSize:
            vs = protocol.GAVariantSet()
            vs.id = self._variantSetIds[j]
            vs.datasetId = "NotImplemented"
            vs.metadata = self._variantSets[vs.id].getMetadata()
            v.append(vs)
            j += 1
        if j < len(self._variantSetIds):
            response.nextPageToken = j
        response.variantSets = v
        return response


class TabixBackend(Backend):
    """
    A class that serves variants indexed by tabix.
    """
    def __init__(self, dataDir):
        Backend.__init__(self, dataDir, TabixDataset)


class WormtableBackend(Backend):
    """
    A backend based on wormtable tables and indexes. This backend is provided
    with a directory within which it can find VCF variant sets converted into
    wormtable format. Each file within the specified directory is assumed to
    be a wormtable, and the variant set ID will be the name of the
    corresponding wormtable directory.
    """
    def __init__(self, dataDir):
        Backend.__init__(self, dataDir, WormtableDataset)
