"""
Server classes for the GA4GH reference implementation.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import random
import datetime
import future
from future.standard_library import hooks
with hooks():
    import http.server

import wormtable as wt

import ga4gh
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
    # This must be a bytes literal for integration with wormtable.
    GENOTYPE_LIKELIHOOD_NAME = b"GL"

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
        self._infoCols = []
        self._firstSamplePosition = -1
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
            # We assume the .GT is the first column for each sample
            elif colName.endswith(".GT"):
                s = colName.split(".")[0]
                self._sampleCols[s] = [(c, "GT")]
                if self._firstSamplePosition == -1:
                    self._firstSamplePosition = c.get_position()
            else:
                # This must be a sample specific column
                s = colName.split(".")
                self._sampleCols[s[0]].append((c, s[1]))

    def convertVariant(self, row, callSetIds):
        """
        Converts the specified wormtable row into a GAVariant object including
        the specified set of callSetIds.
        """
        v = protocol.GAVariant()
        v.variantSetId = self._variantSetId
        v.referenceName = row[self.CHROM_COL]
        v.start = row[self.POS_COL]
        names = row[self.ID_COL]
        if names is not None:
            v.names = names.split(";")
        v.referenceBases = row[self.REF_COL]
        v.end = v.start + len(v.referenceBases)
        v.id = "{0}.{1}.{2}".format(v.variantSetId, v.referenceName, v.start)
        alt = row[self.ALT_COL].split(",")
        v.alternateBases = alt
        v.info = []
        for infoField, col in self._infoCols:
            pos = col.get_position()
            if row[pos] is not None:
                keyValue = protocol.GAKeyValue(infoField, row[pos])
                v.info.append(keyValue)
        v.calls = []
        # All of the remaining values in the row correspond to Calls.
        j = self._firstSamplePosition
        for csid in callSetIds:
            call = protocol.GACall()
            # Why do we have both of these? Wouldn't the ID be sufficient
            # to call back to searchCallSets to get more info?
            call.callSetId = csid
            call.callSetName = csid
            genotype = row[j]
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
            for col, infoName in self._sampleCols[csid][1:]:
                j += 1
                if infoName == self.GENOTYPE_LIKELIHOOD_NAME:
                    call.genotypeLikelihood = row[j]
                else:
                    kv = protocol.GAKeyValue(infoName, str(row[j]))
                    call.info.append(kv)
            j += 1
            v.calls.append(call)
        return v

    def searchVariants(self, request):
        """
        Serves the specified GASearchVariantsRequest and returns a
        GASearchVariantsResponse. If the number of variants to be returned is
        greater than request.maxResults then the nextPageToken is set to a
        non-null value. Subsequent request objects should provide this value in
        the pageToken attribute to obtain the next page of results.
        """
        response = protocol.GASearchVariantsResponse()
        response.variantSetId = self._variantSetId
        # First we get the set of columns for wormtable to retrieve. We don't
        # want to waste time reading in and decoding data for columns that
        # are not needed for the output when we specify callSetIds
        if len(request.callSetIds) == 0:
            readCols = self._table.columns()
            callSetIds = self._sampleCols.keys()
        else:
            readCols = self._table.columns()[:self._firstSamplePosition]
            callSetIds = request.callSetIds
            for csid in request.callSetIds:
                cols = [c for c, name in self._sampleCols[csid]]
                readCols.extend(cols)
        v = []
        # TODO encoding issues here? Wormtable only deals with bytes.
        chrom = request.referenceName.encode()
        if request.variantName is None:
            start = (chrom, request.start)
            end = (chrom, request.end)
            if request.pageToken is not None:
                # TODO we must sanitise the input here! This should be an int.
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
            while r is not None and len(v) < request.maxResults:
                v.append(self.convertVariant(r, callSetIds))
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
            # TODO are there issues with having several names/IDs for a variant?
            if r is not None:
                # The result must still be within the range and must match
                # the specified name exactly. The cursor is positioned at
                # the first row >= the specified key.
                if request.start <= r[self.POS_COL] < request.end \
                        and r[self.ID_COL] == name:
                    v.append(self.convertVariant(r, callSetIds))
        response.variants = v
        return response


class WormtableBackend(object):
    """
    A backend based on wormtable tables and indexes. This backend is provided
    with a directory within which it can find VCF variant sets converted into
    wormtable format. Each file within the specified directory is assumed to
    be a wormtable, and the variant set ID will be the name of the
    corresponding wormtable directory.
    """
    def __init__(self, dataDir):
        self._dataDir = dataDir
        self._variantSets = {}
        # All files in datadir are assumed to correspond to wormtables.
        for vsid in os.listdir(self._dataDir):
            f = os.path.join(self._dataDir, vsid)
            self._variantSets[vsid] = WormtableDataset(vsid, f)

    def searchVariants(self, request):
        # TODO rearrange the code so we can support searching over multiple
        # variantSets.
        assert len(request.variantSetIds) == 1
        ds = self._variantSets[request.variantSetIds[0]]
        return ds.searchVariants(request)


class VariantSimulator(object):
    """
    A class that simulates Variants that can be served by the GA4GH API.
    """
    def __init__(self, seed=0, numCalls=1, variantDensity=0.5):
        self._randomSeed = seed
        self._numCalls = numCalls
        self._variantDensity = variantDensity
        now = protocol.convertDatetime(datetime.datetime.now())
        self._created = now
        self._updated = now

    def generateVariant(self, variantSetId, referenceName, position, rng):
        """
        Generate a random variant for the specified position using the
        specified random number generator. This generator should be seeded
        with a value that is unique to this position so that the same variant
        will always be produced regardless of the order it is generated in.
        """
        v = protocol.GAVariant()
        v.variantSetId = variantSetId
        # The id is the combination of the position, referenceName and variant
        # set id; this allows us to generate the variant from the position and
        # id.
        v.id = "{0}:{1}:{2}".format(
            v.variantSetId, referenceName, position)
        v.referenceName = referenceName
        v.names = []  # What's a good model to generate these?
        v.created = self._created
        v.updated = self._updated
        v.start = position
        v.end = position + 1  # SNPs only for now
        bases = ["A", "C", "G", "T"]
        ref = rng.choice(bases)
        v.referenceBases = ref
        alt = rng.choice([b for b in bases if b != ref])
        v.alternateBases = [alt]
        v.calls = []
        for j in range(self._numCalls):
            c = protocol.GACall()
            # for now, the genotype is either [0,1], [1,1] or [1,0] with equal
            # probability; probably will want to do something more
            # sophisticated later.
            g = rng.choice([[0, 1], [1, 0], [1, 1]])
            c.genotype = g
            # TODO What is a reasonable model for generating these likelihoods?
            # Are these log-scaled? Spec does not say.
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
        response = protocol.GASearchVariantsResponse()
        rng = random.Random()
        v = []
        j = request.start
        if request.pageToken is not None:
            j = request.pageToken
        while j < request.end and len(v) != request.maxResults:
            rng.seed(self._randomSeed + j)
            if rng.random() < self._variantDensity:
                # TODO fix variant set IDS so we can have multiple
                v.append(self.generateVariant(
                    request.variantSetIds[0], request.referenceName, j, rng))
            j += 1
        if j < request.end - 1:
            response.nextPageToken = j
        response.variants = v
        return response


class ProtocolHandler(object):
    """
    Class that handles the GA4GH protocol messages and responses.
    """
    def __init__(self, backend):
        self._backend = backend

    def searchVariants(self, jsonRequest):
        """
        Handles the specified JSON encoding of a GASearchVariantsRequest.
        and returns the corresponding JSON encoded GASearchVariantsResponse
        (in case of success) or GAException (in case of error).
        """
        request = protocol.GASearchVariantsRequest.fromJSON(jsonRequest)
        # TODO wrap this call in try: except and make a GAException object
        # out of the resulting Exception object. Two classes of exception
        # should be identified: those due to input errors and other expected
        # problems the backend must deal with, and other exceptions which
        # indicate a server error. The former type should all subclass a
        # an exception defined in the ga4gh package.
        resp = self._backend.searchVariants(request)
        s = resp.toJSON()
        return s


class HTTPRequestHandler(http.server.BaseHTTPRequestHandler):
    """
    Handler for the HTTP level of the GA4GH protocol.
    """

    def do_POST(self):
        """
        Handle a single POST request.
        """
        h = self.server.ga4ghProtocolHandler
        # TODO read the path and 404 if not correct
        length = int(self.headers['Content-Length'])
        # TODO is this safe encoding-wise? Do we need to specify an
        # explicit encoding?
        jsonRequest = self.rfile.read(length).decode()
        s = h.searchVariants(jsonRequest).encode()
        self.send_response(200)
        self.send_header("Content-type", "application/json")
        self.send_header("Content-Length", len(s))
        self.end_headers()
        self.wfile.write(s)


class HTTPServer(http.server.HTTPServer):
    """
    Basic HTTP server for the GA4GH protocol.
    """
    def __init__(self, serverAddress, backend):
        # Cannot use super() here because of Python 2 issues.
        http.server.HTTPServer.__init__(
            self, serverAddress, HTTPRequestHandler)
        self.ga4ghProtocolHandler = ProtocolHandler(backend)
