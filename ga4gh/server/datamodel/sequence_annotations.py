"""
Module responsible for translating sequence annotation data
into GA4GH native objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import random

import ga4gh.server.datamodel as datamodel
import ga4gh.server.sqlite_backend as sqlite_backend
import ga4gh.server.exceptions as exceptions

import ga4gh.schemas.pb as pb
import ga4gh.schemas.protocol as protocol

# Note to self: There's the Feature ID as understood in a GFF3 file,
# the Feature ID that is its server-assigned compoundId, and the
# ID of the feature's row in the DB FEATURE table.
# I need to be careful about which one is which.

"""
For this implementation, `featureSetId` is required, while `parentId`
is optional, and filters the features within the requested `featureSetId`
by their parent.

Only return features on the reference with this name. Genomic positions
are non-negative integers less than reference length.
Requests spanning the join of circular genomes are represented as two
requests one on each side of the join (position 0) end is also required
If specified, this query matches only annotations which match one of the
provided feature types.
For now do not use the features array in search

GFF3 data is represented by rows in a single table, named FEATURE.
The columns of the FEATURE table correspond to the columns of a GFF3,
with three additional columns prepended representing the ID
of this feature, the ID of its parent (if any), and a whitespace
separated array of its child IDs.

_featureColumns pairs represent the ordered (column_name, column_type).
"""
_featureColumns = [
    ('id', 'TEXT'),  # a synthetic principal key generated on ETL
    ('parent_id', 'TEXT'),
    ('child_ids', 'TEXT'),
    ('reference_name', 'TEXT'),
    ('source', 'TEXT'),
    ('type', 'TEXT'),  # corresponds to featureType, an ontology term
    ('start', 'INT'),
    ('end', 'INT'),
    ('score', 'REAL'),
    ('strand', 'TEXT'),  # limited to one of '+'/'-' or none
    ('name', 'TEXT'),  # the "ID" as found in GFF3, or '' if none
    ('gene_name', 'TEXT'),  # as found in GFF3 attributes
    ('transcript_name', 'TEXT'),  # as found in GFF3 attributes
    ('attributes', 'TEXT')]  # JSON encoding of attributes dict


class Gff3DbBackend(sqlite_backend.SqliteBackedDataSource):
    """
    Notes about the current implementation:
    For this implementation, `featureSetId` is required, while `parentId`
    is optional, and filters the features within the requested `featureSetId`
    by their parent.

    Genomic positions are non-negative integers less than reference length.
    Requests spanning the join of circular genomes are represented as two
    requests one on each side of the join (position 0)
    """

    def __init__(self, dbFile):
        super(Gff3DbBackend, self).__init__(dbFile)
        self.featureColumnNames = [f[0] for f in _featureColumns]
        self.featureColumnTypes = [f[1] for f in _featureColumns]

    def featuresQuery(self, **kwargs):
        """
        Converts a dictionary of keyword arguments into a tuple
        of SQL select statements and the list of SQL arguments
        """
        # TODO: Optimize by refactoring out string concatenation
        sql = ""
        sql_rows = "SELECT * FROM FEATURE WHERE id > 1 "
        sql_args = ()
        if 'name' in kwargs and kwargs['name']:
            sql += "AND name = ? "
            sql_args += (kwargs.get('name'),)
        if 'geneSymbol' in kwargs and kwargs['geneSymbol']:
            sql += "AND gene_name = ? "
            sql_args += (kwargs.get('geneSymbol'),)
        if 'start' in kwargs and kwargs['start'] is not None:
            sql += "AND end > ? "
            sql_args += (kwargs.get('start'),)
        if 'end' in kwargs and kwargs['end'] is not None:
            sql += "AND start < ? "
            sql_args += (kwargs.get('end'),)
        if 'referenceName' in kwargs and kwargs['referenceName']:
            sql += "AND reference_name = ?"
            sql_args += (kwargs.get('referenceName'),)
        if 'parentId' in kwargs and kwargs['parentId']:
            sql += "AND parent_id = ? "
            sql_args += (kwargs['parentId'],)
        if kwargs.get('featureTypes') is not None \
                and len(kwargs['featureTypes']) > 0:
            sql += "AND type IN ("
            sql += ", ".join(["?", ] * len(kwargs.get('featureTypes')))
            sql += ") "
            sql_args += tuple(kwargs.get('featureTypes'))
        sql_rows += sql
        sql_rows += " ORDER BY reference_name, start, end ASC "
        return sql_rows, sql_args

    def searchFeaturesInDb(
            self, startIndex=0, maxResults=None,
            referenceName=None, start=None, end=None,
            parentId=None, featureTypes=None,
            name=None, geneSymbol=None):
        """
        Perform a full features query in database.

        :param startIndex: int representing first record to return
        :param maxResults: int representing number of records to return
        :param referenceName: string representing reference name, ex 'chr1'
        :param start: int position on reference to start search
        :param end: int position on reference to end search >= start
        :param parentId: string restrict search by id of parent node.
        :param name: match features by name
        :param geneSymbol: match features by gene symbol
        :return an array of dictionaries, representing the returned data.
        """
        # TODO: Refactor out common bits of this and the above count query.
        sql, sql_args = self.featuresQuery(
            startIndex=startIndex, maxResults=maxResults,
            referenceName=referenceName, start=start, end=end,
            parentId=parentId, featureTypes=featureTypes,
            name=name, geneSymbol=geneSymbol)
        sql += sqlite_backend.limitsSql(startIndex, maxResults)
        query = self._dbconn.execute(sql, sql_args)
        return sqlite_backend.sqliteRowsToDicts(query.fetchall())

    def getFeatureById(self, featureId):
        """
        Fetch feature by featureID.

        :param featureId: the FeatureID as found in GFF3 records
        :return: dictionary representing a feature object,
            or None if no match is found.
        """
        sql = "SELECT * FROM FEATURE WHERE id = ?"
        query = self._dbconn.execute(sql, (featureId,))
        ret = query.fetchone()
        if ret is None:
            return None
        return sqlite_backend.sqliteRowToDict(ret)


class AbstractFeatureSet(datamodel.DatamodelObject):
    """
    A set of sequence features annotations
    """
    compoundIdClass = datamodel.FeatureSetCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractFeatureSet, self).__init__(parentContainer, localId)
        self._name = localId
        self._sourceUri = ""
        self._referenceSet = None

    def getReferenceSet(self):
        """
        Returns the reference set associated with this FeatureSet.
        """
        return self._referenceSet

    def setReferenceSet(self, referenceSet):
        """
        Sets the reference set associated with this FeatureSet to the
        specified value.
        """
        self._referenceSet = referenceSet

    def toProtocolElement(self):
        """
        Returns the representation of this FeatureSet as the corresponding
        ProtocolElement.
        """
        gaFeatureSet = protocol.FeatureSet()
        gaFeatureSet.id = self.getId()
        gaFeatureSet.dataset_id = self.getParentContainer().getId()
        gaFeatureSet.reference_set_id = pb.string(self._referenceSet.getId())
        gaFeatureSet.name = self._name
        gaFeatureSet.source_uri = self._sourceUri
        attributes = self.getAttributes()
        for key in attributes:
            gaFeatureSet.attributes.attr[key] \
                .values.extend(protocol.encodeValue(attributes[key]))
        return gaFeatureSet

    def getCompoundIdForFeatureId(self, featureId):
        """
        Returns server-style compound ID for an internal featureId.

        :param long featureId: id of feature in database
        :return: string representing ID for the specified GA4GH protocol
            Feature object in this FeatureSet.
        """
        if featureId is not None and featureId != "":
            compoundId = datamodel.FeatureCompoundId(
                self.getCompoundId(), str(featureId))
        else:
            compoundId = ""
        return str(compoundId)


class SimulatedFeatureSet(AbstractFeatureSet):
    """
    Simulated data backend for FeatureSet, used for internal testing.
    """
    def __init__(self, parentContainer, localId, randomSeed=1):
        self._randomSeed = randomSeed
        super(SimulatedFeatureSet, self).__init__(parentContainer, localId)

    def _getRandomfeatureType(self, randomNumberGenerator):
        ontologyTuples = [
            ("gene", "SO:0000704"),
            ("exon", "SO:0000147")]
        term = protocol.OntologyTerm()
        ontologyTuple = randomNumberGenerator.choice(ontologyTuples)
        term.term, term.term_id = ontologyTuple[0], ontologyTuple[1]
        return term

    def _generateSimulatedFeature(self, randomNumberGenerator):
        feature = protocol.Feature()
        feature.feature_set_id = self.getId()
        feature.start = randomNumberGenerator.randint(1000, 2000)
        feature.end = feature.start + randomNumberGenerator.randint(1, 100)
        feature.feature_type.CopyFrom(self._getRandomfeatureType(
            randomNumberGenerator))
        references = ["chr1", "chr2", "chrX"]
        feature.reference_name = randomNumberGenerator.choice(references)
        strands = [protocol.POS_STRAND, protocol.NEG_STRAND]
        feature.strand = randomNumberGenerator.choice(strands)
        attributes = {
            "gene_name": "Frances",
            "gene_type": "mRNA",
            "gene_status": "UNKNOWN"}
        for key, value in attributes.items():
            feature.attributes.attr[key].values.add().string_value = value
        return feature

    def getFeature(self, compoundId):
        """
        Fetches a simulated feature by ID.

        :param compoundId: any non-null string
        :return: A simulated feature with id set to the same value as the
            passed-in compoundId.
        ":raises: exceptions.ObjectWithIdNotFoundException if None is passed
            in for the compoundId.
        """
        if compoundId is None:
            raise exceptions.ObjectWithIdNotFoundException(compoundId)
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(self._randomSeed)
        feature = self._generateSimulatedFeature(randomNumberGenerator)
        feature.id = str(compoundId)
        feature.parent_id = ""  # TODO: Test with nonempty parentIDs?
        return feature

    def getFeatures(self, referenceName=None, start=None, end=None,
                    startIndex=None, maxResults=None,
                    featureTypes=None, parentId=None,
                    name=None, geneSymbol=None, numFeatures=10):
        """
        Returns a set number of simulated features.

        :param referenceName: name of reference to "search" on
        :param start: start coordinate of query
        :param end: end coordinate of query
        :param startIndex: None or int
        :param maxResults: None or int
        :param featureTypes: optional list of ontology terms to limit query
        :param parentId: optional parentId to limit query.
        :param name: the name of the feature
        :param geneSymbol: the symbol for the gene the features are on
        :param numFeatures: number of features to generate in the return.
            10 is a reasonable (if arbitrary) default.
        :return: Yields feature list
        """
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(self._randomSeed)
        for featureId in range(numFeatures):
            gaFeature = self._generateSimulatedFeature(randomNumberGenerator)
            gaFeature.id = self.getCompoundIdForFeatureId(featureId)
            match = (
                gaFeature.start < end and
                gaFeature.end > start and
                gaFeature.reference_name == referenceName and (
                    featureTypes is None or len(featureTypes) == 0 or
                    gaFeature.feature_type in featureTypes))
            if match:
                gaFeature.parent_id = ""  # TODO: Test nonempty parentIDs?
                yield gaFeature


class Gff3DbFeatureSet(AbstractFeatureSet):
    """
    Stub class to directly read sequence annotation features from GFF3 files.
    Tests basic access, not to be used in production.
    """
    def __init__(self, parentContainer, localId):
        super(Gff3DbFeatureSet, self).__init__(parentContainer, localId)
        self._ontology = None
        self._dbFilePath = None
        self._db = None

    def setOntology(self, ontology):
        """
        Sets the Ontology instance used by this FeatureSet to the
        specified value.
        """
        self._ontology = ontology

    def getOntology(self):
        """
        Returns the ontology term map used to translate ontology term names
        to IDs.
        """
        return self._ontology

    def populateFromFile(self, dataUrl):
        """
        Populates the instance variables of this FeatureSet from the specified
        data URL.
        """
        self._dbFilePath = dataUrl
        self._db = Gff3DbBackend(self._dbFilePath)

    def populateFromRow(self, featureSetRecord):
        """
        Populates the instance variables of this FeatureSet from the specified
        DB row.
        """
        self._dbFilePath = featureSetRecord.dataurl
        self.setAttributesJson(featureSetRecord.attributes)
        self._db = Gff3DbBackend(self._dbFilePath)

    def getDataUrl(self):
        """
        Returns the URL providing the data source for this FeatureSet.
        """
        return self._dbFilePath

    def getFeature(self, compoundId):
        """
        Returns a protocol.Feature object corresponding to a compoundId
        :param compoundId: a datamodel.FeatureCompoundId object
        :return: a Feature object.
        :raises: exceptions.ObjectWithIdNotFoundException if invalid
            compoundId is provided.
        """
        featureId = long(compoundId.featureId)
        with self._db as dataSource:
            featureReturned = dataSource.getFeatureById(featureId)

        if featureReturned is None:
            raise exceptions.ObjectWithIdNotFoundException(compoundId)
        else:
            gaFeature = self._gaFeatureForFeatureDbRecord(featureReturned)
            return gaFeature

    def _gaFeatureForFeatureDbRecord(self, feature):
        """
        :param feature: The DB Row representing a feature
        :return: the corresponding GA4GH protocol.Feature object
        """
        gaFeature = protocol.Feature()
        gaFeature.id = self.getCompoundIdForFeatureId(feature['id'])
        if feature.get('parent_id'):
            gaFeature.parent_id = self.getCompoundIdForFeatureId(
                    feature['parent_id'])
        else:
            gaFeature.parent_id = ""
        gaFeature.feature_set_id = self.getId()
        gaFeature.reference_name = pb.string(feature.get('reference_name'))
        gaFeature.start = pb.int(feature.get('start'))
        gaFeature.end = pb.int(feature.get('end'))
        gaFeature.name = pb.string(feature.get('name'))
        if feature.get('strand', '') == '-':
            gaFeature.strand = protocol.NEG_STRAND
        else:
            # default to positive strand
            gaFeature.strand = protocol.POS_STRAND
        gaFeature.child_ids.extend(map(
                self.getCompoundIdForFeatureId,
                json.loads(feature['child_ids'])))
        gaFeature.feature_type.CopyFrom(
            self._ontology.getGaTermByName(feature['type']))
        attributes = json.loads(feature['attributes'])
        # TODO: Identify which values are ExternalIdentifiers and OntologyTerms
        for key in attributes:
            for v in attributes[key]:
                gaFeature.attributes.attr[key].values.add().string_value = v
        if 'gene_name' in attributes and len(attributes['gene_name']) > 0:
            gaFeature.gene_symbol = pb.string(attributes['gene_name'][0])
        return gaFeature

    def getFeatures(self, referenceName=None, start=None, end=None,
                    startIndex=None, maxResults=None,
                    featureTypes=None, parentId=None,
                    name=None, geneSymbol=None):
        """
        method passed to runSearchRequest to fulfill the request
        :param str referenceName: name of reference (ex: "chr1")
        :param start: castable to int, start position on reference
        :param end: castable to int, end position on reference
        :param startIndex: none or castable to int
        :param maxResults: none or castable to int
        :param featureTypes: array of str
        :param parentId: none or featureID of parent
        :param name: the name of the feature
        :param geneSymbol: the symbol for the gene the features are on
        :return: yields a protocol.Feature at a time
        """
        with self._db as dataSource:
            features = dataSource.searchFeaturesInDb(
                startIndex, maxResults,
                referenceName=referenceName,
                start=start, end=end,
                parentId=parentId, featureTypes=featureTypes,
                name=name, geneSymbol=geneSymbol)
            for feature in features:
                gaFeature = self._gaFeatureForFeatureDbRecord(feature)
                yield gaFeature
