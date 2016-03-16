"""
Module responsible for translating sequence annotation data
into GA4GH native objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import random

import ga4gh.protocol as protocol
import ga4gh.datamodel as datamodel
import ga4gh.sqliteBackend as sqliteBackend

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
_featureColumns = [('id', 'TEXT'),
                   ('parent_id', 'TEXT'),
                   ('child_ids', 'TEXT'),
                   ('reference_name', 'TEXT'),
                   ('source', 'TEXT'),
                   ('ontology_term', 'TEXT'),
                   ('start', 'INT'),
                   ('end', 'INT'),
                   ('score', 'REAL'),
                   ('strand', 'TEXT'),  # limited to one of '+'/'-' or none.
                   ('attributes', 'TEXT')  # JSON encoding of attributes dict
                   ]


class Gff3DbBackend(sqliteBackend.SqliteBackedDataSource):
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

    def countFeaturesSearchInDb(self,
                                referenceName=None, start=0, end=0,
                                parentId=None, ontologyTerms=[]):
        """
        Same parameters as searchFeaturesInDb,
        except without the pagetoken/size.
        """
        sql = ("SELECT COUNT(*) FROM FEATURE WHERE "
               "reference_name = ? "
               "AND end > ? "  # compare this to query start
               "AND start < ? "  # and this to query end
               )
        sql_args = (referenceName, start, end)
        if parentId is not None:
            sql += "AND parent_id = ? "
            sql_args += (parentId,)
        if ontologyTerms is not None and len(ontologyTerms) > 0:
            sql += "AND ontology_term IN ("
            sql += ", ".join(["?", ] * len(ontologyTerms))
            sql += ") "
            sql_args += tuple(ontologyTerms)
        query = self._dbconn.execute(sql, sql_args)
        return (query.fetchone())[0]

    def searchFeaturesInDb(self, pageToken=0, pageSize=None,
                           referenceName=None, start=0, end=0,
                           parentId=None, ontologyTerms=[]):
        """
        :param pageToken: int representing first record to return
        :param pageSize: int representing number of records to return
        :param referenceName: string representing reference name, ex 'chr1'
        :param start: int position on reference to start search
        :param end: int position on reference to end search >= start
        :param parentId: string restrict search by id of parent node.
        :return an array of dictionaries, representing the returned data.
        """
        sql = ("SELECT * FROM FEATURE WHERE "
               "reference_name = ? "
               "AND end > ? "  # compare this to query start
               "AND start < ? "  # and this to query end
               )
        sql_args = (referenceName, start, end)
        if parentId is not None:
            sql += "AND parent_id = ? "
            sql_args += (parentId,)
        if ontologyTerms is not None and len(ontologyTerms) > 0:
            sql += "AND ontology_term IN ("
            sql += ", ".join(["?", ] * len(ontologyTerms))
            sql += ") "
            sql_args += tuple(ontologyTerms)
        sql += sqliteBackend.limitsSql(pageToken, pageSize)
        query = self._dbconn.execute(sql, sql_args)
        return sqliteBackend.sqliteRows2dicts(query.fetchall())

    def getFeatureById(self, featureId):
        """
        :param featureId: the FeatureID as found in GFF3 records
        :return: dictionary representing a feature object,
            or None if no match is found.
        """
        sql = ("SELECT * FROM FEATURE WHERE id = ?")
        query = self._dbconn.execute(sql, (featureId,))
        return sqliteBackend.sqliteRow2Dict(query.fetchone())


# Problem: FeatureIds from GFF3 often contain colons,
# and the CompoundId mechanism reserves colons for its special separator
# character. Thus the following pair of functions that mangles and
# de-mangles featureId's to allow them to work in compoundIds.
def _decolonize(featureId):
    """
    Replaces colons in input with CompoundId-safe character
    :param featureId: string potentially with colons in it
    :return: same string with colons replaced with safe character
    """
    return featureId.replace(':', ';')


def _recolonize(mangled_featureId):
    """
    Replaces CompoundId-safe character in input with colons
    :param featureId: string potentially with replacement character
    :return: same string with colons in place.
    """
    return mangled_featureId.replace(';', ':')


class AbstractFeatureSet(datamodel.DatamodelObject):
    """
    A set of sequence features annotations
    """
    compoundIdClass = datamodel.FeatureSetCompoundId

    def __init__(
            self, parentContainer, localId, referenceSetId=None):
        super(AbstractFeatureSet, self).__init__(parentContainer, localId)
        self._referenceSetId = referenceSetId
        self._name = localId
        self._sourceUri = ""
        self._attributes = protocol.Attributes()

    def toProtocolElement(self):
        """
        Returns the representation of this FeatureSet as the corresponding
        ProtocolElement.
        """
        gaFeatureSet = protocol.FeatureSet()
        gaFeatureSet.id = self.getId()
        gaFeatureSet.datasetId = self.getParentContainer().getId()
        gaFeatureSet.referenceSetId = self._referenceSetId
        gaFeatureSet.name = self._name
        gaFeatureSet.sourceUri = self._sourceUri
        gaFeatureSet.attributes = self._attributes
        return gaFeatureSet

    def getCompoundIdForFeatureId(self, featureId):
        """
        :param featureId: string feature ID as reported in GFF3
        :return: string representing ID for the specified GA4GH protocol
            Feature object in this FeatureSet.
        """
        compoundId = datamodel.FeatureCompoundId(
            self.getCompoundId(), _decolonize(featureId))
        return str(compoundId)


class SimulatedFeatureSet(AbstractFeatureSet):
    def __init__(self, parentContainer, localId, randomSeed=1):
        super(SimulatedFeatureSet, self).__init__(
            parentContainer, localId, None)
        self._randomSeed = randomSeed

    def _getRandomOntologyTerm(self, randomNumberGenerator):
        ontologyTuples = [("gene", "SO:0000704"),
                          ("exon", "SO:0000147")]
        term = protocol.OntologyTerm()
        ontologyTuple = randomNumberGenerator.choice(ontologyTuples)
        term.term, term.id = ontologyTuple[0], ontologyTuple[1]
        term.sourceName = "sequenceOntology"
        term.sourceVersion = "0"
        return term

    def _generateSimulatedFeature(self, randomNumberGenerator):
        feature = protocol.Feature()
        # ID and parentId are not generated here.
        feature.featureSetId = self.getId()
        feature.start = randomNumberGenerator.randint(0, 10000)
        feature.end = feature.start + randomNumberGenerator.randint(1, 1000)
        feature.featureType = self._getRandomOntologyTerm(
            randomNumberGenerator)
        references = ["chr1", "chr2", "chrX"]
        feature.referenceName = randomNumberGenerator.choice(references)
        strands = ["+", "-"]
        feature.strand = randomNumberGenerator.choice(strands)
        feature.attributes = {"gene_name": ["Frances", ],
                              "gene_type": ["mRNA", ],
                              "gene_status": ["UNKNOWN", ]}
        return feature

    def getFeature(self, compoundId):
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(self._randomSeed)
        feature = self._generateSimulatedFeature(randomNumberGenerator)
        feature.id = compoundId
        feature.parentId = None  # FIXME
        return feature

    def getFeatures(self, referenceName, start, end,
                    pageToken, pageSize,
                    ontologyTerms=[], parentId=None):
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(self._randomSeed)
        numFeatures = 10
        if pageToken is not None:
            nextPageToken = int(pageToken)
        else:
            nextPageToken = 0
        for featureId in ["feature{}".format(x) for x in xrange(numFeatures)]:
            gaFeature = self._generateSimulatedFeature(randomNumberGenerator)
            gaFeature.id = featureId
            gaFeature.parentId = None  # FIXME
            if nextPageToken < numFeatures - 1:
                nextPageToken += 1
            else:
                nextPageToken = None
            yield gaFeature, nextPageToken


class Gff3DbFeatureSet(AbstractFeatureSet):
    """
    Stub class to directly read sequence annotation features from GFF3 files.
    Tests basic access, not to be used in production.
    """
    def __init__(self, parentContainer, localId, filePath, dataRepository):
        super(Gff3DbFeatureSet, self).__init__(
            parentContainer, localId, None)
        self._sequenceOntology = dataRepository.getOntology(
            'sequence_ontology')
        self._dbFilePath = filePath
        self._dataRepository = dataRepository
        self._db = Gff3DbBackend(self._dbFilePath)

    def getFeature(self, compoundId):
        featureId = _recolonize(compoundId.featureId)
        with self._db as dataSource:
            featureReturned = dataSource.getFeatureById(featureId)

        if featureReturned is not None:
            gaFeature = self._gaFeatureForFeatureDbRecord(featureReturned)
            return gaFeature
        else:
            return None

    def _gaFeatureForFeatureDbRecord(self, feature):
        """
        :param feature: The DB Row representing a feature
        :return: the corresponding GA4GH protocol.Feature object
        """
        gaFeature = protocol.Feature()
        gaFeature.id = self.getCompoundIdForFeatureId(feature['id'])
        if feature.get('parent_id'):
            gaFeature.parentId = self.getCompoundIdForFeatureId(
                    feature['parent_id'])
        else:
            gaFeature.parentId = None
        gaFeature.featureSetId = self.getId()
        gaFeature.referenceName = feature['reference_name']
        gaFeature.start = int(feature['start'])
        gaFeature.end = int(feature['end'])
        gaFeature.childIds = map(
                self.getCompoundIdForFeatureId,
                json.loads(feature['child_ids']))
        gaFeature.featureType = \
            self._sequenceOntology.getGaTermByName(feature['ontology_term'])
        gaFeature.attributes = json.loads(
            feature['attributes'])
        return gaFeature

    def getFeatures(self, referenceName, start, end,
                    pageToken, pageSize,
                    ontologyTerms=[], parentId=None):
        """
        method passed to runSearchRequest to fulfill the request
        :param str referenceName: name of reference (ex: "chr1")
        :param start: castable to int, start position on reference
        :param end: castable to int, end position on reference
        :param pageToken: none or castable to int
        :param pageSize: none or castable to int
        :param ontologyTerms: array of str
        :param parentId: none or featureID of parent
        :return: yields a protocol.Feature at a time, together with
            the corresponding nextPageToken (which is null for the last
            feature served out).
        """
        # parse out the various query parameters from the request.
        start = int(start)
        end = int(end)

        with self._db as dataSource:
            # featuresCount is needed to ensure that once the
            # request is fulfilled, no nextPageTokens past the
            # end of the actual dataset range are returned.
            featuresCount = dataSource.countFeaturesSearchInDb(
                referenceName=referenceName,
                start=start, end=end,
                parentId=parentId, ontologyTerms=ontologyTerms)
            featuresReturned = dataSource.searchFeaturesInDb(
                pageToken, pageSize,
                referenceName=referenceName,
                start=start, end=end,
                parentId=parentId, ontologyTerms=ontologyTerms)

        # pagination logic: None if last feature was returned,
        # else 1 + row number being returned (starting at row 0).
        if pageToken is not None:
            nextPageToken = int(pageToken)
        else:
            nextPageToken = 0
        for featureRecord in featuresReturned:
            gaFeature = self._gaFeatureForFeatureDbRecord(featureRecord)
            if nextPageToken < featuresCount - 1:
                nextPageToken += 1
            else:
                nextPageToken = None
            yield gaFeature, nextPageToken
