"""
Module responsible for translating g2p data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import re
import bisect
import rdflib
from rdflib import RDF

import ga4gh.server.exceptions as exceptions
import ga4gh.server.datamodel.sequence_annotations as sequence_annotations
import ga4gh.server.datamodel.genotype_phenotype as g2p

import ga4gh.schemas.protocol as protocol

# annotation keys
TYPE = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type'
LABEL = 'http://www.w3.org/2000/01/rdf-schema#label'
HAS_QUALITY = 'http://purl.obolibrary.org/obo/BFO_0000159'
FALDO_LOCATION = "http://biohackathon.org/resource/faldo#location"
FALDO_BEGIN = "http://biohackathon.org/resource/faldo#begin"
FALDO_END = "http://biohackathon.org/resource/faldo#end"
FALDO_POSITION = "http://biohackathon.org/resource/faldo#position"
FALDO_REFERENCE = "http://biohackathon.org/resource/faldo#reference"
MEMBER_OF = 'http://purl.obolibrary.org/obo/RO_0002350'
ASSOCIATION = "http://purl.org/oban/association"
HAS_SUBJECT = "http://purl.org/oban/association_has_subject"


class PhenotypeAssociationFeatureSet(
        g2p.G2PUtility, sequence_annotations.Gff3DbFeatureSet):
    """
    An rdf object store.  The cancer genome database
    [Clinical Genomics Knowledge Base]
    (http://nif-crawler.neuinfo.org/monarch/ttl/cgd.ttl),
    published by the Monarch project, was the source of Evidence.
    """

    def __init__(self, parentContainer, localId):
        super(PhenotypeAssociationFeatureSet, self).__init__(
            parentContainer, localId)

    # mimic featureset
    def populateFromRow(self, featureSetRecord):
        """
        Populates the instance variables of this FeatureSet from the specified
        DB row.
        """
        self._dbFilePath = featureSetRecord.dataurl
        self.setAttributesJson(featureSetRecord.attributes)
        self.populateFromFile(self._dbFilePath)

    def populateFromFile(self, dataUrl):
        """
        Populates the instance variables of this FeatureSet from the specified
        data URL.
        Initialize dataset, using the passed dict of sources
        [{source,format}] see rdflib.parse() for more
        If path is set, this backend will load itself
        """
        self._dbFilePath = dataUrl
        # initialize graph
        self._rdfGraph = rdflib.ConjunctiveGraph()
        # save the path
        self._dataUrl = dataUrl
        self._scanDataFiles(self._dataUrl, ['*.ttl'])

        # extract version
        cgdTTL = rdflib.URIRef("http://data.monarchinitiative.org/ttl/cgd.ttl")
        versionInfo = rdflib.URIRef(
            u'http://www.w3.org/2002/07/owl#versionInfo')
        self._version = None
        for _, _, obj in self._rdfGraph.triples((cgdTTL, versionInfo, None)):
            self._version = obj.toPython()

        # setup location cache
        self._initializeLocationCache()

    # mimic featureset
    def getFeature(self, compoundId):
        """
        find a feature and return ga4gh representation, use compoundId as
        featureId
        """
        feature = self._getFeatureById(compoundId.featureId)
        feature.id = str(compoundId)
        return feature

    def _getFeatureById(self, featureId):
        """
        find a feature and return ga4gh representation, use 'native' id as
        featureId
        """
        featureRef = rdflib.URIRef(featureId)
        featureDetails = self._detailTuples([featureRef])
        feature = {}
        for detail in featureDetails:
            feature[detail['predicate']] = []

        for detail in featureDetails:
            feature[detail['predicate']].append(detail['object'])

        pbFeature = protocol.Feature()

        term = protocol.OntologyTerm()
        # Schema for feature only supports one type of `type`
        # here we default to first OBO defined
        for featureType in sorted(feature[TYPE]):
            if "obolibrary" in featureType:
                term.term = self._featureTypeLabel(featureType)
                term.term_id = featureType
                pbFeature.feature_type.MergeFrom(term)
                break

        pbFeature.id = featureId
        # Schema for feature only supports one type of `name` `symbol`
        # here we default to shortest for symbol and longest for name
        feature[LABEL].sort(key=len)
        pbFeature.gene_symbol = feature[LABEL][0]
        pbFeature.name = feature[LABEL][-1]

        pbFeature.attributes.MergeFrom(protocol.Attributes())
        for key in feature:
            for val in sorted(feature[key]):
                pbFeature.attributes.attr[key].values.add().string_value = val

        if featureId in self._locationMap:
            location = self._locationMap[featureId]
            pbFeature.reference_name = location["chromosome"]
            pbFeature.start = location["begin"]
            pbFeature.end = location["end"]

        return pbFeature

    # mimic featureset
    def getFeatures(self, referenceName=None, start=None, end=None,
                    startIndex=None, maxResults=None,
                    featureTypes=None, parentId=None,
                    name=None, geneSymbol=None, numFeatures=10):

        # query to do search
        query = self._filterSearchFeaturesRequest(
            referenceName, geneSymbol, name, start, end)
        featuresResults = self._rdfGraph.query(query)
        featureIds = set()
        try:
            for row in featuresResults.bindings:
                featureIds.add(row['feature'].toPython())
        except re.error:
            raise exceptions.BadFeatureSetSearchRequestRegularExpression()

        if startIndex:
            startPosition = int(startIndex)
        else:
            startPosition = 0
        for i, featureId in enumerate(featureIds):
            if i < startPosition:
                continue
            feature = self._getFeatureById(featureId)
            # _getFeatureById returns native id, cast to compound
            feature.id = self.getCompoundIdForFeatureId(feature.id)
            yield feature

    def _baseQuery(self):
        return """
        PREFIX OBAN: <http://purl.org/oban/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        SELECT DISTINCT
                        ?feature
                        ?feature_label
        WHERE {
                   ?association  a OBAN:association .
                   ?association    OBAN:association_has_subject ?feature .
                   ?feature  rdfs:label ?feature_label  .
                   #%FILTER%
        }
        ORDER BY ?feature
        """

    def _filterSearchFeaturesRequest(self, reference_name, gene_symbol, name,
                                     start, end):
        """
        formulate a sparql query string based on parameters
        """
        filters = []
        query = self._baseQuery()
        filters = []
        location = self._findLocation(reference_name, start, end)
        if location:
            filters.append("?feature = <{}>".format(location))
        if gene_symbol:
            filters.append('regex(?feature_label, "{}")')
        if name:
            filters.append(
                'regex(?feature_label, "{}")'.format(name))
        # apply filters
        filter = "FILTER ({})".format(' && '.join(filters))
        if len(filters) == 0:
            filter = ""
        query = query.replace("#%FILTER%", filter)
        return query

    def _findLocation(self, reference_name, start, end):
        """
        return a location key form the locationMap
        """
        try:
            # TODO - sequence_annotations does not have build?
            return self._locationMap['hg19'][reference_name][start][end]
        except:
            return None

    def _initializeLocationCache(self):
        """
        CGD uses Faldo ontology for locations, it's a bit complicated.
        This function sets up an in memory cache of all locations, which
        can be queried via:
        locationMap[build][chromosome][begin][end] = location["_id"]
        """
        # cache of locations
        self._locationMap = {}
        locationMap = self._locationMap
        triples = self._rdfGraph.triples
        Ref = rdflib.URIRef

        associations = []
        for subj, _, _ in triples((None, RDF.type, Ref(ASSOCIATION))):
            associations.append(subj.toPython())

        locationIds = []
        for association in associations:
            for _, _, obj in triples((Ref(association),
                                      Ref(HAS_SUBJECT), None)):
                locationIds.append(obj.toPython())

        locations = []
        for _id in locationIds:
            location = {}
            location["_id"] = _id
            for subj, predicate, obj in triples((Ref(location["_id"]),
                                                 None, None)):
                if not predicate.toPython() in location:
                    location[predicate.toPython()] = []
                bisect.insort(location[predicate.toPython()], obj.toPython())
                if FALDO_LOCATION in location:
                    locations.append(location)

        for location in locations:
            for _id in location[FALDO_LOCATION]:
                # lookup faldo region, ensure positions are sorted
                faldoLocation = {}
                faldoLocation["_id"] = _id
                for subj, predicate, obj in triples((Ref(faldoLocation["_id"]),
                                                    None, None)):
                    if not predicate.toPython() in faldoLocation:
                        faldoLocation[predicate.toPython()] = []
                    bisect.insort(faldoLocation[predicate.toPython()],
                                  obj.toPython())

                faldoBegins = []

                for _id in faldoLocation[FALDO_BEGIN]:
                    faldoBegin = {}
                    faldoBegin["_id"] = _id
                    for subj, predicate, obj in triples(
                                                (Ref(faldoBegin["_id"]),
                                                    None, None)):
                        faldoBegin[predicate.toPython()] = obj.toPython()
                    faldoBegins.append(faldoBegin)

                faldoReferences = []
                for _id in faldoLocation[FALDO_BEGIN]:
                    faldoReference = {}
                    faldoReference["_id"] = faldoBegin[FALDO_REFERENCE]
                    for subj, predicate, obj in triples(
                                                (Ref(faldoReference["_id"]),
                                                    None, None)):
                        faldoReference[predicate.toPython()] = obj.toPython()
                    faldoReferences.append(faldoReference)

                faldoEnds = []
                for _id in faldoLocation[FALDO_END]:
                    faldoEnd = {}
                    faldoEnd["_id"] = _id
                    for subj, predicate, obj in triples((Ref(faldoEnd["_id"]),
                                                        None, None)):
                        faldoEnd[predicate.toPython()] = obj.toPython()
                    faldoEnds.append(faldoEnd)

                for idx, faldoReference in enumerate(faldoReferences):
                    if MEMBER_OF in faldoReference:
                        build = faldoReference[MEMBER_OF].split('/')[-1]
                        chromosome = faldoReference[LABEL].split(' ')[0]
                        begin = faldoBegins[idx][FALDO_POSITION]
                        end = faldoEnds[idx][FALDO_POSITION]
                        if build not in locationMap:
                            locationMap[build] = {}
                        if chromosome not in locationMap[build]:
                            locationMap[build][chromosome] = {}
                        if begin not in locationMap[build][chromosome]:
                            locationMap[build][chromosome][begin] = {}
                        if end not in locationMap[build][chromosome][begin]:
                            locationMap[build][chromosome][begin][end] = {}
                        locationMap[build][chromosome][begin][end] = \
                            location["_id"]
                        locationMap[location["_id"]] = {
                            "build": build,
                            "chromosome": chromosome,
                            "begin": begin,
                            "end": end,
                        }
