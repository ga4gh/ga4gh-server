"""
Module responsible for translating g2p data into GA4GH native
objects.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import rdflib
import urlparse
import sets
import os
import types

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


class G2PDataset:
    """
    An rdf object store.  The cancer genome database
    [Clinical Genomics Knowledge Base]
    (http://nif-crawler.neuinfo.org/monarch/ttl/cgd.ttl)
    published by the Monarch project was the source of Evidence.
    """

    def __init__(self, setName=None, relativePath=None, dataDir=None):
        """
        Initialize dataset, using the passed dict of sources
        [{source,format}] see rdflib.parse() for more
        If path is set, this backend will load itself
        """
        self._sources = None
        if setName is not None:
            self._sources = [os.path.join(relativePath, f)
                             for f in os.listdir(relativePath)]
        self._searchQuery = """
        PREFIX OBAN: <http://purl.org/oban/>
        PREFIX OBO: <http://purl.obolibrary.org/obo/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX faldo: <http://biohackathon.org/resource/faldo#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        SELECT  %PROPERTIES%
            WHERE {
                ?s  a OBAN:association .
                ?s  OBAN:association_has_subject ?drug .
                ?s  OBAN:association_has_object  ?disease .
                ?s  OBO:RO_has_approval_status ?evidence .
                ?disease OBO:RO_has_biomarker ?location .
                ?drug  rdfs:label ?drug_label  .
                ?disease  rdfs:label ?disease_label  .
                ?disease rdf:type ?disease_type .
                ?evidence  rdfs:label ?evidence_label  .
                ?location  rdfs:label ?location_label  .
                ?disease_type  rdfs:label ?disease_type_label  .
                %FILTER%
        }
        ORDER BY ?s
            """

        # initialize graph
        self._rdfGraph = rdflib.ConjunctiveGraph()

        # load with data
        if self._sources is None:
            self._rdfGraph.parse('tests/data/g2pDatasets/cgd/cgd-test.ttl',
                                 format='n3')
        else:
            for source in self._sources:
                if source.endswith('.ttl'):
                    self._rdfGraph.parse(source, format='n3')
                else:
                    self._rdfGraph.parse(source, format='xml')

        # log queries that take more than N seconds
        # import time
        # self.RESPONSETIME_LOGGING_THRESHOLD = 2

    def queryLabels(
            self, location=None, drug=None, disease=None, pageSize=None,
            offset=0):

        """
        This query is the main search mechanism.
        It queries the graph for annotations that match the
        AND of [location,drug,disease].
        """
        query = self._formatQuery(location, drug, disease)

        query = query.replace("%PROPERTIES%", "".join(
            ["distinct ?s ?location ?location_label ",
             "?disease ?disease_label ?drug ?drug_label ",
             "?evidence ?evidence_label"]))

        # query += ("LIMIT {} OFFSET {} ".format(pageSize, offset))

        # now = time.time()
        results = self._rdfGraph.query(query)
        # Depending on the cardinality this query can return multiple rows
        # per annotations.  Here we reduce it to a list of unique annotations
        # URIs
        uniqueAnnotations = sets.Set()
        for row in results:
            uniqueAnnotations.add("<{}>".format(row['s'].toPython()))

        annotations = AnnotationFactory(self._rdfGraph,
                                        uniqueAnnotations).annotations()

        # responseTime = time.time()-now
        # if responseTime > self.RESPONSETIME_LOGGING_THRESHOLD:
        #     print('queryLabels', responseTime)
        #     print('len(annotations)', len(annotations))

        return annotations

    def _formatQuery(self, location=None, drug=None, disease=None):
        """
        Generate a formatted sparql query with appropriate filters
        """
        query = self._searchQuery
        if location is None and drug is None and disease is None:
            raise exceptions.NotImplementedException(
               "At least one of [location, drug, disease] must be specified")
        filters = []

        if location and isinstance(location, basestring):
            filters.append('regex(?location_label, "{}")'.format(location))
        if drug and isinstance(drug, basestring):
            filters.append('regex(?drug_label, "{}")'.format(drug))
        if disease and isinstance(disease, basestring):
            filters.append('regex(?disease_label, "{}")'.format(disease))

        locationClause = ""
        if isinstance(location, dict):
            locations = []
            for id in location['ids']:
                    locations.append('?location = <{}> '.format
                                     (id['database'] + id['identifier']))
            locationClause = "({})".format(" || ".join(locations))
            filters.append(locationClause)
            locationClause = "?l  faldo:location ?location .\n"

        if isinstance(drug, dict):
            drugs = []
            for id in drug['ids']:
                    drugs.append('?drug = <{}> '.format
                                 (id['database'] + id['identifier']))
            drugsClause = "({})".format(" || ".join(drugs))

            filters.append(drugsClause)

        if isinstance(disease, dict):
            diseases = []
            for id in disease['ids']:
                    diseases.append('?disease = <{}> '.format
                                    (id['database'] + id['identifier']))
            diseasesClause = "({})".format(" || ".join(diseases))
            filters.append(diseasesClause)

        filter = "FILTER ({})".format(' && '.join(filters))
        query = query.replace("%FILTER%", filter)
        query = query.replace("%LOCATION%", locationClause)
        return query


class AnnotationFactory:
    """
    Given a RDF query result set, return corresponding set of ProtocolElements
    """

    def __init__(self, graph, uniqueAnnotations):
        self._rdfGraph = graph

        # now fetch the details on the annotation
        self._annotations = []
        for annotation in uniqueAnnotations:
            self._annotations.append(
                self._toGA4GH(self._query(annotation)))

        # add a callback to return ProtocolElement
        # since we have already transformed it into ProtocolElement
        # we just return self
        def toProtocolElement(self):
            return self

        for annotation in self._annotations:
            annotation.toProtocolElement = \
                types.MethodType(toProtocolElement, annotation)

    def annotations(self):
        """
        return annotations built in constructor
        """
        return self._annotations

    def _query(self, subject=''):
        """
        This is the 'detail' query

        Return a list of dictionaries.
        Each dict is an [annotation](http://www.w3.org/ns/oa#Annotation)
        All keys in the dict are predicates of the annotation.
        All cells in the dict are predicate->object converted to strings.

        If an annotation has a <http://purl.org/oban/association_has_subject>
        predicate that class is appended to the annotation dict in the
        'location' property
        """

        annotationQuery = """
        PREFIX OBAN: <http://purl.org/oban/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        SELECT distinct *
          WHERE {
          %s ?p ?o .
          OPTIONAL {?o  rdfs:label ?label .}
        }
        """ % subject

        # now = time.time()
        results = self._rdfGraph.query(annotationQuery)
        rows = [row.asdict() for row in results]

        # responseTime = time.time()-now
        # if responseTime > self.RESPONSETIME_LOGGING_THRESHOLD:
        #     print('annotationQuery', responseTime)
        #     print(annotationQuery)

        for row in rows:
            for k in row:
                row[k] = row[k].toPython()
            row['s'] = subject

        locationQueryTemplate = """
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
             PREFIX OBO: <http://purl.obolibrary.org/obo/>
        SELECT distinct *
        WHERE {
        %SUBJECT%    a  OBO:OMIM_606764 .
        %SUBJECT%   ?p ?o .
        OPTIONAL {?o  rdfs:label ?label .} .
        }
        """

        locationRows = []
        uniqueLocations = sets.Set()
        for row in rows:
            if row['p'] == 'http://purl.org/oban/association_has_object':
                uniqueLocations.add("<" + row['o'] + ">")

        for location in uniqueLocations:
            locationQuery = locationQueryTemplate.replace(
                "%SUBJECT%", location)
            results = self._rdfGraph.query(locationQuery)
            locationRows = [locationRow.asdict() for locationRow in results]
            for locationRow in locationRows:
                for k in locationRow:
                    locationRow[k] = locationRow[k].toPython()
                locationRow['s'] = location
            # if responseTime > self.RESPONSETIME_LOGGING_THRESHOLD:
            #     print('locationQuery', responseTime)
            #     print(locationQuery)

        annotation = self._flatten(rows)
        location = self._flatten(locationRows)
        annotation['location'] = location
        return annotation

    def _flatten(self, dict):
        """
        Given a dict of dicts,
        flatten it to a single dict using predicate as keys
        For multiple occurrences of a predicate, create an array
        Each value in the dict is an object {val:'x', label:'y'}
        The value of 's' (subject) is copied to the 'id' property
        """
        a = {}
        for row in dict:
            obj = {'val': row['o']}
            if 'label' in row:
                obj['label'] = row['label']

            if row['p'] in a and \
                    a[row['p']].__class__.__name__ != "list":
                asIs = a[row['p']]
                a[row['p']] = []
                a[row['p']].append(asIs)

            if row['p'] in a:
                a[row['p']].append(obj)
            else:
                a[row['p']] = obj

            a['id'] = row['s']
        return a

    def _toGA4GH(self, annotation):
        """
        given an annotation dict, return a protocol.FeaturePhenotypeAssociation
        """

        fpa = None

        # annotation keys
        location = 'location'
        hasObject = 'http://purl.org/oban/association_has_object'
        has_approval_status = \
            'http://purl.obolibrary.org/obo/RO_has_approval_status'
        # location keys
        GENO_0000408 = 'http://purl.obolibrary.org/obo/GENO_0000408'

        location = annotation['location']
        if GENO_0000408 in location:
            id_, ontologySource = self.namespaceSplit(
                                        location[GENO_0000408]['val'])
            name = location[GENO_0000408]['label']
        else:
            id_, ontologySource = self.namespaceSplit(location['id'])
            name = location['id']

        f = protocol.Feature()
        f.featureType = protocol.OntologyTerm.fromJsonDict({
            "name": name,
            "id": id_,
            "ontologySource": ontologySource})
        f.id = annotation['id']
        f.featureSetId = ''
        f.parentIds = []
        f.attributes = protocol.Attributes.fromJsonDict({"vals": {}})

        # # debugger example how to validate and capture validation errors
        # if not protocol.Feature.validate(f.toJsonDict()):
        #     e = exceptions.RequestValidationFailureException(
        #         f.toJsonDict(),protocol.Feature)
        #     print(e.message)
        #     from IPython.core.debugger import Pdb ;        Pdb().set_trace()

        id_, ontologySource = self.namespaceSplit(
                                       annotation[hasObject]['val'])

        fpa = protocol.FeaturePhenotypeAssociation()
        fpa.id = annotation['id']
        fpa.features = [f]
        fpa.description = None
        fpa.evidence = []
        fpa.environmentalContexts = []

        phenotypeInstance = protocol.PhenotypeInstance()
        phenotypeInstance.type = protocol.OntologyTerm.fromJsonDict({
            "name": annotation[hasObject]['label'],
            "id": id_,
            "ontologySource": ontologySource})
        fpa.phenotype = phenotypeInstance

        #  ECO or OBI is recommended
        if has_approval_status in annotation:
            approval_status = annotation[has_approval_status]
            evidence = protocol.Evidence()
            evidence.evidenceType = protocol.OntologyTerm()
            id_, ontology_source = self.namespaceSplit(approval_status['val'])
            evidence.evidenceType.ontologySource = ontology_source
            evidence.evidenceType.id = id_

            evidence.evidenceType.name = ''
            if 'label' in approval_status:
                evidence.evidenceType.name = approval_status['label']
                fpa.evidence.append(evidence)
                if not protocol.Evidence.validate(evidence.toJsonDict()):
                    raise exceptions.RequestValidationFailureException(
                        evidence.toJsonDict(), protocol.Evidence)

        return fpa

    def namespaceSplit(self, url, separator='/'):
        """
        given a url return the id of the resource and the ontology source
        """
        o = urlparse.urlparse(url)
        _id = o.path.split(separator)[-1]
        ontologySource = urlparse.urlunsplit([o[0],
                                              o[1],
                                              o[2].replace(_id, ''),
                                              o[4], ''])
        return _id, ontologySource
