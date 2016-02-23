"""
Module responsible for translating g2p data into GA4GH native
objects.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import rdflib
import urlparse
import os
import types

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
import ga4gh.datamodel as datamodel


class AbstractG2PDataset(datamodel.DatamodelObject):

    def __init__(self, setName=None):
        self._searchQuery = ""

    def toProtocolElement(self):
        # TODO remove this... added to follow pattern of datadriven tests
        fpa = protocol.FeaturePhenotypeAssociation()
        fpa.features = []
        fpa.id = "TEST"
        fpa.evidence = []
        fpa.phenotype = protocol.PhenotypeInstance()
        fpa.phenotype.type = protocol.OntologyTerm()
        fpa.phenotype.type.id = "fakeid"
        fpa.phenotype.type.ontologySource = "fakesource"
        fpa.environmentalContexts = []
        return fpa


class SimulatedG2PDataset(AbstractG2PDataset):
    def __init__(self, setName=None, relativePath=None, dataDir=None):
        pass

    def queryLabels(self, location=None, drug=None,
                    disease=None, pageSize=None, offset=0):
        fpa = self.toProtocolElement()
        fpas = []
        if location or drug or disease:
            fpas.append(fpa)
        # TODO this seems like a hack

        def toProtocolElement(self):
            return self

        for fpa in fpas:
            fpa.toProtocolElement = \
                types.MethodType(toProtocolElement, fpa)
        return fpas


class G2PDataset(AbstractG2PDataset):
    """
    An rdf object store.  The cancer genome database
    [Clinical Genomics Knowledge Base]
    (http://nif-crawler.neuinfo.org/monarch/ttl/cgd.ttl),
    published by the Monarch project, was the source of Evidence.
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

    def _formatQuery(self, location=None, drug=None, disease=None):
        """
        Generate a formatted sparql query with appropriate filters
        """
        query = self._searchQuery
        if location is None and drug is None and disease is None:
            # TODO is this really the exception we want to throw?
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

            for _id in location['ids']:
                    locations.append('?location = <{}> '.format
                                     (_id['database'] + _id['identifier']))
            locationClause = "({})".format(" || ".join(locations))
            filters.append(locationClause)
            locationClause = "?l  faldo:location ?location .\n"

        if isinstance(drug, dict):
            drugs = []
            for _id in drug['ids']:
                    drugs.append('?drug = <{}> '.format
                                 (_id['database'] + _id['identifier']))
            drugsClause = "({})".format(" || ".join(drugs))

            filters.append(drugsClause)

        if isinstance(disease, dict):
            diseases = []
            for _id in disease['ids']:
                    diseases.append('?disease = <{}> '.format
                                    (_id['database'] + _id['identifier']))
            diseasesClause = "({})".format(" || ".join(diseases))
            filters.append(diseasesClause)

        filter = "FILTER ({})".format(' && '.join(filters))
        query = query.replace("%FILTER%", filter)
        query = query.replace("%LOCATION%", locationClause)
        return query

    def queryLabels(
            self, location=None, drug=None, disease=None, pageSize=None,
            offset=0):

        """
        This query is the main search mechanism.
        It queries the graph for annotations that match the
        AND of [location,drug,disease].
        """
        query = self._formatQuery(location, drug, disease)
        # TODO refactor to testable function
        query = query.replace("%PROPERTIES%", "".join(
            ["distinct ?s ?location ?location_label ",
             "?disease ?disease_label ?drug ?drug_label ",
             "?evidence ?evidence_label"]))

        # TODO why is this commented out?
        # query += ("LIMIT {} OFFSET {} ".format(pageSize, offset))

        results = self._rdfGraph.query(query)
        # Depending on the cardinality this query can return multiple rows
        # per annotations.  Here we reduce it to a list of unique annotations
        # URIs
        uniqueAnnotations = set()
        for row in results:
            uniqueAnnotations.add("<{}>".format(row['s'].toPython()))

        annotations = AnnotationFactory(self._rdfGraph,
                                        uniqueAnnotations).annotations()

        return annotations

    def toProtocolElement(self):
        # TODO remove this... added to follow pattern of datadriven tests
        fpa = protocol.FeaturePhenotypeAssociation()
        fpa.features = []
        fpa.id = "TEST"
        fpa.evidence = []
        fpa.phenotype = protocol.PhenotypeInstance()
        fpa.phenotype.type = protocol.OntologyTerm()
        fpa.phenotype.type.id = "fakeid"
        fpa.phenotype.type.ontologySource = "fakesource"
        fpa.environmentalContexts = []
        return fpa


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

        results = self._rdfGraph.query(annotationQuery)
        rows = [row.asdict() for row in results]

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
        uniqueLocations = set()
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

        annotation = self._flatten(rows)
        location = self._flatten(locationRows)
        annotation['location'] = location
        return annotation

    def _flatten(self, dictList):
        """
        Given a list of dicts, with form
        [
            {
             'p': predicate,
             's': subject,
             'o': object,
             'label': label  # optional
            }
        ]

        flatten it to a single dict using predicate as keys
        For multiple occurrences of a predicate, create an array
        Each value in the dict is an object {val:'x', label:'y'}
        The value of 's' (subject) is copied to the 'id' property
        """
        # TODO what if 's' is multiply valued?
        a = {}
        for row in dictList:
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
        locationKey = 'location'
        hasObject = 'http://purl.org/oban/association_has_object'
        has_approval_status = \
            'http://purl.obolibrary.org/obo/RO_has_approval_status'
        # location keys
        GENO_0000408 = 'http://purl.obolibrary.org/obo/GENO_0000408'

        location = annotation[locationKey]
        if GENO_0000408 in location:
            id_, ontologySource = self.namespaceSplit(
                                        location[GENO_0000408]['val'])
            name = location[GENO_0000408]['label']
        else:
            id_, ontologySource = self.namespaceSplit(location['id'])
            name = location['id']

        f = protocol.Feature()
        f.featureType = protocol.OntologyTerm.fromJsonDict({
            "term": name,
            "id": id_,
            "sourceName": ontologySource})
        f.id = annotation['id']  # TODO connect with real feature Ids
        f.featureSetId = ''
        f.parentIds = []
        f.attributes = protocol.Attributes.fromJsonDict({"vals": {}})

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
            "term": annotation[hasObject]['label'],
            "id": id_,
            "sourceName": ontologySource})
        fpa.phenotype = phenotypeInstance

        #  ECO or OBI is recommended
        if has_approval_status in annotation:
            approval_status = annotation[has_approval_status]
            evidence = protocol.Evidence()
            evidence.evidenceType = protocol.OntologyTerm()
            id_, ontology_source = self.namespaceSplit(approval_status['val'])
            evidence.evidenceType.sourceName = ontology_source
            evidence.evidenceType.id = id_

            evidence.evidenceType.term = ''
            if 'label' in approval_status:
                evidence.evidenceType.term = approval_status['label']
                fpa.evidence.append(evidence)
                if not protocol.Evidence.validate(evidence.toJsonDict()):
                    raise exceptions.RequestValidationFailureException(
                        evidence.toJsonDict(), protocol.Evidence)
        return fpa

    def namespaceSplit(self, url, separator='/'):
        """
        given a url return the id of the resource and the ontology source
        """
        url = url.strip("<").strip(">")
        o = urlparse.urlparse(url)
        _id = o.path.split(separator)[-1]
        ontologySource = urlparse.urlunsplit([o[0],
                                              o[1],
                                              o[2].replace(_id, ''),
                                              o[4], ''])
        return _id, ontologySource
