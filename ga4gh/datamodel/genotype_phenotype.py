"""
Module responsible for translating g2p data into GA4GH native
objects.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import rdflib
import urlparse

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
import ga4gh.datamodel as datamodel


class AbstractPhenotypeAssociationSet(datamodel.DatamodelObject):
    compoundIdClass = datamodel.PhenotypeAssociationSetCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractPhenotypeAssociationSet, self).__init__(
            parentContainer, localId)

    def toProtocolElement(self):
        pas = protocol.PhenotypeAssociationSet()
        pas.name = self.getLocalId()
        pas.id = self.getId()
        pas.datasetId = self.getParentContainer().getId()
        pas.info = self.getInfo()
        return pas

    def getInfo(self):
        return {}


class SimulatedPhenotypeAssociationSet(AbstractPhenotypeAssociationSet):
    def __init__(self, parentContainer, localId, randomSeed):
        super(SimulatedPhenotypeAssociationSet, self).__init__(
            parentContainer, localId)

    def getAssociations(
            self, location=None, drug=None, disease=None, pageSize=None,
            offset=0):
        if location or drug or disease:
            fpa = protocol.FeaturePhenotypeAssociation()
            fpa.id = "test"
            fpa.features = []
            fpa.evidence = []
            fpa.environmentalContexts = []
            fpa.phenotype = protocol.PhenotypeInstance()
            fpa.phenotype.type = protocol.OntologyTerm()
            fpa.phenotype.type.id = "test"
            fpa.phenotypeAssociationSetId = self.getId()
            return [fpa]
        else:
            return []


class PhenotypeAssociationSet(AbstractPhenotypeAssociationSet):
    """
    An rdf object store.  The cancer genome database
    [Clinical Genomics Knowledge Base]
    (http://nif-crawler.neuinfo.org/monarch/ttl/cgd.ttl),
    published by the Monarch project, was the source of Evidence.
    """

    def __init__(self, parentContainer, localId, dataDir):
        super(PhenotypeAssociationSet, self).__init__(parentContainer, localId)
        """
        Initialize dataset, using the passed dict of sources
        [{source,format}] see rdflib.parse() for more
        If path is set, this backend will load itself
        """

        # initialize graph
        self._rdfGraph = rdflib.ConjunctiveGraph()
        try:
            self._scanDataFiles(dataDir, ['*.ttl', '*.xml'])
        except AttributeError:
            pass

    def getAssociations(
            self, location=None, drug=None, disease=None, pageSize=None,
            offset=0):
        """
        This query is the main search mechanism.
        It queries the graph for annotations that match the
        AND of [location,drug,disease].
        """
        query = self._formatFilterQuery(location, drug, disease)

        results = self._rdfGraph.query(query)
        # Depending on the cardinality this query can return multiple rows
        # per annotations.  Here we reduce it to a list of unique annotations
        # URIs

        # TODO whys it a set? we called distinct
        # TODO paging using SPARQL?
        uniqueAnnotations = set()
        for row in results:
            uniqueAnnotations.add("<{}>".format(row['s'].toPython()))

        annotations = []
        for annotation in uniqueAnnotations:
            annotations.append(
                self._toGA4GH(self._detailQuery(annotation)))

        return annotations

    def _addDataFile(self, filename):
        if filename.endswith('.ttl'):
            self._rdfGraph.parse(filename, format='n3')
        else:
            self._rdfGraph.parse(filename, format='xml')

    def _formatFilterQuery(self, location=None, drug=None, disease=None):
        """
        Generate a formatted sparql query with appropriate filters
        """
        query = """
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
        if location is None and drug is None and disease is None:
            # TODO is this really the exception we want to throw?
            raise exceptions.NotImplementedException(
               "At least one of [location, drug, disease] must be specified")
        filters = []

        if location and isinstance(location, basestring):
            filters.append('regex(?location_label, "{}")'.format(location))
            # TODO check this regex (versus string match)
        if drug and isinstance(drug, basestring):
            filters.append('regex(?drug_label, "{}")'.format(drug))
        if disease and isinstance(disease, basestring):
            filters.append('regex(?disease_label, "{}")'.format(disease))

        locationClause = ""
        if isinstance(location, dict):
            locations = []

            for id_ in location['ids']:
                    locations.append('?location = <{}> '.format
                                     (id_['database'] + id_['identifier']))
            locationClause = "({})".format(" || ".join(locations))
            filters.append(locationClause)
            locationClause = "?l  faldo:location ?location .\n"

        if isinstance(drug, dict):
            drugs = []
            for id_ in drug['ids']:
                    drugs.append('?drug = <{}> '.format
                                 (id_['database'] + id_['identifier']))
            drugsClause = "({})".format(" || ".join(drugs))

            filters.append(drugsClause)

        if isinstance(disease, dict):
            diseases = []
            for id_ in disease['ids']:
                    diseases.append('?disease = <{}> '.format
                                    (id_['database'] + id_['identifier']))
            diseasesClause = "({})".format(" || ".join(diseases))
            filters.append(diseasesClause)

        filter = "FILTER ({})".format(' && '.join(filters))
        query = query.replace("%FILTER%", filter)
        query = query.replace("%LOCATION%", locationClause)

        # TODO this is just static why not put in straight?
        query = query.replace("%PROPERTIES%", "".join(
            ["distinct ?s ?location ?location_label ",
             "?disease ?disease_label ?drug ?drug_label ",
             "?evidence ?evidence_label"]))

        # TODO why is this commented out?
        # query += ("LIMIT {} OFFSET {} ".format(pageSize, offset))
        return query

    def _detailQuery(self, subject=''):
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
        # annotation keys
        locationKey = 'location'
        hasObject = 'http://purl.org/oban/association_has_object'
        has_approval_status = \
            'http://purl.obolibrary.org/obo/RO_has_approval_status'
        # location keys
        GENO_0000408 = 'http://purl.obolibrary.org/obo/GENO_0000408'

        # Build feature element
        location = annotation[locationKey]
        if GENO_0000408 in location:
            id_, ontologySource = self.namespaceSplit(
                                        location[GENO_0000408]['val'])
            name = location[GENO_0000408]['label']
        else:
            id_, ontologySource = self.namespaceSplit(location['id'])
            name = location['id']

        featureTypeTerm = protocol.OntologyTerm().fromJsonDict({
                "term": name,
                "id": id_,
                "sourceName": ontologySource})

        feature = protocol.Feature().fromJsonDict({
                "id": annotation['id'],
                "featureType": featureTypeTerm.toJsonDict()})

        # TODO ontology term backing should be switchable
        # not hardcoded in the g2p turtle file?

        # Build phenotype element
        phenotypeId, phenotypeOntologySource = self.namespaceSplit(
                                       annotation[hasObject]['val'])

        phenotypeOntologyTerm = protocol.OntologyTerm().fromJsonDict({
            "term": annotation[hasObject]['label'],
            "sourceName": phenotypeOntologySource})

        phenotypeInstance = protocol.PhenotypeInstance().fromJsonDict({
            "id": phenotypeId,
            "type": phenotypeOntologyTerm.toJsonDict()})

        # Build evidence element
        #  ECO or OBI is recommended
        if has_approval_status in annotation:
            approval_status = annotation[has_approval_status]
            id_, ontologySource = self.namespaceSplit(approval_status['val'])
            evidenceTypeOntologyTerm = protocol.OntologyTerm().fromJsonDict({
                "sourceName": ontologySource,
                "id": id_})
            evidence = protocol.Evidence().fromJsonDict({
                "evidenceType": evidenceTypeOntologyTerm.toJsonDict()})
            evidence.evidenceType.term = ''
            if 'label' in approval_status:
                evidence.evidenceType.term = approval_status['label']
                # TODO why are we throwing request validation here?
                # is this a data quality error?
                if not protocol.Evidence.validate(evidence.toJsonDict()):
                    raise exceptions.RequestValidationFailureException(
                        evidence.toJsonDict(), protocol.Evidence)
        # TODO what about multiple features? sequence annotation

        fpa = protocol.FeaturePhenotypeAssociation().fromJsonDict({
            "id": annotation['id'],
            "features": [feature.toJsonDict()],
            "phenotypeAssociationSetId": self.getId(),
            "evidence": [evidence.toJsonDict()],
            "phenotype": phenotypeInstance.toJsonDict(),
            "environmentalContexts": []})
        return fpa

    def namespaceSplit(self, url, separator='/'):
        """
        given a url return the id of the resource and the ontology source
        """
        url = url.strip("<").strip(">")
        o = urlparse.urlparse(url)
        id_ = o.path.split(separator)[-1]
        ontologySource = urlparse.urlunsplit([o[0],
                                              o[1],
                                              o[2].replace(id_, ''),
                                              o[4], ''])
        return id_, ontologySource
