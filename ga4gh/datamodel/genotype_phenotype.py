"""
Module responsible for translating g2p data into GA4GH native
objects.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import rdflib
import urlparse
import types

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

    def _generateFeature(self, featureId="", childIds=[], parentId="",
                         featureSetId="", referenceName="", start=None,
                         end=None, strand=None, _term="", _termId="",
                         sourceName="", sourceVersion=""):

        feature = protocol.Feature()
        feature.id = featureId
        feature.childIds = []
        feature.parentId = parentId
        feature.featureSetId = featureSetId
        feature.referenceName = referenceName
        feature.start = start
        feature.end = end
        feature.strand = strand
        feature.featureType = self._generateOntologyTerm(_term,
                                                         _termId,
                                                         sourceName,
                                                         sourceVersion)
        return feature

    def _generatePhenotypeInstance(self, _id="", _term="",
                                   _termId="", sourceName="",
                                   sourceVersion="", qualifierTerms=None,
                                   onsetTerm=None, description=""):
        phenotypeInstance = protocol.PhenotypeInstance()
        phenotypeInstance.id = _id
        phenotypeInstance.type = self._generateOntologyTerm(_term,
                                                            _termId,
                                                            sourceName,
                                                            sourceVersion)
        phenotypeInstance.qualifier = qualifierTerms
        phenotypeInstance.ageOfOnset = onsetTerm
        phenotypeInstance.description = description
        return phenotypeInstance

    def _generateOntologyTerm(self, _term="", _id="",
                              sourceName="", sourceVersion=""):
        term = protocol.OntologyTerm()
        term.term = _term
        term.id = _id
        term.sourceName = sourceName
        term.sourceVersion = sourceVersion
        return term

    def _generateEvidence(self, _term="", _id="",
                          sourceName="", sourceVersion="", description=""):
        evidence = protocol.Evidence()
        evidence.evidenceType = self._generateOntologyTerm(_term,
                                                           _id,
                                                           sourceName,
                                                           sourceVersion)
        evidence.description = description
        return evidence

    def _generateFeaturePhenotypeAssociation(self, _id="",
                                             phenotypeAssociationSetId="",
                                             features=[], evidence=[],
                                             phenotype=None, description="",
                                             environmentalContexts=[]):
        fpa = protocol.FeaturePhenotypeAssociation()
        fpa.id = _id
        fpa.phenotypeAssociationSetId = phenotypeAssociationSetId
        fpa.features = features
        fpa.evidence = evidence
        fpa.phenotype = phenotype
        fpa.description = description
        fpa.environmentalContexts = environmentalContexts
        return fpa


class SimulatedPhenotypeAssociationSet(AbstractPhenotypeAssociationSet):
    def __init__(self, parentContainer, localId, randomSeed):
        super(SimulatedPhenotypeAssociationSet, self).__init__(
            parentContainer, localId)


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

        self._scanDataFiles(dataDir, ['*.ttl', '*.xml'])

    def toProtocolElement(self):
        pas = protocol.PhenotypeAssociationSet()
        pas.name = self.getLocalId()
        pas.id = self.getId()
        pas.datasetId = self.getParentContainer().getId()
        pas.info = self.getInfo()
        return pas

    def _addDataFile(self, filename):
        if filename.endswith('.ttl'):
            self._rdfGraph.parse(filename, format='n3')
        else:
            self._rdfGraph.parse(filename, format='xml')

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
            # TODO check this regex (versus string match)
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
        # TODO this is just static why not put in straight?
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

        # TODO whys it a set? we called distinct
        uniqueAnnotations = set()
        for row in results:
            uniqueAnnotations.add("<{}>".format(row['s'].toPython()))

        annotations = AnnotationFactory(self._rdfGraph,
                                        uniqueAnnotations,
                                        self.getId()).annotations()

        return annotations


# TODO refactor into single class
class AnnotationFactory(AbstractPhenotypeAssociationSet):
    """
    Given a RDF query result set, return corresponding set of ProtocolElements
    """

    def __init__(self, graph, uniqueAnnotations, parentId):
        self._rdfGraph = graph
        self._parentId = parentId

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

        f = self._generateFeature(
            featureId=annotation['id'],  # TODO connect with real feature Ids
            childIds=[],
            parentId="",
            featureSetId="",
            referenceName="",
            start=None,
            end=None,
            strand=None,
            _term=name,
            _termId=id_,
            sourceName=ontologySource,
            sourceVersion="")

        # TODO ontology term backing should be switchable
        # not hardcoded in the g2p turtle file?

        phenotypeId, phenotypeOntologySource = self.namespaceSplit(
                                       annotation[hasObject]['val'])

        phenotypeInstance = self._generatePhenotypeInstance(
            _id=phenotypeId,
            _term=annotation[hasObject]['label'],
            _termId="",
            sourceName=phenotypeOntologySource,
            sourceVersion="",
            qualifierTerms=None,
            onsetTerm=None,
            description="")

        #  ECO or OBI is recommended
        if has_approval_status in annotation:
            approval_status = annotation[has_approval_status]
            id_, ontologySource = self.namespaceSplit(approval_status['val'])
            evidence = self._generateEvidence(
                description="",
                _id=id_,
                sourceName=ontologySource)
            evidence.evidenceType.term = ''
            if 'label' in approval_status:
                evidence.evidenceType.term = approval_status['label']
                # TODO why are we throwing request validation here?
                # is this a data quality error?
                if not protocol.Evidence.validate(evidence.toJsonDict()):
                    raise exceptions.RequestValidationFailureException(
                        evidence.toJsonDict(), protocol.Evidence)

        fpa = self._generateFeaturePhenotypeAssociation(
            _id=annotation['id'],  # TODO why is this the ID?
            phenotypeAssociationSetId=self._parentId,  # TODO remove postrefac
            features=[f],  # TODO what about multiple features?
            evidence=[evidence],
            phenotype=phenotypeInstance,
            description="",
            environmentalContexts=[])

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
