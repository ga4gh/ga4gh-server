"""
Module responsible for translating g2p data into GA4GH native
objects.
"""

# no futures
# std lib
import rdflib
import urlparse
import time
# locals
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


class G2PDataset:
    """
    An rdf object store.  The cancer genome database
    [Clinical Genomics Knowledge Base]
    (http://nif-crawler.neuinfo.org/monarch/ttl/cgd.ttl)
    published by the Monarch project was the source of Evidence.
    """

    def __init__(self, sources):
        """
        Initialize dataset, using the passed dict of sources
        [{source,format}] see rdflib.parse() for more
        """

        self._searchQuery = """
            PREFIX OBAN: <http://purl.org/oban/>
            PREFIX OBO: <http://purl.obolibrary.org/obo/>
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            PREFIX faldo: <http://biohackathon.org/resource/faldo#>
            PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            SELECT %PROPERTIES%
                WHERE {
                    ?s    a OBAN:association .
                    ?s  OBAN:association_has_subject ?l .
                    ?l rdfs:label ?location_label  .
                    %LOCATION%
                    ?s  OBO:RO_has_environment  ?drug .
                    ?drug  rdfs:label ?drug_label  .
                    ?s  OBAN:association_has_object  ?d .
                    ?d  rdfs:label ?disease_label  .
                    ?d rdf:type ?disease .
                    ?s  OBAN:association_has_object_property  ?evidence .
                    OPTIONAL {  ?evidence  rdfs:label ?evidence_label } .
              %FILTER%
            }
            """

        # initialize graph
        self._rdfGraph = rdflib.ConjunctiveGraph()

        # load with data
        for source in sources:
            if not source['format']:
                self._rdfGraph.parse(source['source'])
            else:
                self._rdfGraph.parse(source['source'],
                                     format=source['format'])

        # TODO is this necessary?
        self.associationsLength = 0

    def _search(self, request):
        offset = request.offset

        associations = self.queryLabels(
            request.feature, request.evidence, request.phenotype,
            request.pageSize, offset)

        self.associationsLength = len(associations)
        for association in associations:
            yield association

    def queryLabels(
         self, location=None, drug=None, disease=None, pageSize=None,
         offset=0):

        """
        This query is the main search mechanism.
        It queries the graph for annotations that match the
        AND of [location,drug,disease]
        """
        query = self.formatQuery(location, drug, disease)

        query = query.replace("%PROPERTIES%",
                              "distinct ?s  ?location ?location_label " +
                              "?disease ?disease_label ?drug  ?drug_label")

        query += ("LIMIT {} OFFSET {} ".format(pageSize, offset))

        # print(query)

        results = self._rdfGraph.query(query)

        annotations = []
        for row in results:
            annotations.append(
                self.toGA4GH(
                              self.query("<%s>" %  row['s'].toPython())))

        return annotations

    def formatQuery(self, location=None, drug=None, disease=None):
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

    def query(self, subject=''):
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
          %SUBJECT% ?p ?o .
          OPTIONAL {?o  rdfs:label ?label .}
        }
        """

        annotationQuery = annotationQuery.replace("%SUBJECT%", subject)
        results = self._rdfGraph.query(annotationQuery)
        now = time.time()
        rows = [row.asdict() for row in results]
        #print('annotationQuery',time.time()-now)

        for row in rows:
            for k in row:
                row[k] = row[k].toPython()
            row['s'] = subject

        locationQuery = """
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
             PREFIX OBO: <http://purl.obolibrary.org/obo/>
        SELECT distinct *
        WHERE {
        %SUBJECT%    a  OBO:SO_0001059  .
        %SUBJECT%   ?p ?o .
        OPTIONAL {?o  rdfs:label ?label .} .
        }
        """

        locationRows = []
        for row in rows:
            if row['p'] == 'http://purl.org/oban/association_has_subject':
                location = "<" + row['o'] + ">"
                locationQuery = locationQuery.replace(
                    "%SUBJECT%",location)
                # print(locationQuery)
                now = time.time()
                results = self._rdfGraph.query(locationQuery)
                locationRows = [row.asdict() for row in results]
                for row in locationRows:
                    for k in row:
                        row[k] = row[k].toPython()
                    row['s'] = location
                #print('locationQuery',time.time()-now)

        annotation = self.flatten(rows)
        location = self.flatten(locationRows)
        annotation['location'] = location
        return annotation

    def flatten(self, dict):
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

    def toGA4GH(self, annotation):
        """
        given an annotation dict, return a protocol.FeaturePhenotypeAssociation
        """

        fpa = None

        # annotation keys
        source = 'http://purl.org/dc/elements/1.1/source'
        location = 'location'
        evidenceURI = 'http://purl.obolibrary.org/obo/RO_0002558'
        hasObject = 'http://purl.org/oban/association_has_object'
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
        if source in annotation:
            if not isinstance(annotation[source], list):
                annotation[source] = [annotation[source]]
            for src in annotation[source]:
                evidence = protocol.Evidence()
                evidence.evidenceType = protocol.OntologyTerm()
                id_, ontologySource = self.namespaceSplit(src['val'])
                evidence.evidenceType.ontologySource = ontologySource
                evidence.evidenceType.id = id_

                evidence.evidenceType.name = ''
                if 'label' in annotation[evidenceURI]:
                    evidence.evidenceType.name = \
                        annotation[evidenceURI]['label']
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
