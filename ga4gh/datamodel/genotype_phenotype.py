"""
Module responsible for translating g2p data into GA4GH native
objects.
"""

import rdflib
import urlparse
import sets
import os
import types
import sys
import traceback

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
            self._rdfGraph.parse(format='n3', data=self.testData())
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

# flake8: noqa
    def testData(self):
        return """
@prefix ns1: <http://purl.org/oban/> .
@prefix ns2: <http://purl.obolibrary.org/obo/> .
@prefix ns3: <http://purl.org/dc/elements/1.1/> .
@prefix ns4: <http://biohackathon.org/resource/faldo#> .
@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix xml: <http://www.w3.org/XML/1998/namespace> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .
@prefix owl: <http://www.w3.org/2002/07/owl#> .
@prefix OBAN: <http://purl.org/oban/> .
@prefix OBO: <http://purl.obolibrary.org/obo/> .
@prefix DrugBank: <http://www.drugbank.ca/drugs/> .
@prefix CGD: <http://ohsu.edu/cgd/> .
@prefix faldo: <http://biohackathon.org/resource/faldo#> .

<http://ohsu.edu/cgd/0531e5c5> a OBAN:association ;
    OBO:RO_has_approval_status <http://ohsu.edu/cgd/489e03d2> ;
    OBAN:association_has_object <http://ohsu.edu/cgd/1b6442aa> ;
    OBAN:association_has_object_property OBO:RO_0002606 ;
    OBAN:association_has_subject DrugBank:DB00619 .

CGD:bed32b89 a OBAN:association ;
    OBO:RO_has_approval_status CGD:dc7d284a ;
    OBAN:association_has_object <http://ohsu.edu/cgd/9a0c7b6e> ;
    OBAN:association_has_object_property OBO:RO_0002606 ;
    OBAN:association_has_subject DrugBank:DB01254 .

CGD:c2ed4da0 a OBAN:association ;
    OBO:RO_has_approval_status CGD:c703f7ab ;
    OBAN:association_has_object <http://ohsu.edu/cgd/37da8697> ;
    OBAN:association_has_object_property OBO:RO_0002606 ;
    OBAN:association_has_subject DrugBank:DB00398 .

CGD:fae52a14 a OBAN:association ;
    OBO:RO_has_approval_status <http://ohsu.edu/cgd/489e03d2> ;
    OBAN:association_has_object CGD:bff18fe6 ;
    OBAN:association_has_object_property OBO:RO_0002606 ;
    OBAN:association_has_subject DrugBank:DB01268 .

<http://ohsu.edu/cgd/489e03d2> a owl:Class ;
    rdfs:label "late trials" .

CGD:dc7d284a a owl:Class ;
    rdfs:label "preclinical" .

CGD:c703f7ab a owl:Class ;
    rdfs:label "early trials" .

<http://ohsu.edu/cgd/1b6442aa> a OBO:OMIM_606764 ;
    rdfs:label "GIST with biomarker KIT  wild type no mutation" ;
    OBO:RO_has_biomarker <http://ohsu.edu/cgd/27d2169c> .

<http://ohsu.edu/cgd/9a0c7b6e> a OBO:OMIM_606764 ;
    rdfs:label "GIST with biomarker KIT  wild type no mutation" ;
    OBO:RO_has_biomarker CGD:d212cbf9 .

<http://ohsu.edu/cgd/37da8697> a OBO:OMIM_606764 ;
    rdfs:label "GIST with biomarker KIT  wild type no mutation" ;
    OBO:RO_has_biomarker <http://ohsu.edu/cgd/4841bf74> .

CGD:bff18fe6 a OBO:OMIM_606764 ;
    rdfs:label "GIST with biomarker KIT  wild type no mutation" ;
    OBO:RO_has_biomarker <http://ohsu.edu/cgd/24228ad6> .

<http://ohsu.edu/cgd/27d2169c> a OBO:SO_0001059 ;
    rdfs:label "KIT  wild type no mutation" ;
    OBO:RO_0002200 <http://ohsu.edu/cgd/87752f6c> .

<http://ohsu.edu/cgd/24228ad6> a OBO:SO_0001059 ;
    rdfs:label "KIT  wild type no mutation" ;
    OBO:RO_0002200 <http://ohsu.edu/cgd/87795e43> .

<http://ohsu.edu/cgd/87752f6c> a OBO:OMIM_606764 ;
    rdfs:label "GIST with decreased sensitivity to therapy" ;
    OBO:BFO_0000159 CGD:decreased_sensitivity .

<http://ohsu.edu/cgd/5c895709> a OBO:OMIM_606764 ;
    rdfs:label "GIST with sensitivity to therapy" ;
    OBO:BFO_0000159 CGD:sensitivity .

<http://ohsu.edu/cgd/87795e43> a OBO:OMIM_606764 ;
    rdfs:label "GIST with response to therapy" ;
    OBO:BFO_0000159 CGD:response .

OBO:OMIM_606764 a owl:Class ;
    rdfs:label "GIST" ;
    rdfs:subClassOf OBO:DOID_4 .

CGD:d212cbf9 a OBO:SO_0001059 ;
    rdfs:label "KIT  wild type no mutation" ;
    OBO:RO_0002200 <http://ohsu.edu/cgd/5c895709> .

<http://ohsu.edu/cgd/4841bf74> a OBO:SO_0001059 ;
    rdfs:label "KIT  wild type no mutation" ;
    OBO:RO_0002200 <http://ohsu.edu/cgd/505c855a> .

DrugBank:DB00619 a owl:Class ;
    rdfs:label "imatinib" ;
    OBO:RO_0002606 <http://ohsu.edu/cgd/0348822e>,
        <http://ohsu.edu/cgd/0663d7aa>,
        <http://ohsu.edu/cgd/10027e2b>,
        <http://ohsu.edu/cgd/1141037f>,
        <http://ohsu.edu/cgd/16ed1e21>,
        <http://ohsu.edu/cgd/1b49fa54>,
        <http://ohsu.edu/cgd/1b6442aa>,
        <http://ohsu.edu/cgd/3660373b>,
        <http://ohsu.edu/cgd/396e17ca>,
        <http://ohsu.edu/cgd/4387ad11>,
        <http://ohsu.edu/cgd/49773dd9>,
        <http://ohsu.edu/cgd/4c9df69d>,
        <http://ohsu.edu/cgd/593ba90d>,
        <http://ohsu.edu/cgd/61a0ae6e>,
        <http://ohsu.edu/cgd/6200b434>,
        <http://ohsu.edu/cgd/622f5d87>,
        <http://ohsu.edu/cgd/66956de7>,
        <http://ohsu.edu/cgd/6f842dcb>,
        <http://ohsu.edu/cgd/6fe97f66>,
        <http://ohsu.edu/cgd/782a5adf>,
        <http://ohsu.edu/cgd/80857771>,
        <http://ohsu.edu/cgd/829ecfe4>,
        <http://ohsu.edu/cgd/83589dc1>,
        <http://ohsu.edu/cgd/8d3ef6b9>,
        CGD:a05b065c,
        CGD:a28d6b63,
        CGD:a52d56bd,
        CGD:ad7fe4a2,
        CGD:b0523fef,
        CGD:b40a93d7,
        CGD:b8059d5f,
        CGD:b82bb601,
        CGD:bd5ec05d,
        CGD:c3f8dbe4,
        CGD:c585206a,
        CGD:c5e74da6,
        CGD:cd634e7c,
        CGD:cfa27e95,
        CGD:cfaf14e7,
        CGD:d467afaa,
        CGD:d634d636,
        CGD:daa8ce89,
        CGD:de0c523d,
        CGD:e4e2d39b,
        CGD:ebd2e0a0,
        CGD:ed09a007,
        CGD:f1c0bea0,
        CGD:f52055a9,
        CGD:f5c5ae17,
        CGD:fe0f3b26 ;
    rdfs:subClassOf OBO:CHEBI_23888 .

DrugBank:DB01254 a owl:Class ;
    rdfs:label "dasatinib" ;
    OBO:RO_0002606 <http://ohsu.edu/cgd/0270d49e>,
        <http://ohsu.edu/cgd/06d1537a>,
        <http://ohsu.edu/cgd/07be695d>,
        <http://ohsu.edu/cgd/204c5c4d>,
        <http://ohsu.edu/cgd/211a6bfd>,
        <http://ohsu.edu/cgd/24fc6515>,
        <http://ohsu.edu/cgd/312a7c01>,
        <http://ohsu.edu/cgd/316ba70b>,
        <http://ohsu.edu/cgd/37ba4794>,
        <http://ohsu.edu/cgd/39b33923>,
        <http://ohsu.edu/cgd/4214ffb0>,
        <http://ohsu.edu/cgd/53477124>,
        <http://ohsu.edu/cgd/5730ad5d>,
        <http://ohsu.edu/cgd/5d496096>,
        <http://ohsu.edu/cgd/60d258f3>,
        <http://ohsu.edu/cgd/6339e561>,
        <http://ohsu.edu/cgd/79c33766>,
        <http://ohsu.edu/cgd/7b37c3d3>,
        <http://ohsu.edu/cgd/7f7b6467>,
        <http://ohsu.edu/cgd/80412f7f>,
        <http://ohsu.edu/cgd/83c53f39>,
        <http://ohsu.edu/cgd/85722ca9>,
        <http://ohsu.edu/cgd/88afcbf8>,
        <http://ohsu.edu/cgd/90988c75>,
        <http://ohsu.edu/cgd/9a0c7b6e>,
        <http://ohsu.edu/cgd/9a3b8f86>,
        <http://ohsu.edu/cgd/9ec35520>,
        CGD:a1313702,
        CGD:a60218cf,
        CGD:b0148c82,
        CGD:b0a03037,
        CGD:b12958da,
        CGD:b5bacf6d,
        CGD:c602b6f4,
        CGD:caf00fdf,
        CGD:cc04d9d5,
        CGD:d006d388,
        CGD:df43c04f,
        CGD:ef8944b8,
        CGD:f84840a9 ;
    rdfs:subClassOf OBO:CHEBI_23888 .

DrugBank:DB00398 a owl:Class ;
    rdfs:label "sorafenib" ;
    OBO:RO_0002606 <http://ohsu.edu/cgd/1165021d>,
        <http://ohsu.edu/cgd/37da8697>,
        <http://ohsu.edu/cgd/47b7a3ba>,
        <http://ohsu.edu/cgd/7aca8eca>,
        <http://ohsu.edu/cgd/7b12ee2a>,
        <http://ohsu.edu/cgd/87cabd66>,
        <http://ohsu.edu/cgd/8b978c84>,
        <http://ohsu.edu/cgd/8f795083>,
        <http://ohsu.edu/cgd/90fadb4d>,
        CGD:a7b840b5,
        CGD:acaf1e33,
        CGD:af2fb811,
        CGD:b8293b8b,
        CGD:bb66b1e9,
        CGD:bf53fae5,
        CGD:c56941f8,
        CGD:d90f6f8c,
        CGD:f21e69ef,
        CGD:f4b538b0 ;
    rdfs:subClassOf OBO:CHEBI_23888 .

DrugBank:DB01268 a owl:Class ;
    rdfs:label "sunitinib" ;
    OBO:RO_0002606 <http://ohsu.edu/cgd/11ba3f14>,
        <http://ohsu.edu/cgd/202aab8b>,
        <http://ohsu.edu/cgd/330fdb9d>,
        <http://ohsu.edu/cgd/392fb86a>,
        <http://ohsu.edu/cgd/3f99c173>,
        <http://ohsu.edu/cgd/4da27469>,
        <http://ohsu.edu/cgd/54039374>,
        <http://ohsu.edu/cgd/5c4e4aa5>,
        <http://ohsu.edu/cgd/5ea8d37b>,
        <http://ohsu.edu/cgd/65b688b1>,
        <http://ohsu.edu/cgd/71fe9f0f>,
        <http://ohsu.edu/cgd/7e7ac65b>,
        <http://ohsu.edu/cgd/86415d40>,
        <http://ohsu.edu/cgd/88e4c778>,
        <http://ohsu.edu/cgd/983a1528>,
        CGD:bff18fe6,
        CGD:c9ce40f4,
        CGD:d8f69729,
        CGD:d9a4077f,
        CGD:e901ee1e,
        CGD:f7303413 ;
    rdfs:subClassOf OBO:CHEBI_23888 .

CGD:d8c2d551 a OBO:SO_0001059,
        OBO:SO_0001583,
        owl:NamedIndividual ;
    rdfs:label "KIT V559I, H687Y, T670, V654A, A829P, D816, N822, Y823D missense mutation" ;
    faldo:location <http://www.monarchinitiative.org/_654654UniProtKB:P10721#P10721-1Region> ;
    OBO:GENO_0000408 <http://www.ncbi.nlm.nih.gov/gene/3815> ;
    OBO:GENO_reference_amino_acid "V" ;
    OBO:GENO_results_in_amino_acid_change "A" ;
    OBO:RO_0002200 <http://ohsu.edu/cgd/4c7ff2c7> ;
    OBO:RO_0002205 <http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=3496.1> .

"""


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
        source = 'http://purl.org/dc/elements/1.1/source'
        location = 'location'
        hasObject = 'http://purl.org/oban/association_has_object'
        has_approval_status = 'http://purl.obolibrary.org/obo/RO_has_approval_status'
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
