"""
Module responsible for translating g2p data into GA4GH native
objects.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import rdflib

import ga4gh.datamodel as datamodel
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol

# annotation keys
TYPE = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type'
LABEL = 'http://www.w3.org/2000/01/rdf-schema#label'
HAS_QUALITY = 'http://purl.obolibrary.org/obo/BFO_0000159'


class AbstractPhenotypeAssociationSet(datamodel.DatamodelObject):
    compoundIdClass = datamodel.PhenotypeAssociationSetCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractPhenotypeAssociationSet, self).__init__(
            parentContainer, localId)

    def toProtocolElement(self):
        pas = protocol.PhenotypeAssociationSet()
        pas.name = self.getLocalId()
        pas.id = self.getId()
        pas.dataset_id = self.getParentContainer().getId()
        return pas


class SimulatedPhenotypeAssociationSet(AbstractPhenotypeAssociationSet):
    def __init__(self, parentContainer, localId, randomSeed):
        super(SimulatedPhenotypeAssociationSet, self).__init__(
            parentContainer, localId)

    def getAssociations(
            self, request=None,
            pageSize=None, offset=0):
        if request:
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


class G2PUtility(object):

    def _featureTypeLabel(self, featureType):
        """
        return a label for known types
        """
        featureTypes = {'http://purl.obolibrary.org/obo/SO_0001059':
                        'sequence_alteration',
                        'http://purl.obolibrary.org/obo/SO_0001583':
                        'missense_variant',
                        'http://purl.obolibrary.org/obo/SO_0000147': 'exon'}
        if featureType in featureTypes:
            return featureTypes[featureType]
        return featureType

    def _extractAssociationsDetails(self, associations):
        """
        given a set of results from our search query, return a
        list of the URIRef for each `detail` (feature,environment,phenotype)
        """
        detailedURIRef = []
        for row in associations.bindings:
            # empty set [{}]
            if 'feature' in row:
                detailedURIRef.append(row['feature'])
                detailedURIRef.append(row['environment'])
                detailedURIRef.append(row['phenotype'])

        return detailedURIRef

    def _detailTuples(self, uriRefs):
        """
        given a list of uriRefs, return a list of dicts:
        {'subject': s, 'predicate': p, 'object': o }
        all values are strings
        """
        details = []
        for uriRef in uriRefs:
            for s, p, o in self._rdfGraph.triples((uriRef, None, None)):
                details.append({
                    'subject': s.toPython(),
                    'predicate': p.toPython(),
                    'object': o.toPython()
                })
        return details

    def _bindingsToDict(self, bindings):
        """
        given a binding from the sparql query result,
        create a dict of plain text
        """
        myDict = {}
        for key, val in bindings.iteritems():
            myDict[key.toPython().replace('?', '')] = val.toPython()
        return myDict

    def _addDataFile(self, filename):
        """ given a filename, add it to the graph """
        if filename.endswith('.ttl'):
            self._rdfGraph.parse(filename, format='n3')
        else:
            self._rdfGraph.parse(filename, format='xml')

    def _getDetails(self, uriRef, associations_details):
        """
        given a uriRef, return a dict of all the details for that Ref
        use the uriRef as the 'id' of the dict
        """
        associationDetail = {}
        for detail in associations_details:
            if detail['subject'] == uriRef:
                associationDetail[detail['predicate']] = detail['object']
            associationDetail['id'] = uriRef
        return associationDetail

    def _formatExternalIdentifiers(self, element, element_type):
        """
        Formats several external identifiers for query
        """
        elementClause = None
        elements = []
        if not issubclass(element.__class__, dict):
            element = protocol.toJsonDict(element)
        if element['externalIdentifiers']:
            for _id in element['externalIdentifiers']:
                elements.append(self._formatExternalIdentifier(_id,
                                                               element_type))
            elementClause = "({})".format(" || ".join(elements))
        return elementClause

    def _formatExternalIdentifier(self, element, element_type):
        """
        Formats a single external identifier for query
        """
        if "http" not in element['database']:
            term = "{}:{}".format(element['database'], element['identifier'])
            namespaceTerm = self._toNamespaceURL(term)
        else:
            namespaceTerm = "{}{}".format(element['database'],
                                          element['identifier'])
        comparison = '?{} = <{}> '.format(element_type, namespaceTerm)
        return comparison

    def _formatOntologyTerm(self, element, element_type):
        """
        Formats the ontology terms for query
        """
        elementClause = None
        if isinstance(element, dict) and element.get('terms'):
            elements = []
            for _term in element['terms']:
                if _term.get('id'):
                    elements.append('?{} = <{}> '.format(
                        element_type, _term['id']))
                else:
                    elements.append('?{} = <{}> '.format(
                        element_type, self._toNamespaceURL(_term['term'])))
            elementClause = "({})".format(" || ".join(elements))
        return elementClause

    def _formatOntologyTermObject(self, terms, element_type):
        """
        Formats the ontology term object for query
        """
        elementClause = None
        if not isinstance(terms, collections.Iterable):
            terms = [terms]
        elements = []
        for term in terms:
            if not issubclass(term.__class__, dict):
                term = protocol.toJsonDict(term)
            if term['id']:
                elements.append('?{} = <{}> '.format(
                    element_type, term['id']))
            else:
                elements.append('?{} = <{}> '.format(
                    element_type, self._toNamespaceURL(term['term'])))
        elementClause = "({})".format(" || ".join(elements))
        return elementClause

    def _formatFeatureClause(self, feature):
        """
        Format a feature for lookup by gene_symbol
        """
        return 'regex(?feature_label, "{}")'.format(feature.gene_symbol)

    def _formatId(self, element, element_type):
        """
        Formats an identifiers for query
        """
        return '?{} = <{}> '.format(
            element_type, element)

    def _formatIds(self, element, element_type):
        """
        Formats a set of identifiers for query
        """
        elementClause = None
        if isinstance(element, collections.Iterable):
            elements = []
            for _id in element:
                elements.append('?{} = <{}> '.format(
                    element_type, _id))
            elementClause = "({})".format(" || ".join(elements))
        return elementClause

    def _formatEvidence(self, elements):
        elementClause = None
        filters = []
        for evidence in elements:
            if evidence.description:
                elementClause = 'regex(?{}, "{}")'.format(
                    'environment_label', evidence.description)
            if hasattr(evidence, 'externalIdentifiers') \
                    and evidence.externalIdentifiers:
                # TODO will this pick up > 1 externalIdentifiers ?
                for externalIdentifier in evidence['externalIdentifiers']:
                    exid_clause = self._formatExternalIdentifier(
                        externalIdentifier, 'environment')
                    # cleanup parens from _formatExternalIdentifier method
                    elementClause = exid_clause[1:-1]
            if elementClause:
                filters.append(elementClause)
        elementClause = "({})".format(" || ".join(filters))
        return elementClause

    def _baseQuery(self):
        return """
            PREFIX OBAN: <http://purl.org/oban/>
            PREFIX OBO: <http://purl.obolibrary.org/obo/>
            PREFIX dc: <http://purl.org/dc/elements/1.1/>
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            PREFIX BFO: <http://purl.obolibrary.org/obo/BFO_>
            PREFIX owl: <http://www.w3.org/2002/07/owl#>
            PREFIX faldo: <http://biohackathon.org/resource/faldo#>
            SELECT
                ?association
                ?environment
                ?environment_label
                ?feature
                ?feature_label
                ?phenotype
                ?phenotype_label
                (GROUP_CONCAT(?source; separator="|") AS ?sources)
                ?evidence_type
                ?external_id
                ?phenotype_quality
                WHERE {
                    ?association  a OBAN:association .
                    ?association    OBO:RO_0002558 ?evidence_type .
                    ?association    OBO:RO_has_environment ?environment   .
                    OPTIONAL { ?association  dc:source ?source } .
                    ?association    OBAN:association_has_subject ?feature .
                    ?association    OBAN:association_has_object ?phenotype .
                    ?environment  rdfs:label ?environment_label  .
                    ?phenotype  rdfs:label ?phenotype_label  .
                    ?feature  rdfs:label ?feature_label  .
                    OPTIONAL { ?feature owl:sameAs  ?external_id } .
                    OPTIONAL { ?phenotype BFO:0000159 ?phenotype_quality } .
                    #%FILTER%
                    }
            GROUP  BY ?association
            ORDER  BY ?association
"""

    def _toNamespaceURL(self, term):
        """
        Given an ontologyterm.term return namespace identifier.
        Leverages prefixes already in graph namespace
        Ex.  "DrugBank:DB01268" -> "http://www.drugbank.ca/drugs/DB01268"
        """
        (termPrefix, termId) = term.split(':')
        for prefix, namespace in self._rdfGraph.namespaces():
            if prefix == termPrefix:
                return "".join([namespace, termId])
        raise exceptions.NotImplementedException(
           "Term has a prefix not found in this instance. {}".format(term))

    def _getIdentifier(self, url):
        """
        Given a url identifier return identifier portion
        Leverages prefixes already in graph namespace
        Returns None if no match
        Ex.  "http://www.drugbank.ca/drugs/DB01268" -> "DB01268"
        """
        for prefix, namespace in self._rdfGraph.namespaces():
            if namespace in url:
                return(url.replace(namespace, ''))

    def _getPrefix(self, url):
        """
        Given a url return namespace prefix.
        Leverages prefixes already in graph namespace
        Ex.  "http://www.drugbank.ca/drugs/" -> "Drugbank"
        """
        for prefix, namespace in self._rdfGraph.namespaces():
            if namespace.toPython() == url or namespace == url:
                return prefix
        raise exceptions.NotImplementedException(
           "No namespace found for url {}".format(url))

    def _getPrefixURL(self, url):
        """
        Given a url return namespace prefix.
        Leverages prefixes already in graph namespace
        Ex.  "http://www.drugbank.ca/drugs/DDD"
            -> "http://www.drugbank.ca/drugs/"
        """
        for prefix, namespace in self._rdfGraph.namespaces():
            if namespace.toPython() in url:
                return(namespace)

    def _toGA4GH(self, association):
        """
        given an association dict,
        return a protocol.FeaturePhenotypeAssociation
        """
        # TODO: This method needs to broken in several parts. It's
        # doing too much things
        fpa = None

        # The association dict has the keys: environment, environment
        # label, evidence, feature label, phenotype and sources. Each
        # key's value is a dict with the RDF predicates as keys and
        # subject as values

        # useful
        # ECO_0000033 traceable author statement
        # RO_0002558 has evidence
        # RO_0002200 has phenotype
        # RO_0002606 is substance that treats
        # SO_0001059 sequence_alteration
        # BFO_0000159 has quality
        # OMIM_606764
        # OBO:SO_0000147 exon

        feature = association['feature']

        fpa = protocol.FeaturePhenotypeAssociation()
        fpa.id = association['id']
        # fpa.features.extend([f])
        fpa.feature_ids.extend([feature['id']])

        msg = 'Association: genotype:[{}] phenotype:[{}] environment:[{}] ' \
              'evidence:[{}] publications:[{}]'
        fpa.description = msg.format(
            association['feature_label'],
            association['phenotype_label'],
            association['environment_label'],
            self._getIdentifier(association['evidence']),
            association['sources']
            )
        evidence = protocol.Evidence()
        phenotype = association['phenotype']

        term = protocol.OntologyTerm()
        term.term = association['evidence_type']
        term.id = phenotype['id']
        term.source_version = self._version
        term.source_name = self._getPrefix(
            self._getPrefixURL(association['id']))
        evidence.evidence_type.MergeFrom(term)

        evidence.description = self._getIdentifier(association['evidence'])
        # TODO there is nowhere in evidence to place list of sources?
        fpa.evidence.extend([evidence])

        # map environment (drug)
        environmentalContext = protocol.EnvironmentalContext()
        environment = association['environment']
        environmentalContext.id = environment['id']
        environmentalContext.description = association['environment_label']

        term = protocol.OntologyTerm()
        term.term = environment['id']
        term.id = 'http://purl.obolibrary.org/obo/RO_0002606'
        term.source_version = self._version
        term.source_name = self._getPrefix(
            self._getPrefixURL(association['id']))
        environmentalContext.environment_type.MergeFrom(term)

        fpa.environmental_contexts.extend([environmentalContext])

        phenotypeInstance = protocol.PhenotypeInstance()
        term = protocol.OntologyTerm()
        term.term = phenotype[TYPE]
        term.id = phenotype['id']
        term.source_version = self._version
        term.source_name = self._getPrefix(
            self._getPrefixURL(association['id']))
        phenotypeInstance.type.MergeFrom(term)

        phenotypeInstance.description = phenotype[LABEL]
        phenotypeInstance.id = phenotype['id']
        fpa.phenotype.MergeFrom(phenotypeInstance)
        fpa.phenotype_association_set_id = self.getId()
        return fpa


class PhenotypeAssociationSet(G2PUtility, AbstractPhenotypeAssociationSet):
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
        # save the path
        self._dataUrl = dataDir

        try:
            self._scanDataFiles(dataDir, ['*.ttl', '*.xml'])
        except AttributeError:
            pass

        # extract version
        cgdTTL = rdflib.URIRef("http://data.monarchinitiative.org/ttl/cgd.ttl")
        versionInfo = rdflib.URIRef(
            u'http://www.w3.org/2002/07/owl#versionInfo')
        self._version = None
        for s, p, o in self._rdfGraph.triples((cgdTTL, versionInfo, None)):
            self._version = o.toPython()

    def getAssociations(self, request=None, pageSize=None, offset=0,
                        featureSets=[]):
        """
        This query is the main search mechanism.
        It queries the graph for annotations that match the
        AND of [feature,environment,phenotype].
        """
        # query to do search
        query = self._formatFilterQuery(request, featureSets)
        associations = self._rdfGraph.query(query)
        # associations is now a dict with rdflib terms with variable and
        # URIrefs or literals

        # given get the details for the feature,phenotype and environment
        associations_details = self._detailTuples(
                                    self._extractAssociationsDetails(
                                        associations))

        # association_details is now a list of {subject,predicate,object}
        # for each of the association detail
        # http://nmrml.org/cv/v1.0.rc1/doc/doc/objectproperties/BFO0000159___-324347567.html
        # label "has quality at all times" (en)
        associationList = []
        for assoc in associations.bindings:
            if '?feature' in assoc:
                association = self._bindingsToDict(assoc)
                association['feature'] = self._getDetails(
                    association['feature'],
                    associations_details)
                association['environment'] = self._getDetails(
                    association['environment'],
                    associations_details)
                association['phenotype'] = self._getDetails(
                    association['phenotype'],
                    associations_details)
                association['evidence'] = association['phenotype'][HAS_QUALITY]
                association['id'] = association['association']
                associationList.append(association)

                # our association list is now a list of dicts with the
                # elements of an association: environment, evidence,
                # feature, phenotype and sources, each with their
                # references (labels or URIrefs)

        # create GA4GH objects
        associations = [self._toGA4GH(assoc) for
                        assoc in associationList]
        return associations

    def _formatFilterQuery(self, request=None, featureSets=[]):
        """
        Generate a formatted sparql query with appropriate filters
        """
        query = self._baseQuery()

        filters = []

        if issubclass(request.__class__,
                      protocol.SearchGenotypePhenotypeRequest):
            filters += self._filterSearchGenotypePhenotypeRequest(request,
                                                                  featureSets)

        if issubclass(request.__class__, protocol.SearchGenotypesRequest):
            filters += self._filterSearchGenotypesRequest(request)

        if issubclass(request.__class__, protocol.SearchPhenotypesRequest):
            filters += self._filterSearchPhenotypesRequest(request)

        # apply filters
        filter = "FILTER ({})".format(' && '.join(filters))
        if len(filters) == 0:
            filter = ""
        query = query.replace("#%FILTER%", filter)

        return query

    def _filterSearchGenotypePhenotypeRequest(self, request, featureSets):
        filters = []
        if request.feature_ids:
            featureFilters = []
            for featureId in request.feature_ids:
                for featureSet in featureSets:
                    try:
                        compoundId = datamodel.FeatureCompoundId. \
                                                parse(featureId)
                        # we have a compoundId, so use it to lookup
                        # import ipdb; ipdb.set_trace()
                        if compoundId.feature_set == self.getLocalId():
                            featureFilters.append(self._formatId(
                                                  compoundId.featureId,
                                                  'feature'))
                            break
                        else:
                            feature = featureSet.getFeature(compoundId)
                            if feature:
                                featureFilters.append(self.
                                                      _formatFeatureClause(
                                                        feature))
                                break
                    except Exception as e:
                        featureFilters.append(self._formatId(
                                              "NO-FIND", 'feature'))
            if len(featureFilters) > 0:
                filters.append("({})".format(" || ".join(featureFilters)))

        # if request.genotype_ids:
        #     featureClause = self._formatIds(request.genotype_ids,
        #                                     'feature')
        #     if featureClause:
        #         filters.append(featureClause)

        if request.evidence:
            evidenceClause = self._formatEvidence(request.evidence)
            if evidenceClause:
                filters.append(evidenceClause)

        if request.phenotype_ids:
            phenotypeClause = self._formatIds(request.phenotype_ids,
                                              'phenotype')
            filters.append(phenotypeClause)

        return filters

    def _filterSearchGenotypesRequest(self, request):
        """
        Filters the request for genotype search requests
        """
        filters = []
        if request.id:
            filters.append("?feature = <{}>".format(request.id))

        if request.reference_name:
            filters.append('regex(?feature_label, "{}")'
                           .format(request.reference_name))

        featureClause = self._formatExternalIdentifiers(request, 'external_id')
        if featureClause:
            filters.append(featureClause)
        # OntologyTerms
        featureOntologytermsClause = self._formatOntologyTerm(request,
                                                              'feature')
        if featureOntologytermsClause:
            filters.append(featureOntologytermsClause)

        return filters

    def _filterSearchPhenotypesRequest(self, request):
        """
        Filters request for phenotype search requests
        """
        filters = []
        if request.id:
            filters.append("?phenotype = <{}>".format(request.id))

        if request.description:
            filters.append('regex(?phenotype_label, "{}")'
                           .format(request.description))
        # OntologyTerms
        # TODO: refactor this repetitive code
        if hasattr(request.type, 'id') and request.type.id:
            ontolgytermsClause = self._formatOntologyTermObject(
                request.type, 'phenotype')
            if ontolgytermsClause:
                filters.append(ontolgytermsClause)
        if len(request.qualifiers) > 0:
            ontolgytermsClause = self._formatOntologyTermObject(
                request.qualifiers, 'phenotype_quality')
            if ontolgytermsClause:
                filters.append(ontolgytermsClause)
        if hasattr(request.age_of_on_set, 'id') and request.age_of_on_set.id:
            ontolgytermsClause = self._formatOntologyTermObject(
                request.age_of_on_set, 'phenotype_quality')
            if ontolgytermsClause:
                filters.append(ontolgytermsClause)
        return filters
