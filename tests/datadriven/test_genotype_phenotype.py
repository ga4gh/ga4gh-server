"""
Data-driven tests for g2p.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import rdflib

import ga4gh.server.datamodel.genotype_phenotype as genotype_phenotype
import ga4gh.server.datamodel.datasets as datasets
import tests.datadriven as datadriven
import tests.paths as paths

import ga4gh.schemas.protocol as protocol


def testG2P():
    testDataDir = os.path.join(
        paths.testDataDir, 'datasets/dataset1/phenotypes')
    for test in datadriven.makeTests(testDataDir, PhenotypeAssociationSetTest):
        yield test


class PhenotypeAssociationSetTest(datadriven.DataDrivenTest):
    def __init__(self, localId, baseDir):
        self._dataset = datasets.Dataset("ds")
        super(PhenotypeAssociationSetTest, self).__init__(localId, baseDir)
        self.phenotypeAssocationSet = self.getDataModelInstance(
            localId, baseDir)

    def getDataModelInstance(self, localId, dataPath):
        return genotype_phenotype.RdfPhenotypeAssociationSet(
            self._dataset, localId, dataPath)

    def getProtocolClass(self):
        return protocol.PhenotypeAssociationSet

    def testDetailTuples(self):
        test_uriRefs = [
            rdflib.term.URIRef(u'http://ohsu.edu/cgd/27d2169c'),
            rdflib.term.URIRef(u'http://ohsu.edu/cgd/87752f6c')
        ]
        details = self.phenotypeAssocationSet._detailTuples(test_uriRefs)
        self.assertEqual(len(details), 6)
        sample1 = {
            u'predicate': u'http://www.w3.org/2000/01/rdf-schema#label',
            u'object': u'GIST with decreased sensitivity to therapy',
            u'subject': u'http://ohsu.edu/cgd/87752f6c',
        }
        sample2 = {
            u'predicate': u'http://purl.obolibrary.org/obo/RO_0002200',
            u'object': u'http://ohsu.edu/cgd/87752f6c',
            u'subject': u'http://ohsu.edu/cgd/27d2169c',
        }
        self.assertIn(sample1, details)
        self.assertIn(sample2, details)

    def testBindingsToDict(self):
        bindings = {
            rdflib.term.Variable(u'association'):
            rdflib.term.URIRef(u'http://ohsu.edu/cgd/fe484b5c'),
            rdflib.term.Variable(u'feature'):
            rdflib.term.URIRef(
                u'http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=965'
            ),
            rdflib.term.Variable(u'phenotype_label'):
            rdflib.term.Literal(
                u'Papillary thyroid carcinoma with sensitivity to therapy'
            ),
            rdflib.term.Variable(u'environment_label'):
            rdflib.term.Literal(u'sunitinib'),
            rdflib.term.Variable(u'environment'):
            rdflib.term.URIRef(u'http://www.drugbank.ca/drugs/DB01268'),
            rdflib.term.Variable(u'evidence_type'):
            rdflib.term.URIRef(u'http://purl.obolibrary.org/obo/ECO_0000033'),
            rdflib.term.Variable(u'sources'):
            rdflib.term.Literal(
                u'http://www.ncbi.nlm.nih.gov/pubmed/21470995|'
                'http://www.ncbi.nlm.nih.gov/pubmed/21470995'),
            rdflib.term.Variable(u'phenotype'):
            rdflib.term.URIRef(u'http://ohsu.edu/cgd/30ebfd1a'),
            rdflib.term.Variable(u'feature_label'):
            rdflib.term.Literal(u'COSM965'),
        }
        myDict = self.phenotypeAssocationSet._bindingsToDict(bindings)
        sampleDict = {
            u'environment_label': u'sunitinib',
            u'feature_label': u'COSM965',
            u'evidence_type': u'http://purl.obolibrary.org/obo/ECO_0000033',
            u'feature':
            u'http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=965',
            u'environment': u'http://www.drugbank.ca/drugs/DB01268',
            u'sources':
            u'http://www.ncbi.nlm.nih.gov/pubmed/21470995|'
            'http://www.ncbi.nlm.nih.gov/pubmed/21470995',
            u'phenotype': u'http://ohsu.edu/cgd/30ebfd1a',
            u'phenotype_label':
            u'Papillary thyroid carcinoma with sensitivity to therapy',
            u'association': u'http://ohsu.edu/cgd/fe484b5c',
        }
        self.assertEqual(myDict, sampleDict)

    def testGetDetails(self):
        uriRef = 'http://www.drugbank.ca/drugs/DB01268'
        associations_details = [
            {u'predicate': u'http://purl.obolibrary.org/obo/RO_0002606',
             u'object': u'http://ohsu.edu/cgd/54039374',
             u'subject': u'http://www.drugbank.ca/drugs/DB01268'},
            {u'predicate': u'http://purl.obolibrary.org/obo/RO_0002606',
             u'object': u'http://ohsu.edu/cgd/983a1528',
             u'subject': u'http://www.drugbank.ca/drugs/DB01268'},
            {u'predicate': u'http://purl.obolibrary.org/obo/RO_0002606',
             u'object': u'http://ohsu.edu/cgd/71fe9f0f',
             u'subject': u'http://www.drugbank.ca/drugs/DB01268'},
            {u'predicate': u'http://www.w3.org/2000/01/rdf-schema#subClassOf',
             u'object': u'http://purl.obolibrary.org/obo/CHEBI_23888',
             u'subject': u'http://www.drugbank.ca/drugs/DB01268'}
        ]
        sample_details = {
            u'http://purl.obolibrary.org/obo/RO_0002606':
            u'http://ohsu.edu/cgd/71fe9f0f',
            u'http://www.w3.org/2000/01/rdf-schema#subClassOf':
            u'http://purl.obolibrary.org/obo/CHEBI_23888',
            u'id': u'http://www.drugbank.ca/drugs/DB01268'}
        details = self.phenotypeAssocationSet._getDetails(
            uriRef, associations_details)
        self.assertEqual(details, sample_details)

    def testToNamespaceURL(self):
        sample_term = 'DrugBank:DB01268'
        result = self.phenotypeAssocationSet._toNamespaceURL(sample_term)
        self.assertEqual('http://www.drugbank.ca/drugs/DB01268', result)

    def testGetIdentifier(self):
        sample_url = 'http://www.drugbank.ca/drugs/DB01268'
        result = self.phenotypeAssocationSet._getIdentifier(sample_url)
        self.assertEqual('DB01268', result)

    def testGetPrefix(self):
        sample_url = 'http://www.drugbank.ca/drugs/'
        result = self.phenotypeAssocationSet._getPrefix(sample_url)
        self.assertEqual('DrugBank', result)

    def testGetPrefixURL(self):
        sample_url = 'http://www.drugbank.ca/drugs/DDD'
        result = self.phenotypeAssocationSet._getPrefixURL(sample_url)
        self.assertEqual('http://www.drugbank.ca/drugs/', str(result))

    def testExtractAssociationsDetails(self):
        sample_query = """
        PREFIX OBAN: <http://purl.org/oban/>
            PREFIX OBO: <http://purl.obolibrary.org/obo/>
            PREFIX dc: <http://purl.org/dc/elements/1.1/>
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
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
                    FILTER ((?feature = <http://ohsu.edu/cgd/27d2169c> ))
                    }
            GROUP  BY ?association
            ORDER  BY ?association
        """
        sample_associations = \
            self.phenotypeAssocationSet._rdfGraph.query(sample_query)
        result = self.phenotypeAssocationSet._extractAssociationsDetails(
            sample_associations)
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0].toPython(), 'http://ohsu.edu/cgd/27d2169c')
        self.assertEqual(
            result[1].toPython(),
            'http://www.drugbank.ca/drugs/DB00619')
        self.assertEqual(result[2].toPython(), 'http://ohsu.edu/cgd/87752f6c')

    def testToGA4GH(self):
        sample_associations = {
            u'environment_label': u'sunitinib',
            u'feature_label': u'RET M918T missense mutation',
            u'evidence_type': u'http://purl.obolibrary.org/obo/ECO_0000033',
            u'feature': {
                u'http://purl.obolibrary.org/obo/GENO_0000408':
                u'http://www.ncbi.nlm.nih.gov/gene/5979',
                u'http://purl.obolibrary.org/obo/GENO_reference_amino_acid':
                u'M',
                u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
                u'http://purl.obolibrary.org/obo/SO_0001059',
                u'http://biohackathon.org/resource/faldo#location':
                u'http://www.monarchinitiative.org/_918918UniProtKB:'
                'P07949#P07949-1Region',
                u'http://purl.obolibrary.org/obo/GENO_reference_nucleotide':
                u'T',
                u'http://purl.obolibrary.org/obo/'
                'GENO_results_in_amino_acid_change':
                u'T',
                u'http://purl.obolibrary.org/obo/RO_0002200':
                u'http://ohsu.edu/cgd/3774b1d2',
                u'http://purl.obolibrary.org/obo/RO_0002205':
                u'http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?'
                'REQUEST=CCDS&DATA=7200.1',
                u'http://purl.obolibrary.org/obo/GENO_altered_nucleotide':
                u'C',
                u'http://www.w3.org/2000/01/rdf-schema#label':
                u'RET M918T missense mutation',
                u'id': u'http://cancer.sanger.ac.uk/cosmic/mutation/'
                'overview?id=965',
                u'http://www.w3.org/2002/07/owl#sameAs':
                u'http://www.ncbi.nlm.nih.gov/SNP/74799832',
            },
            u'evidence': u'http://ohsu.edu/cgd/sensitivity',
            u'environment': {
                u'http://purl.obolibrary.org/obo/RO_0002606':
                u'http://ohsu.edu/cgd/71fe9f0f',
                u'http://www.w3.org/2000/01/rdf-schema#subClassOf':
                u'http://purl.obolibrary.org/obo/CHEBI_23888',
                u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
                u'http://www.w3.org/2002/07/owl#Class',
                u'http://www.w3.org/2000/01/rdf-schema#label': u'sunitinib',
                u'id': u'http://www.drugbank.ca/drugs/DB01268',
            },
            u'sources':
            u'http://www.ncbi.nlm.nih.gov/pubmed/21470995|'
            'http://www.ncbi.nlm.nih.gov/pubmed/21470995',
            u'phenotype': {
                u'http://purl.obolibrary.org/obo/BFO_0000159':
                u'http://ohsu.edu/cgd/sensitivity',
                u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
                u'http://purl.obolibrary.org/obo/DOID_3969',
                u'http://www.w3.org/2000/01/rdf-schema#label':
                u'Papillary thyroid carcinoma with sensitivity to therapy',
                u'id': u'http://ohsu.edu/cgd/30ebfd1a',
            },
            u'phenotype_label':
            u'Papillary thyroid carcinoma with sensitivity to therapy',
            u'id': u'http://ohsu.edu/cgd/fe484b5c',
            u'association': u'http://ohsu.edu/cgd/fe484b5c',
        }
        result = self.phenotypeAssocationSet._toGA4GH(sample_associations)
        self.assertEqual(
            result.__class__.__name__, 'FeaturePhenotypeAssociation')
        fpa_dict = protocol.toJsonDict(result)
        description = 'Association: genotype:[RET M918T missense mutation]' \
                      ' phenotype:[Papillary thyroid carcinoma with ' \
                      'sensitivity to therapy] environment:[sunitinib]' \
                      ' evidence:[sensitivity] publications:' \
                      '[http://www.ncbi.nlm.nih.gov/pubmed/21470995|' \
                      'http://www.ncbi.nlm.nih.gov/pubmed/21470995]'

        self.assertEqual(fpa_dict['description'], description)
        self.assertIn('featureIds', fpa_dict.keys())
        self.assertIn('evidence', fpa_dict.keys())
        self.assertIn('environmentalContexts', fpa_dict.keys())
        self.assertEqual(len(fpa_dict['featureIds']), 1)
        self.assertEqual(len(fpa_dict['evidence']), 1)
        self.assertEqual(len(fpa_dict['environmentalContexts']), 1)
