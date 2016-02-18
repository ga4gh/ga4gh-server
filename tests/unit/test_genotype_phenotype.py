"""
Unit tests for genotypephenotype objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel.genotype_phenotype as g2p
import ga4gh.exceptions as exceptions


class TestG2PDataset(unittest.TestCase):
    def setUp(self):
        self.g2pDataset = g2p.G2PDataset()

    def testSearchQuery(self):
        # This is hardcoded during __init__
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
            """  # watch for trailing blank space
        self.assertEqual(self.g2pDataset._searchQuery, query)

    def testFormatQuery(self):
        #  "At least one of [location, drug, disease] must be specified"
        self.assertRaises(
            exceptions.NotImplementedException,
            self.g2pDataset._formatQuery)
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
                FILTER (regex(?location_label, "*************"))
        }
        ORDER BY ?s
            """
        self.assertEqual(
            self.g2pDataset._formatQuery(
                location="*************"), query)
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
                FILTER (regex(?location_label, "*LOCATION*")"""
        # because of flake :(
        query += """ && regex(?drug_label, "***DRUG***")"""
        query += """ && regex(?disease_label, "***DISEASE***"))
        }
        ORDER BY ?s
            """
        self.assertEqual(
            self.g2pDataset._formatQuery(
                location="*LOCATION*",
                drug="***DRUG***",
                disease="***DISEASE***"), query)


class TestAnnotationFactory(unittest.TestCase):
    def setUp(self):
        self.annFac = g2p.AnnotationFactory({}, [])

    def testNamespaceSplit(self):
        # I think this is broken
        given = "<http://ohsu.edu/cgd/bff18fe6>"
        expect = ("bff18fe6", "http://ohsu.edu/cgd/")
        self.assertEqual(self.annFac.namespaceSplit(given), expect)
        given = "http://ohsu.edu/cgd/bff18fe6"
        expect = ("bff18fe6", "http://ohsu.edu/cgd/")
        self.assertEqual(self.annFac.namespaceSplit(given), expect)

    def testToGA4GH(self):
        # print(self.annFac._toGA4GH(
        # {"location": "***LOCATION***", "id": "***ID***"}))
        pass

    def testAnnotations(self):
        # No annotations are loaded should be empty
        self.assertEqual(self.annFac.annotations(), [])

    def testFlatten(self):
        # TODO which value of `s` should become the id?
        d = [
            {"p": "predicate1",
             "s": "subject1",
             "o": "object1",
             "label": "label1"
             },
            {"p": "predicate2",
             "s": "subject2",
             "o": "object2",
             "label": "label2"
             },
            {"p": "predicate1",
             "s": "subject3",
             "o": "object3",
             "label": "label3"
             }
        ]

        expect = {u'predicate2':
                  {u'val': u'object2', u'label': u'label2'},
                  u'predicate1': [{u'val': u'object1', u'label': u'label1'},
                                  {u'val': u'object3', u'label': u'label3'}],
                  u'id': u'subject3'}
        self.assertEqual(self.annFac._flatten(d), expect)
