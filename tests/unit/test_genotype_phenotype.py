"""
Unit tests for genotypephenotype objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel.genotype_phenotype as g2p
import ga4gh.datamodel.datasets as datasets


class TestPhenotypeAssociationSet(unittest.TestCase):
    def setUp(self):
        ds = datasets.AbstractDataset("test")
        self.phenotypeAssocationSet = g2p.PhenotypeAssociationSet(
            ds, "test", None)

    def testFormatQuery(self):
        #  "At least one of [location, drug, disease] must be specified"
        query = u'\n        PREFIX OBAN: <http://purl.org/oban/>\n' \
                u'        PREFIX OBO: <http://purl.obolibrary.org/obo/>\n' \
                u'        PREFIX rdfs: ' \
                u'<http://www.w3.org/2000/01/rdf-schema#>\n' \
                u'        PREFIX faldo: ' \
                u'<http://biohackathon.org/resource/faldo#>\n' \
                u'        PREFIX rdf: ' \
                u'<http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n' \
                u'        SELECT  distinct ' \
                u'?s ?location ?location_label ?disease ?' \
                u'disease_label ?drug ?drug_label ' \
                u'?evidence ?evidence_label\n' \
                u'            WHERE {\n' \
                u'                ?s  a OBAN:association .\n' \
                u'                ?s  OBAN:association_has_subject ?drug .\n' \
                u'                ?s  ' \
                u'OBAN:association_has_object  ?disease .\n' \
                u'                ?s  ' \
                u'OBO:RO_has_approval_status ?evidence .\n' \
                u'                ?disease ' \
                u'OBO:RO_has_biomarker ?location .\n' \
                u'                ?drug  ' \
                u'rdfs:label ?drug_label  .\n' \
                u'                ?disease  ' \
                u'rdfs:label ?disease_label  .\n' \
                u'                ?disease ' \
                u'rdf:type ?disease_type .\n' \
                u'                ?evidence  ' \
                u'rdfs:label ?evidence_label  .\n' \
                u'                ?location  ' \
                u'rdfs:label ?location_label  .\n' \
                u'                ?disease_type  ' \
                u'rdfs:label ?disease_type_label  .\n' \
                u'                FILTER ' \
                u'(regex(?location_label, "*LOCATION*") ' \
                u'&& regex(?drug_label, "***DRUG***") ' \
                u'&& regex(?disease_label, "***DISEASE***"))\n' \
                u'        }\n        ORDER BY ?s\n            '
        self.assertEqual(
            self.phenotypeAssocationSet._formatFilterQuery(
                location="*LOCATION*",
                drug="***DRUG***",
                disease="***DISEASE***"), query)

    def testNamespaceSplit(self):
        # I think this is broken
        given = "<http://ohsu.edu/cgd/bff18fe6>"
        expect = ("bff18fe6", "http://ohsu.edu/cgd/")
        self.assertEqual(
            self.phenotypeAssocationSet.namespaceSplit(given), expect)
        given = "http://ohsu.edu/cgd/bff18fe6"
        expect = ("bff18fe6", "http://ohsu.edu/cgd/")
        self.assertEqual(
            self.phenotypeAssocationSet.namespaceSplit(given), expect)

    def testToGA4GH(self):
        pass

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
        self.assertEqual(self.phenotypeAssocationSet._flatten(d), expect)
