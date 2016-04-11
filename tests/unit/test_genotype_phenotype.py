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

# from ga4gh.protocol import FeaturePhenotypeAssociation


class TestPhenotypeAssociationSet(unittest.TestCase):
    def setUp(self):
        ds = datasets.AbstractDataset("test")
        self.phenotypeAssocationSet = g2p.PhenotypeAssociationSet(
            ds, "test", None)

    def testFormatQuery(self):
        """
        At least one of [location, environment, phenotype] must be specified
        """
        query = """
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
                    #%FILTER%
                    }
            GROUP  BY ?association
            ORDER  BY ?association
"""
        query = query.replace("#%FILTER%",
                              'FILTER (regex(?feature_label, "*LOCATION*") ' +
                              '&& regex(?environment_label, "***DRUG***") ' +
                              '&& regex(?phenotype_label, "***DISEASE***"))'
                              )
        # see all differences
        self.maxDiff = None
        self.assertEqual(
            self.phenotypeAssocationSet._formatFilterQuery(
                feature="*LOCATION*",
                environment="***DRUG***",
                phenotype="***DISEASE***"), query)
