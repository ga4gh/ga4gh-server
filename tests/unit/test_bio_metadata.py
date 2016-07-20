"""
Tests the biodata module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol
import ga4gh.registry as registry


@unittest.skip("TODO Proper tests for Biometadata")
class TestIndividuals(unittest.TestCase):
    """
    Tests the Individuals class
    """
    def testToProtocolElement(self):
        term = protocol.OntologyTerm()
        term.term = "male genotypic sex"
        term.id = "PATO:0020001"
        term.source_name = "PATO"
        term.source_version = "2015-11-18"
        # Write out a valid input
        validIndividual = protocol.Individual(
            name="test",
            created="2016-05-19T21:00:19Z",
            updated="2016-05-19T21:00:19Z",
            sex=term)
        validIndividual.info['test'].values.add().string_value = 'test-info'
        # pass through protocol creation
        individual = registry.Individual()
        individual.populateFromJson(protocol.toJson(validIndividual))
        gaIndividual = individual.toProtocolElement()
        # Verify elements exist
        self.assertEqual(gaIndividual.created, validIndividual.created)
        self.assertEqual(gaIndividual.updated, validIndividual.updated)
        # Invalid input
        invalidIndividual = '{"bad:", "json"}'
        individual = registry.Individual("test")
        # Should fail
        self.assertRaises(
            exceptions.InvalidJsonException,
            individual.populateFromJson,
            invalidIndividual)


@unittest.skip("TODO Proper tests for Biometadata")
class TestBioSamples(unittest.TestCase):
    """
    Tests the BioSamples class
    """
    def testToProtocolElement(self):
        # Write out a valid input
        validBioSample = protocol.BioSample(
            name="test",
            created="2016-05-19T21:00:19Z",
            updated="2016-05-19T21:00:19Z")
        validBioSample.info['test'].values.add().string_value = 'test-info'
        # pass through protocol creation
        bioSample = registry.BioSample("test")
        bioSample.populateFromJson(protocol.toJson(validBioSample))
        gaBioSample = bioSample.toProtocolElement()
        # Verify elements exist
        self.assertEqual(gaBioSample.created, validBioSample.created)
        self.assertEqual(gaBioSample.updated, validBioSample.updated)
        # Invalid input
        invalidBioSample = '{"bad:", "json"}'
        bioSample = registry.Individual("test")
        # Should fail
        self.assertRaises(
            exceptions.InvalidJsonException,
            bioSample.populateFromJson,
            invalidBioSample)
