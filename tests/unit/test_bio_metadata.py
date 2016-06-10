"""
Tests the biodata module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel.datasets as datasets
import ga4gh.exceptions as exceptions
import ga4gh.datamodel.bio_metadata as bioMetadata
import ga4gh.protocol as protocol
from ga4gh import pb


class TestIndividuals(unittest.TestCase):
    """
    Tests the Individuals class
    """
    def testToProtocolElement(self):
        dataset = datasets.Dataset('dataset1')
        term = protocol.OntologyTerm()
        term.term = "male genotypic sex"
        term.id = "PATO:0020001"
        term.source_name = "PATO"
        term.source_version = pb.string("2015-11-18")
        # Write out a valid input
        print(protocol.toJsonDict(term))
        validIndividual = protocol.Individual(
            name="test",
            created="2016-05-19T21:00:19Z",
            updated="2016-05-19T21:00:19Z",
            sex=term)
        validIndividual.info['test'].values.add().string_value = 'test-info'
        # pass through protocol creation
        individual = bioMetadata.Individual(
            dataset, "test")
        individual.populateFromJson(protocol.toJson(validIndividual))
        gaIndividual = individual.toProtocolElement()
        # Verify elements exist
        self.assertEqual(gaIndividual.created, validIndividual.created)
        self.assertEqual(gaIndividual.updated, validIndividual.updated)
        # Invalid input
        invalidIndividual = '{"bad:", "json"}'
        individual = bioMetadata.Individual(dataset, "test")
        # Should fail
        self.assertRaises(
            exceptions.InvalidJsonException,
            individual.populateFromJson,
            invalidIndividual)


class TestBioSamples(unittest.TestCase):
    """
    Tests the BioSamples class
    """
    def testToProtocolElement(self):
        dataset = datasets.Dataset('dataset1')
        # Write out a valid input
        validBioSample = protocol.BioSample(
            name="test",
            created="2016-05-19T21:00:19Z",
            updated="2016-05-19T21:00:19Z")
        validBioSample.info['test'].values.add().string_value = 'test-info'
        # pass through protocol creation
        bioSample = bioMetadata.BioSample(
            dataset, "test")
        bioSample.populateFromJson(protocol.toJson(validBioSample))
        gaBioSample = bioSample.toProtocolElement()
        # Verify elements exist
        self.assertEqual(gaBioSample.created, validBioSample.created)
        self.assertEqual(gaBioSample.updated, validBioSample.updated)
        # Invalid input
        invalidBioSample = '{"bad:", "json"}'
        bioSample = bioMetadata.Individual(dataset, "test")
        # Should fail
        self.assertRaises(
            exceptions.InvalidJsonException,
            bioSample.populateFromJson,
            invalidBioSample)
