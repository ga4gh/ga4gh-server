"""
Tests the biodata module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import json

import ga4gh.datamodel.datasets as datasets
import ga4gh.exceptions as exceptions
import ga4gh.datamodel.biodata as biodata
import ga4gh.protocol as protocol


class TestIndividuals(unittest.TestCase):
    """
    Tests the Individuals class
    """
    def testToProtocolElement(self):
        dataset = datasets.AbstractDataset('dataset1')
        # Write out a valid input
        validIndividual = protocol.Individual(
            name="test",
            created="2016-05-19T21:00:19Z",
            updated="2016-05-19T21:00:19Z")
        # pass through protocol creation
        individual = biodata.Individual("test")
        individual.populateFromJson(validIndividual.toJsonDict())
        # Verify elements exist
        self.assertEqual(individual.created, validIndividual.created)
        self.assertEqual(individual.updated, validIndividual.updated)
        # Invalid input
        invalidIndividual = json.dumps({"name": "test"})
        individual = biodata.Individual(dataset, "test")
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
        dataset = datasets.AbstractDataset('dataset1')
        # Write out a valid input
        validBioSample = protocol.BioSample(
            name="test",
            created="2016-05-19T21:00:19Z",
            updated="2016-05-19T21:00:19Z")
        # pass through protocol creation
        bioSample = biodata.bioSample(dataset, "test")
        bioSample.populateFromJson(validBioSample.toJsonDict())
        # Verify elements exist
        self.assertEqual(bioSample.created, validBioSample.created)
        self.assertEqual(bioSample.updated, validBioSample.updated)
        # Invalid input
        invalidIndividual = json.dumps({"name": "test"})
        individual = biodata.Individual("test")
        # Should fail
        self.assertRaises(
            exceptions.InvalidJsonException,
            individual.populateFromJson,
            invalidIndividual)
