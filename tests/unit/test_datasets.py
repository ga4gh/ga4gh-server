"""
Tests the datasets module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.server.datamodel.datasets as datasets


class TestDatasets(unittest.TestCase):
    """
    Tests the datasets class
    """
    def testToProtocolElement(self):
        datasetId = 'ds1'
        dataset = datasets.SimulatedDataset(datasetId, 1, 2, 3, 4, 5)
        dataset.setAttributes({"test": "test"})
        gaDataset = dataset.toProtocolElement()
        self.assertIsNotNone(gaDataset.attributes.attr)
        self.assertEqual(
            gaDataset.attributes.attr['test'].values[0].string_value, "test")
        self.assertEqual(dataset.getId(), gaDataset.id)
