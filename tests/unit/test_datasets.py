"""
Tests the datasets module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel.datasets as datasets


class TestDatasets(unittest.TestCase):
    """
    Tests the datasets class
    """
    def testToProtocolElement(self):
        dataset_id = 'ds1'
        dataset = datasets.SimulatedDataset(dataset_id, 1, 2, 3, 4)
        gaDataset = dataset.toProtocolElement()
        self.assertEqual(dataset.getId(), gaDataset.id)
