"""
Tests the datasets module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
from datetime import datetime

import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.exceptions as exceptions


class TestDatasets(unittest.TestCase):
    """
    Tests the datasets class
    """
    def testToProtocolElement(self):
        datasetId = 'ds1'
        dataset = datasets.SimulatedDataset(datasetId, 1, 2, 3, 4, 5)
        dataset.setInfo({"test": "test"})
        createDT = '2015-01-01T00:00:00'  # set in simulated
        updateDT = datetime.utcnow().isoformat()
        dataset.setUpdateDateTime(updateDT)
        gaDataset = dataset.toProtocolElement()
        self.assertIsNotNone(gaDataset.info)
        self.assertEqual(gaDataset.info['test'].values[0].string_value, "test")
        self.assertEqual(dataset.getId(), gaDataset.id)
        self.assertEqual(gaDataset.createDateTime, createDT)
        self.assertEqual(gaDataset.updateDateTime, updateDT)

    def testBadDateTime(self):
        dataset = datasets.SimulatedDataset('ds1', 1, 2, 3, 4, 5)
        badDatetime = '215010100:00:00'
        self.assertRaises(
                exceptions.TimestampFormatException,
                dataset.setUpdateDateTime, badDatetime)
