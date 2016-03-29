"""
Tests the biodata module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import json
import os

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets
import ga4gh.exceptions as exceptions
import ga4gh.datamodel.biodata as biodata
import ga4gh.protocol as protocol
import ga4gh.avrotools as avrotools

import tests.paths as paths


class TestIndividuals(unittest.TestCase):
    """
    Tests the JSON Individuals class
    """
    dataset = datasets.AbstractDataset('dataset1')

    def testToProtocolElement(self):
        for localId in os.listdir(paths.testIndividualsDataDir):
            jsonFilename = os.path.join(paths.individualsDir, localId)
            jsonDict = {}
            try:
                with open(jsonFilename) as data:
                    jsonDict = json.load(data)
            except (ValueError, IOError):
                # Poorly formed JSON throws expected exception
                self.assertRaises(
                    exceptions.FileOpenFailedException,
                    biodata.JsonIndividual,
                    self.dataset, localId, jsonFilename)
            if jsonDict != {}:
                jsonDict['id'] = str(
                    datamodel.IndividualCompoundId(
                        self.dataset._compoundId, localId))
                validator = avrotools.Validator(protocol.Individual)
                if validator.getInvalidFields(jsonDict) != {}:
                    self.assertRaises(
                        exceptions.FileOpenFailedException,
                        biodata.JsonIndividual,
                        self.dataset,
                        localId, jsonFilename)
                else:
                    self.assertTrue(protocol.Individual.validate(jsonDict))
                    individual = biodata.JsonIndividual(
                        self.dataset, localId, jsonFilename)
                    gaIndividual = individual.toProtocolElement()
                    for key in jsonDict:
                        self.assertEqual(
                            jsonDict[key],
                            gaIndividual.__getattribute__(key))


class TestBioSamples(unittest.TestCase):
    """
    Tests the JSON BioSamples class
    """
    dataset = datasets.AbstractDataset('dataset1')

    def testToProtocolElement(self):
        for localId in os.listdir(paths.testBioSamplesDataDir):
            jsonFilename = os.path.join(paths.bioSamplesDir, localId)
            jsonDict = {}
            try:
                with open(jsonFilename) as data:
                    jsonDict = json.load(data)
            except (ValueError, IOError):
                # Poorly formed JSON throws expected exception
                self.assertRaises(
                    exceptions.FileOpenFailedException,
                    biodata.JsonBioSample, self.dataset, localId, jsonFilename)
            if jsonDict != {}:
                jsonDict['id'] = str(
                    datamodel.BioSampleCompoundId(
                        self.dataset._compoundId, localId))
                validator = avrotools.Validator(protocol.BioSample)
                if validator.getInvalidFields(jsonDict) != {}:
                    self.assertRaises(
                        exceptions.FileOpenFailedException,
                        biodata.JsonBioSample,
                        self.dataset,
                        localId, jsonFilename)
                else:
                    self.assertTrue(protocol.BioSample.validate(jsonDict))
                    bioSample = biodata.JsonBioSample(
                        self.dataset, localId, jsonFilename)
                    gaBioSample = bioSample.toProtocolElement()
                    for key in jsonDict:
                        self.assertEqual(
                            jsonDict[key],
                            gaBioSample.__getattribute__(key))
