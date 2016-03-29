"""
Biodata objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import os
import json

import ga4gh.datamodel as datamodel
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


def readJsonMetadata(filepath):
    jsonFilename = os.path.splitext(filepath)[0] + ".json"
    metadata = {}
    try:
        with open(jsonFilename) as data:
            metadata = json.load(data)
    except (ValueError, IOError):
        raise exceptions.FileOpenFailedException(jsonFilename)
    return metadata


class AbstractBioSample(datamodel.DatamodelObject):
    """
    This class represents an abstract BioSample object.
    It sets default values and getters, as well as the
    toProtocolElement function.
    """
    compoundIdClass = datamodel.BioSampleCompoundId

    def __init__(self, parentContainer, localId, individualName=None):
        super(AbstractBioSample, self).__init__(parentContainer, localId)
        self._createDateTime = datetime.datetime.now().isoformat()
        self._updateDateTime = datetime.datetime.now().isoformat()
        self._description = ""
        self._disease = None
        self._info = {}
        self._name = localId
        if individualName is not None:
            datasetId = self.getParentContainer().getCompoundId()
            compoundId = datamodel.IndividualCompoundId(
                datasetId, individualName)
            self._individualId = str(compoundId)
        else:
            self._individualId = None

    def toProtocolElement(self):
        bioSample = protocol.BioSample()
        bioSample.createDateTime = self.getCreateDateTime()
        bioSample.description = self.getDescription()
        bioSample.disease = self.getDisease()
        bioSample.id = self.getId()
        bioSample.individualId = self.getIndividualId()
        bioSample.info = self.getInfo()
        bioSample.name = self.getName()
        bioSample.updateDateTime = self.getUpdateDateTime()
        return bioSample

    def getIndividualId(self):
        return self._individualId

    def getCreateDateTime(self):
        return self._createDateTime

    def getUpdateDateTime(self):
        return self._updateDateTime

    def getDescription(self):
        return self._description

    def getDisease(self):
        return self._disease

    def getInfo(self):
        return self._info

    def getName(self):
        return self._name


class JsonBioSample(AbstractBioSample):
    """
    Allows a BioSample object to be populated using JSON
    metadata information by passing a filepath during
    instantiation.
    b = JsonBioSample(dataset, name, filepath)
    """
    def __init__(self, parentContainer, localId, filePath):
        super(JsonBioSample, self).__init__(parentContainer, localId)
        self._filePath = filePath
        metadata = readJsonMetadata(filePath)
        self._createDateTime = metadata.get(
            "createDateTime", self._createDateTime)
        self._updateDateTime = metadata.get(
            "updateDateTime", self._updateDateTime)
        self._description = metadata.get("description", self._description)
        self._disease = metadata.get("disease", self._disease)
        self._info = metadata.get("info", self._info)
        self._name = metadata.get("name", self._name)
        self._individualId = metadata.get("individualId", self._individualId)
        self.validate()

    def validate(self):
        if not protocol.BioSample.validate(
                self.toProtocolElement().toJsonDict()):
            raise exceptions.FileOpenFailedException(self._filePath)


class AbstractIndividual(datamodel.DatamodelObject):
    """
    This class represents an abstract Individual object.
    It sets default values and getters, as well as the
    toProtocolElement function.
    """
    compoundIdClass = datamodel.IndividualCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractIndividual, self).__init__(parentContainer, localId)
        self._createDateTime = datetime.datetime.now().isoformat()
        self._updateDateTime = datetime.datetime.now().isoformat()
        self._description = ""
        self._species = None
        self._sex = None
        self._info = {}
        self._name = localId

    def toProtocolElement(self):
        gaIndividual = protocol.Individual()
        gaIndividual.createDateTime = self.getCreateDateTime()
        gaIndividual.updateDateTime = self.getUpdateDateTime()
        gaIndividual.description = self.getDescription()
        gaIndividual.species = self.getSpecies()
        gaIndividual.sex = self.getSex()
        gaIndividual.id = self.getId()
        gaIndividual.info = self.getInfo()
        gaIndividual.name = self.getName()
        return gaIndividual

    def getCreateDateTime(self):
        return self._createDateTime

    def getUpdateDateTime(self):
        return self._updateDateTime

    def getDescription(self):
        return self._description

    def getSpecies(self):
        return self._species

    def getSex(self):
        return self._sex

    def getInfo(self):
        return self._info

    def getName(self):
        return self._name


class JsonIndividual(AbstractIndividual):
    """
    Allows a Individual object to be populated using JSON
    metadata information by passing a filepath during
    instantiation.
    b = JsonIndividual(dataset, name, filepath)
    """
    def __init__(self, parentContainer, localId, filePath):
        super(JsonIndividual, self).__init__(parentContainer, localId)
        self._filePath = filePath
        metadata = readJsonMetadata(filePath)
        self._createDateTime = metadata.get(
            "createDateTime", self._createDateTime)
        self._updateDateTime = metadata.get(
            "updateDateTime", self._updateDateTime)
        self._description = metadata.get(
            "description", self._description)
        self._species = metadata.get(
            "species", self._species)
        self._sex = metadata.get("sex", self._sex)
        self._info = metadata.get("info", self._info)
        self._name = metadata.get("name", self._name)
        self.validate()

    def validate(self):
        if not protocol.Individual.validate(
                self.toProtocolElement().toJsonDict()):
            raise exceptions.FileOpenFailedException(self._filePath)

    def getSpecies(self):
        if self._species is not None:
            return protocol.OntologyTerm().fromJsonDict(
                self._species)
        else:
            return None

    def getSex(self):
        if self._sex is not None:
            return protocol.OntologyTerm().fromJsonDict(
                self._sex)
        return None
