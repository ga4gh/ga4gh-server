"""
Biodata objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import json

import ga4gh.server.datamodel as datamodel
import ga4gh.server.exceptions as exceptions

import ga4gh.schemas.protocol as protocol


class Biosample(datamodel.DatamodelObject):
    """
    This class represents an abstract Biosample object.
    It sets default values and getters, as well as the
    toProtocolElement function.
    """
    compoundIdClass = datamodel.BiosampleCompoundId

    def __init__(self, parentContainer, localId):
        super(Biosample, self).__init__(parentContainer, localId)
        self._created = datetime.datetime.now().isoformat()
        self._updated = datetime.datetime.now().isoformat()
        self._description = None
        self._disease = None
        self._name = localId
        self._individualId = None
        self._datasetId = parentContainer.getId()
        self._individualAgeAtCollection = None

    def toProtocolElement(self):
        disease = None
        if self.getDisease():
            disease = protocol.fromJson(
                json.dumps(self.getDisease()), protocol.OntologyTerm)
        individualAgeAtCollection = None
        if self.getIndividualAgeAtCollection():
            individualAgeAtCollection = protocol.fromJson(
                json.dumps(self.getIndividualAgeAtCollection()), protocol.Age)
        biosample = protocol.Biosample(
            dataset_id=self._datasetId,
            created=self.getCreated(),
            updated=self.getUpdated(),
            description=self.getDescription(),
            id=self.getId(),
            individual_id=self.getIndividualId(),
            name=self.getName(),
            disease=disease,
            individual_age_at_collection=individualAgeAtCollection)
        self.serializeAttributes(biosample)
        return biosample

    def populateFromJson(self, jsonString):
        try:
            parsed = protocol.fromJson(jsonString, protocol.Biosample)
        except:
            raise exceptions.InvalidJsonException(jsonString)
        self._created = parsed.created
        self._updated = parsed.updated
        self._description = parsed.description
        self._disease = protocol.toJsonDict(parsed.disease)
        self._individualId = parsed.individual_id
        self._individualAgeAtCollection = protocol.toJsonDict(
                                           parsed.individual_age_at_collection)
        attributes = {}
        for key in parsed.attributes.attr:
            attributes[key] = {
                "values": protocol.toJsonDict(parsed.attributes.attr[key])}
        self.setAttributes(attributes)
        return self

    def populateFromRow(self, biosampleRecord):
        # TODO coerce to types
        self._created = biosampleRecord.created
        self._updated = biosampleRecord.updated
        self._description = biosampleRecord.description
        self._disease = json.loads(biosampleRecord.disease)
        self._individualId = biosampleRecord.individualid
        self.setAttributesJson(biosampleRecord.attributes)
        self._individualAgeAtCollection = json.loads(
                                biosampleRecord.individualAgeAtCollection)
        return self

    def setIndividualId(self, individualId):
        self._individualId = individualId

    def getIndividualId(self):
        return self._individualId

    def getCreated(self):
        return self._created

    def getUpdated(self):
        return self._updated

    def getDescription(self):
        return self._description

    def getDisease(self):
        if self._disease is not {}:
            return self._disease
        else:
            return None

    def getName(self):
        return self._name

    def getIndividualAgeAtCollection(self):
        if self._individualAgeAtCollection is not {}:
            return self._individualAgeAtCollection
        else:
            return None


class Individual(datamodel.DatamodelObject):
    """
    This class represents an abstract Individual object.
    It sets default values and getters, as well as the
    toProtocolElement function.
    """
    compoundIdClass = datamodel.IndividualCompoundId

    def __init__(self, parentContainer, localId):
        super(Individual, self).__init__(parentContainer, localId)
        self._created = datetime.datetime.now().isoformat()
        self._updated = datetime.datetime.now().isoformat()
        self._description = None
        self._species = None
        self._sex = None
        self._name = localId
        self._datasetId = parentContainer.getId()

    def toProtocolElement(self):
        species = None
        sex = None
        if self.getSpecies():
            species = protocol.fromJson(
                json.dumps(self.getSpecies()), protocol.OntologyTerm)
        if self.getSex():
            sex = protocol.fromJson(
                json.dumps(self.getSex()), protocol.OntologyTerm)
        gaIndividual = protocol.Individual(
            dataset_id=self._datasetId,
            created=self.getCreated(),
            updated=self.getUpdated(),
            description=self.getDescription(),
            id=self.getId(),
            name=self.getName(),
            species=species,
            sex=sex)
        self.serializeAttributes(gaIndividual)
        return gaIndividual

    def populateFromRow(self, individualRecord):
        # TODO coerce to types
        self._name = individualRecord.name
        self._created = individualRecord.created
        self._updated = individualRecord.updated
        self._description = individualRecord.description
        self._species = json.loads(individualRecord.species)
        self._sex = json.loads(individualRecord.sex)
        self.setAttributesJson(individualRecord.attributes)
        return self

    def populateFromJson(self, jsonString):
        # TODO validate
        try:
            parsed = protocol.fromJson(jsonString, protocol.Individual)
        except:
            raise exceptions.InvalidJsonException(jsonString)
        self._created = parsed.created
        self._updated = parsed.updated
        self._description = parsed.description
        self._species = protocol.toJsonDict(parsed.species)
        self._sex = protocol.toJsonDict(parsed.sex)
        attributes = {}
        for key in parsed.attributes.attr:
            attributes[key] = {
                "values": protocol.toJsonDict(parsed.attributes.attr[key])}
        self.setAttributes(attributes)
        return self

    def getCreated(self):
        return self._created

    def getUpdated(self):
        return self._updated

    def getDescription(self):
        return self._description

    def getSpecies(self):
        return self._species

    def getSex(self):
        return self._sex

    def getName(self):
        return self._name
