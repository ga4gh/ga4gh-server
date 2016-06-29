"""
Biodata objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import json

import ga4gh.datamodel as datamodel
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


class BioSample(datamodel.DatamodelObject):
    """
    This class represents an abstract BioSample object.
    It sets default values and getters, as well as the
    toProtocolElement function.
    """
    compoundIdClass = datamodel.BioSampleCompoundId

    def __init__(self, parentContainer, localId):
        super(BioSample, self).__init__(parentContainer, localId)
        self._created = datetime.datetime.now().isoformat()
        self._updated = datetime.datetime.now().isoformat()
        self._description = None
        self._disease = None
        self._info = {}
        self._name = localId
        self._individualId = None
        self._datasetId = parentContainer.getId()

    def toProtocolElement(self):
        disease = None
        if self.getDisease():
            disease = protocol.fromJson(
                json.dumps(self.getDisease()), protocol.OntologyTerm)
        bioSample = protocol.BioSample(
            dataset_id=self._datasetId,
            created=self.getCreated(),
            updated=self.getUpdated(),
            description=self.getDescription(),
            id=self.getId(),
            individual_id=self.getIndividualId(),
            name=self.getName(),
            disease=disease)
        for key in self.getInfo():
            for value in self.getInfo()[key]['values']:
                bioSample.info[key].values.add().string_value = value
        return bioSample

    def populateFromJson(self, jsonString):
        try:
            parsed = protocol.fromJson(jsonString, protocol.BioSample)
        except:
            raise exceptions.InvalidJsonException(jsonString)
        self._created = parsed.created
        self._updated = parsed.updated
        self._description = parsed.description
        self._disease = protocol.toJsonDict(parsed.disease)
        self._individualId = parsed.individual_id
        self._info = {}
        for key in parsed.info:
            self._info[key] = {"values": protocol.toJsonDict(parsed.info[key])}
        return self

    def populateFromRow(self, row):
        # TODO coerce to types
        self._created = row[b'created']
        self._updated = row[b'updated']
        self._description = row[b'description']
        self._disease = json.loads(row[b'disease'])
        self._individualId = row[b'individualId']
        self._info = json.loads(row[b'info'])
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

    def getInfo(self):
        return self._info

    def getName(self):
        return self._name


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
        self._info = {}
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
        for key in self.getInfo():
            for value in self.getInfo()[key]['values']:
                gaIndividual.info[key].values.add().string_value = value
        return gaIndividual

    def populateFromRow(self, row):
        # TODO coerce to types
        self._name = row[b'name']
        self._created = row[b'created']
        self._updated = row[b'updated']
        self._description = row[b'description']
        self._species = json.loads(row[b'species'])
        self._sex = json.loads(row[b'sex'])
        self._info = json.loads(row[b'info'])
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
        self._info = {}
        for key in parsed.info:
            self._info[key] = {"values": protocol.toJsonDict(parsed.info[key])}
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

    def getInfo(self):
        return self._info

    def getName(self):
        return self._name
