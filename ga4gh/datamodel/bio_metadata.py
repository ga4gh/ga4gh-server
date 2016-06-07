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
        self._info = None
        self._name = localId
        self._individualId = None

    def toProtocolElement(self):
        bioSample = protocol.BioSample(
            created=self.getCreated(),
            updated=self.getUpdated(),
            description=self.getDescription(),
            disease=self.getDisease(),
            id=self.getId(),
            individual_id=self.getIndividualId(),
            info=self.getInfo(),
            name=self.getName())
        return bioSample

    def populateFromJson(self, jsonString):
        try:
            parsed = json.loads(jsonString)
        except:
            raise exceptions.InvalidJsonException(jsonString)
        # TODO validate
        if 'created' in parsed:
            self._created = parsed['created']
        if 'updated' in parsed:
            self._updated = parsed['updated']
        if 'description' in parsed:
            self._description = parsed['description']
        if 'disease' in parsed:
            self._disease = protocol.fromJson(
                json.dumps(parsed['disease']),
                protocol.OntologyTerm)
        if 'individualId' in parsed:
            self._individualId = parsed['individualId']
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
        self._info = None
        self._name = localId

    def toProtocolElement(self):
        gaIndividual = protocol.Individual(
            created=self.getCreated(),
            updated=self.getUpdated(),
            description=self.getDescription(),
            species=self.getSpecies(),
            sex=self.getSex(),
            id=self.getId(),
            info=self.getInfo(),
            name=self.getName())
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
            parsed = json.loads(jsonString)
        except:
            raise exceptions.InvalidJsonException(jsonString)
        if 'created' in parsed:
            self._created = parsed['created']
        if 'updated' in parsed:
            self._updated = parsed['updated']
        if 'description' in parsed:
            self._description = parsed['description']
        if 'species' in parsed:
            self._species = protocol.fromJson(
                json.dumps(parsed['species']),
                protocol.OntologyTerm)
        if 'sex' in parsed:
            self._sex = protocol.fromJson(
                json.dumps(parsed['sex']),
                protocol.OntologyTerm)
        if 'info' in parsed:
            self._info = parsed['info']
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
