"""
Dataset objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime

import ga4gh.datamodel as datamodel
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


class AbstractBioSample(datamodel.DatamodelObject):
    """
    This class represents an abstract BioSample object.
    It sets default values and getters, as well as the
    toProtocolElement function.
    """
    compoundIdClass = datamodel.BioSampleCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractBioSample, self).__init__(parentContainer, localId)
        self._createDateTime = datetime.datetime.now().isoformat()
        self._updateDateTime = datetime.datetime.now().isoformat()
        self._description = ""
        self._disease = None
        self._info = {}
        self._name = localId

    def toProtocolElement(self):
        bioSample = protocol.BioSample()
        bioSample.createDateTime = self.getCreateDateTime()
        bioSample.description = self.getDescription()
        bioSample.disease = self.getDisease()
        bioSample.id = self.getId()
        bioSample.info = self.getInfo()
        bioSample.name = self.getName()
        bioSample.updateDateTime = self.getUpdateDateTime()
        return bioSample

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


class JsonBioSample(AbstractBioSample, datamodel.MetadataSidecarMixin):
    """
    Allows a BioSample object to be populated using JSON
    metadata information by passing a filepath during
    instantiation.
    b = JsonBioSample(dataset, name, filepath)
    """
    def __init__(self, parentContainer, localId, filePath):
        super(JsonBioSample, self).__init__(parentContainer, localId)
        self._filePath = filePath
        self.loadSidecar(filePath)
        self.validate()

    def validate(self):
        if not protocol.BioSample.validate(
                self.toProtocolElement().toJsonDict()):
            raise exceptions.FileOpenFailedException(self._filePath)

    def _getField(self, fieldName):
        if self.sidecar(fieldName):
            return self.sidecar(fieldName)
        else:
            return self.__getattribute__('_' + fieldName)

    def getCreateDateTime(self):
        return self._getField('createDateTime')

    def getUpdateDateTime(self):
        return self._getField('updateDateTime')

    def getDescription(self):
        return self._getField('description')

    def getDisease(self):
        return self._getField('disease')

    def getInfo(self):
        return self._getField('info')

    def getName(self):
        return self._getField('name')
