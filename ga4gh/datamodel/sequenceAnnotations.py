"""
Module responsible for translating sequence annotation data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy

class Attributes(dict):
    """Type defining a collection of attributes associated with various protocol
    records."""
    def add(self, name, value):
        """added a value or list or tuple of values for the specific attribute"""
        attrValues = self.get(name)
        if attrValues == None:
            self[name] = attrValues = []
        if isinstance(value, list) or isinstance(value, tuple):
            attrValues.extend(value)
        else:
            attrValues.append(value)

    def toProtocolElement(self):
        """
        Returns the representation of this Attributes as the corresponding
        ProtocolElement.
        """
        gaAttributes = protocol.Attributes()
        # TODO: only handles string values noe
        for name in self.iterkeys():
            gaAttributes[name] = self[name]
        return gaAttributes

class Feature(object):
    "Class representing a Feature annotations"
    def __init__(self, id, parentIds, featureSetId, referenceName, start, end, featureType, attributes):
        self.__id = id
        self.__parentIds = parentIds
        self.__featureSetId = featureSetId
        self.__referenceName = referenceName
        self.__start = start
        self.__end = end
        self.__featureType = featureType
        self.__attributes = attributes # TODO: deepcopy?

    def toProtocolElement(self):
        """
        Returns the representation of this Feature as the corresponding
        ProtocolElement.
        """
        gaFeature = protocol.Feature()
        gaFeature.id = self.__id
        gaFeature.parentIds = copy.copy(self.__parentIds)
        gaFeature.featureSetId = self.__featureSetId
        gaFeature.referenceName = self.__referenceName
        gaFeature.start = self.__start
        gaFeature.end = self.__end
        gaFeature.featureType = self.__featureType
        gaFeature.attributes = self.__attributes.toProtocolElement()
        return gaFeature

class FeatureSet(object):
    "A set of sequence features annotations"
    def __init__(self, id, datasetId, referenceSetId, name, sourceUri, attributes):
        self.__id = id
        self.__datasetId = datasetId
        self.__referenceSetId = referenceSetId
        self.__name = name
        self.__sourceUri = sourceUri
        self.__attributes = attributes

    def toProtocolElement(self):
        """
        Returns the representation of this FeatureSet as the corresponding
        ProtocolElement.
        """
        gaFeatureSet = protocol.FeatureSet()
        gaFeatureSet.id = self.__id
        gaFeatureSet.datasetId = self..__datasetId
        gaFeatureSet.referenceSetId = self..__referenceSetId
        gaFeatureSet.name = self..__name
        gaFeatureSet.sourceUri = self..__sourceUri
        gaFeatureSet.attributes = self..__attributes.toProtocolElement()
        return gaFeatureSet
