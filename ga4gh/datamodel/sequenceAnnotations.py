"""
Module responsible for translating sequence annotation data
into GA4GH native objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy
import os

import ga4gh.protocol as protocol
import ga4gh.gff3 as gff3

""" TODO: the search as described in the schema
    The annotation set to search within. Either `featureSetId` or
    `parentId` must be non-empty.
    Only return features with on the reference with this name.  One of this
    field or `referenceId` is required.  (case-sensitive, exact match)
    Required. The beginning of the window (0-based, inclusive) for which
    overlapping features should be returned.  Genomic positions are
    non-negative integers less than reference length.  Requests spanning the
    join of circular genomes are represented as two requests one on each side
    of the join (position 0)
    end is also required
    If specified, this query matches only annotations which match one of the
    provided feature types.
    for now do not use the features array in search
"""


class SequenceAnnotation(object):
    def __init__(self, gffFileName, gffFilePath):
        self._gffFile = os.path.join(gffFilePath, gffFileName)
        self._annotationId = os.path.splitext(gffFileName)
        self._gff3Parser = gff3.Gff3Parser(self._gffFile)

    def getFeature(self, start, stop, referenceName=None, referenceId=None,
                   featureSetId=None, parentId=None):
        """
        Parse the gff file, extract Feature data and return as ga4gh object
        """
        gff3Data = self._gff3Parser.parse()
        # if parentId and featureSetId are both missing return nothing or
        # throw exception?  return nothing for now
        featureList = []
        if parentId is not None:
            featureBranch = gff3Data.byFeatureId.get(parentId)
            featureList = featureBranch[0].children
        elif featureSetId is not None:
            # TODO: get by featureSetId
            featureList = []
        # TODO: restrict by references
        for featureData in featureList:
            id = featureData.featureId()
            # TODO: how to generate featureSetIds?
            featureSetId = None
            featureReference = featureData.seqName
            if ((start >= int(featureData.start)) and (stop <= int(
                                                       featureData.stop))):
                gffFeature = Feature(id, copy.copy(featureData.parents),
                                     featureSetId, featureReference,
                                     featureData.start, featureData.end,
                                     featureData.type, featureData.attributes)
                yield gffFeature.toProtocolElement()

    def getFeatureSet(self):
        """
        Parse the gff file, extract FeatureSet data and return as ga4gh object
        """
        # TODO: not implemented
        pass
        # gff3Data = self._gff3Parser.parse()
        # TODO: extract what we need
        # gffFeatureSet = FeatureSet(id, datasetId, referenceSetId, name,
        #                            sourceUri, attributes)
        # yield gffFeatureSet.toProtocolElement()


class Attributes(dict):
    """
    Type defining a collection of attributes
    associated with various protocol records.
    """
    def add(self, name, value):
        """
        added a value or list or tuple of values
        for the specific attribute
        """
        attrValues = self.get(name)
        if attrValues is None:
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
    """
    Class representing a Feature annotations
    """
    def __init__(
            self, id, parentIds, featureSetId,
            referenceName, start, end, featureType, attributes):
        self._id = id
        self._parentIds = parentIds
        self._featureSetId = featureSetId
        self._referenceName = referenceName
        self._start = start
        self._end = end
        self._featureType = featureType
        self._attributes = attributes  # TODO: deepcopy?

    def toProtocolElement(self):
        """
        Returns the representation of this Feature as the corresponding
        ProtocolElement.
        """
        gaFeature = protocol.Feature()
        gaFeature.id = self._id
        gaFeature.parentIds = copy.copy(self._parentIds)
        gaFeature.featureSetId = self._featureSetId
        gaFeature.referenceName = self._referenceName
        gaFeature.start = self._start
        gaFeature.end = self._end
        gaFeature.featureType = self._featureType
        gaFeature.attributes = self._attributes.toProtocolElement()
        return gaFeature


class FeatureSet(object):
    """
    A set of sequence features annotations
    """
    def __init__(
            self, id, datasetId, referenceSetId,
            name, sourceUri, attributes):
        self._id = id
        self._datasetId = datasetId
        self._referenceSetId = referenceSetId
        self._name = name
        self._sourceUri = sourceUri
        self._attributes = attributes

    def toProtocolElement(self):
        """
        Returns the representation of this FeatureSet as the corresponding
        ProtocolElement.
        """
        gaFeatureSet = protocol.FeatureSet()
        gaFeatureSet.id = self._id
        gaFeatureSet.datasetId = self._datasetId
        gaFeatureSet.referenceSetId = self._referenceSetId
        gaFeatureSet.name = self._name
        gaFeatureSet.sourceUri = self._sourceUri
        gaFeatureSet.attributes = self._attributes.toProtocolElement()
        return gaFeatureSet
