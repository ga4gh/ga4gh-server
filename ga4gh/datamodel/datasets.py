"""
Module responsible for translating variant data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.datamodel as datamodel
import ga4gh.datamodel.variants as variants

class Dataset(datamodel.DataObject):
    """
    Class representing a single dataset in the GA4GH data model. A
    dataset corresponds to a single directory in the file system
    that contains a set of reads and variants.
    """
    def __init__(self, dataDir):
        super(Dataset, self).__init__(dataDir)
        self._variantSetIdMap = {}
        # Allocate the variantSets
        variantSetsDir = os.path.join(dataDir, "variants")
        for entry in os.listdir(variantSetsDir):
            relativePath = os.path.join(variantSetsDir, entry)
            variantSet = variants.variantSetFactory(relativePath, self)
            variantSetId = variantSet.getId()
            self._variantSetIdMap[variantSetId] = variantSet
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

