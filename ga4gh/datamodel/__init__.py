"""
The GA4GH data model. Defines all the methods required to translate
data in existing formats into GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

class DataObject(object):
    """
    Superclass of all data objects in the GA4GH server. A data object
    correspond to a directory in the file system and represents the
    aggregation of all of the data in this file hierarchy.
    """
    def __init__(self, dataDir, parent=None):
        self._dataDir = dataDir
        self._parent = parent
        # The ID of this object is the concatentation of the parent id
        # and the directory name.
        self._id = os.path.split(dataDir)[1]
        if parent is not None:
            self._id = "{}:{}".format(parent.getId(), self._id)
        print("Allocated new dataobject:", self._id)

    def getId(self):
        """
        Returns the ID for this data object.
        """
        return self._id

    def getParent(self):
        """
        Returns the parent of the DataObject.
        """
        return self._parent
