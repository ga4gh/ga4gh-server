"""
Data containers aggregate data from a hierarchy of sources.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import bintrees

class DataContainer(object):
    """
    The superclass of all data containers in the GA4GH server.
    A data container is some aggregation of data, arranged in a
    hierarchy. We can view the data model as a tree, in which
    each node contains the IDs and links to child containers.
    Leaf containers are direct interfaces to the data in question,
    which may be represented by a file in the file system, a
    database connection, etc.
    """
    def __init__(self):
        # A map of IDs to child containers.
        self._children = bintrees.AVLTree()
        self._id = None

    def getId(self):
        """
        Returns the ID for this data container.
        """
        return self._id

    def addChild(self, dataContainer):
        """
        Adds the specified child container to this data container.
        """
        childId = dataContainer.getId()
        if childId in self._children:
            currentChild = self._children[childId]
            # TODO shutdown currentChild cleanly.
        self._children[childId] = dataContainer

class FileSystemDataContainer(DataContainer):
    """
    Superclass of all data containers backed by a directory in the
    file system. This corresponds to a directory in the file
    system and represents the aggregation of all of the data
    in this file hierarchy.
    """
    def __init__(self, dataDir, parent):
        self._dataDir = dataDir
        self._parent = parent
        # The ID of this object is the concatentation of the parent id
        # and the directory name.
        self._id = os.path.split(dataDir)[1]
        if parent is not None:
            self._id = "{}:{}".format(parent.getId(), self._id)
        print("Allocated fs new dataobject:", self._id)

    def getParent(self):
        """
        Returns the parent of the DataObject.
        """
        return self._parent
