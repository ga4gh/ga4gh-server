"""
Ontology objects
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import fnmatch
import os


class FileSystemOntology(object):
    """
    The base class of storing an ontology
    At this moment an "Ontology" is just a map between names and IDs (e.g.
    in Sequence Ontology we would have "SO:0001583 <-> missense_variant")
    This is a tempotrary solution adn must be replaced by more comprehensive
    ontology object.
    """

    def __init__(self):
        self._nameIdMap = dict()
        self._idNameMap = dict()

    def add(self, id_, name):
        self._nameIdMap[id_] = name
        self._idNameMap[name] = id_

    def getId(self, name):
        return self._idNameMap[name]

    def getName(self, id_):
        return self._nameIdMap[id_]

    def readOntology(self, filename):
        with open(filename) as f:
            for line in f:
                # File format: id \t name
                fields = line.rstrip().split("\t")
                self.add(fields[0], fields[1])


class FileSystemOntologies(object):
    """
    An ontology read from the file system.
    This implementation uses a tab separated TXT file: "id\tname"
    """

    def __init__(self, localId, dataDir, backend):
        self._ontologyNameMap = dict()
        self.readOntologies(dataDir)

    def add(self, ontologyName, ontology):
        self._ontologyNameMap[ontologyName] = ontology

    def get(self, ontologyName):
        return self._ontologyNameMap[ontologyName]

    def keys(self):
        return self._ontologyNameMap.keys()

    def len(self):
        return len(self._ontologyNameMap)

    def readOntologies(self, dataDir):
        self._dataDir = dataDir
        # Find TXT files
        for filename in os.listdir(dataDir):
            if fnmatch.fnmatch(filename, '*.txt'):
                ontologyName, _ = os.path.splitext(filename)
                path = os.path.join(dataDir, filename)
                ontology = FileSystemOntology()
                ontology.readOntology(path)
                self.add(ontologyName, ontology)
