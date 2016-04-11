"""
Support for Ontologies.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.protocol as protocol


class Ontology(object):
    """
    At this moment an "Ontology" is just a map between names and IDs (e.g.
    in Sequence Ontology we would have "SO:0001583 <-> missense_variant")
    This is a tempotrary solution and must be replaced by more comprehensive
    ontology object.

    This implementation uses a tab separated TXT file: "id\tname"
    """
    def __init__(self, localId):
        self._localId = localId
        self._dataUrl = None
        self._nameIdMap = dict()
        self._idNameMap = dict()

    def populateFromFile(self, dataUrl):
        """
        Populates this ontology map from the specified dataUrl.
        """
        self._dataUrl = dataUrl
        self.readOntology(self._dataUrl)

    def populateFromRow(self, row):
        """
        Populates this Ontology using values in the specified DB row.
        """
        self._dataUrl = row[b'dataUrl']
        self.readOntology(self._dataUrl)

    def add(self, id_, name):
        """
        Adds an ontology term (id, name pair)

        :param id_: ontology term ID (ex "SO:0000704")
        :param name: corresponding ontology term name (ex "gene")
        """
        self._nameIdMap[id_] = name
        self._idNameMap[name] = id_

    def getDataUrl(self):
        return self._dataUrl

    def getLocalId(self):
        return self._localId

    def getId(self, name, default=""):
        return self._idNameMap.get(name, default)

    def hasId(self, id_):
        return id_ in self._nameIdMap

    def hasName(self, name):
        return name in self._idNameIdMap

    def getName(self, id_, default=""):
        return self._nameIdMap.get(id_, default)

    def getGaTermByName(self, name):
        """
        Returns a GA4GH OntologyTerm object by name.

        :param name: name of the ontology term, ex. "gene".
        :return: GA4GH OntologyTerm object.
        """
        term = protocol.OntologyTerm()
        term.term = name
        term.id = self.getId(name)
        # TODO set source name smarter
        term.sourceName = self._sourceName
        # TODO how do we get the right version?
        term.sourceVersion = None
        return term

    def getGaTermById(self, id_):
        """
        Returns a GA4GH OntologyTerm object by its ontology ID.

        :param name: name of the ontology term, ex. "SO:0000704"
            is the ID for "gene" in the Sequence ontology.
        :return: GA4GH OntologyTerm object.
        """
        term = protocol.OntologyTerm()
        term.term = self.getName(id_)
        term.id = id_
        term.sourceName = self._sourceName
        # TODO how do we get the right version?
        term.sourceVersion = None
        return term

    def readOntology(self, filename):
        """
        Extracts ontology maps (ID's to names and vice versa)
        from a file.

        :param filename: the name of the file with ID, name pairs.
        """
        self._sourceName = filename
        with open(filename) as f:
            for line in f:
                # File format: id \t name
                fields = line.rstrip().split("\t")
                self.add(fields[0], fields[1])
