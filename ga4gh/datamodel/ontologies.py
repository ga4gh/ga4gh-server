"""
Support for Ontologies.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.protocol as protocol


class OntologyTermMap(object):
    """
    A bidectional map between ontology names and IDs (e.g. in Sequence
    Ontology we would have "SO:0001583 <-> missense_variant"). This
    implementation uses a tab separated TXT file: "id\tname"

    This is a temporary solution and must be replaced by more comprehensive
    ontology object.
    """
    def __init__(self, localId):
        self._localId = localId
        self._dataUrl = None
        self._sourceName = None
        self._nameIdMap = dict()
        self._idNameMap = dict()

    def populateFromFile(self, dataUrl):
        """
        Populates this ontology map from the specified dataUrl.
        This reads the ontology term name and ID pairs from the
        specified file.
        """
        self._dataUrl = dataUrl
        # TODO shouldn't sourceName be the name of the ontology, (e.g,
        # sequence_ontology?)
        self._sourceName = self._dataUrl
        self._readFile()

    def populateFromRow(self, row):
        """
        Populates this Ontology using values in the specified DB row.
        """
        self._dataUrl = row[b'dataUrl']
        # TODO shouldn't sourceName be the name of the ontology, (e.g,
        # sequence_ontology?)
        self._sourceName = self._dataUrl
        self._readFile()

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

    def _readFile(self):
        self._sourceName = self._dataUrl
        with open(self._dataUrl) as f:
            for line in f:
                # File format: id \t name
                fields = line.rstrip().split("\t")
                self.add(fields[0], fields[1])
