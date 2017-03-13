"""
Support for Ontologies.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import os.path

import ga4gh.server.exceptions as exceptions
import ga4gh.server.datamodel.obo_parser as obo_parser

import ga4gh.schemas.protocol as protocol


SEQUENCE_ONTOLOGY_PREFIX = "SO"


class OboReader(obo_parser.OBOReader):
    """
    We extend the OBOReader class to allow us throw a custom exception
    so that it will be handled correctly. The default implementation
    throws an Exception instance on error, which cannot be caught
    without masking pretty much any kind of error.
    """
    def _die(self, message, line):
        raise exceptions.OntologyFileFormatException(
            self.obo_file, "Error at line {}: {}".format(line, message))


class Ontology(object):
    """
    A bidectional map between ontology names and IDs (e.g. in Sequence
    Ontology we would have "SO:0001583 <-> missense_variant") derived
    from an OBO file.
    """
    def __init__(self, name):
        self._id = None
        self._name = name
        self._sourceVersion = None
        self._ontologyPrefix = None
        self._dataUrl = None
        # There can be duplicate names, so we need to store a list of IDs.
        self._nameIdMap = collections.defaultdict(list)

    def _readFile(self):
        if not os.path.exists(self._dataUrl):
            raise exceptions.FileOpenFailedException(self._dataUrl)
        reader = OboReader(obo_file=self._dataUrl)
        ids = set()
        for record in reader:
            if record.id in ids:
                raise exceptions.OntologyFileFormatException(
                    self._dataUrl, "Duplicate ID {}".format(record.id))
            ids.add(record.id)
            self._nameIdMap[record.name].append(record.id)
        self._sourceVersion = reader.format_version
        if len(ids) == 0:
            raise exceptions.OntologyFileFormatException(
                self._dataUrl, "No valid records found.")
        # To get prefix, pull out an ID and parse it.
        self._ontologyPrefix = record.id.split(":")[0]
        self._sourceVersion = reader.data_version

    def populateFromFile(self, dataUrl):
        """
        Populates this ontology map from the specified dataUrl.
        This reads the ontology term name and ID pairs from the
        specified file.
        """
        self._dataUrl = dataUrl
        self._readFile()

    def populateFromRow(self, ontologyRecord):
        """
        Populates this Ontology using values in the specified DB row.
        """
        self._id = ontologyRecord.id
        self._dataUrl = ontologyRecord.dataurl
        self._readFile()
        # TODO sanity check the stored values against what we have just read.

    def getId(self):
        """
        Returns the ID of this Ontology. This is an internal
        identifier.
        """
        return self._id

    def getOntologyPrefix(self):
        """
        Returns the ontology prefix string, i.e. "SO" for a sequence
        ontology and "GO" for gene a ontology.
        """
        return self._ontologyPrefix

    def getDataUrl(self):
        return self._dataUrl

    def getName(self):
        """
        Returns the name of this ontology.
        """
        return self._name

    def getTermIds(self, termName):
        """
        Returns the list of ontology IDs scorresponding to the specified term
        name. If the term name is not found, return the empty list.
        """
        return self._nameIdMap[termName]

    def getGaTermByName(self, name):
        """
        Returns a GA4GH OntologyTerm object by name.

        :param name: name of the ontology term, ex. "gene".
        :return: GA4GH OntologyTerm object.
        """
        # TODO what is the correct value when we have no mapping??
        termIds = self.getTermIds(name)
        if len(termIds) == 0:
            termId = ""
            # TODO add logging for missed term translation.
        else:
            # TODO what is the correct behaviour here when we have multiple
            # IDs matching a given name?
            termId = termIds[0]
        term = protocol.OntologyTerm()
        term.term = name
        term.term_id = termId
        return term
