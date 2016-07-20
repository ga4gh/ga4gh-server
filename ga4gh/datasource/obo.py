"""
Support for Ontologies.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import os.path

import sqlalchemy

import ga4gh.protocol as protocol
import ga4gh.registry as registry
import ga4gh.exceptions as exceptions
import ga4gh.parser.obo as obo_parser

SEQUENCE_ONTOLOGY_PREFIX = "SO"


# To make sure we don't parse the same ontology file many times we
# use a cache mapping file names to OntologyFile objects. Note:
# this is just a workaround to avoid reparsing the same ontology
# many times when SqlAlchemy is in charge of creating new instances.
# There may be a better way to do this within SqlAlchemy.
_file_cache = {}


class OboOntology(registry.Ontology):
    """
    A bidectional map between ontology names and IDs (e.g. in Sequence
    Ontology we would have "SO:0001583 <-> missense_variant") derived
    from an OBO file.
    """

    __tablename__ = 'OboOntology'

    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('Ontology.id'),
        primary_key=True)
    data_url = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'OboOntology',
    }

    def __init__(self, name, data_url):
        self.name = name
        self.data_url = data_url
        ontology_file = self.__get_ontology_file()
        self.ontology_prefix = ontology_file.getOntologyPrefix()

    def __get_ontology_file(self):
        if self.data_url not in _file_cache:
            _file_cache[self.data_url] = OntologyFile(self.name, self.data_url)
        return _file_cache[self.data_url]

    def get_ga_term_by_name(self, name):
        """
        Returns a GA4GH OntologyTerm object by name.

        :param name: name of the ontology term, ex. "gene".
        :return: GA4GH OntologyTerm object.
        """
        ontology_file = self.__get_ontology_file()
        return ontology_file.getGaTermByName(name)


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


# TODO refactor this code and change to snake_case

class OntologyFile(object):
    """
    A bidectional map between ontology names and IDs (e.g. in Sequence
    Ontology we would have "SO:0001583 <-> missense_variant") derived
    from an OBO file.
    """
    def __init__(self, name, dataUrl):
        self._id = None
        self._dataUrl = dataUrl
        self._name = name
        self._sourceName = name
        self._sourceVersion = None
        self._ontologyPrefix = None
        # There can be duplicate names, so we need to store a list of IDs.
        self._nameIdMap = collections.defaultdict(list)
        self._readFile()

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
        self._sourceVersion = ""
        if reader.data_version is not None:
            self._sourceVersion = reader.data_version

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

    def getSourceVersion(self):
        """
        The version of the ontology derived from the OBO file.
        """
        return self._sourceVersion

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
        term.id = termId
        term.source_name = self._sourceName
        term.source_version = self._sourceVersion
        return term
