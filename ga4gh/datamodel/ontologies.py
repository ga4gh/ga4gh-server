"""
Support for Ontologies.
FIXME: this is currently hard-coded for the GENCODE set as a short-term hack
to support BRCA changeling.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.protocol as protocol
import ga4gh.datamodel as datamodel


class OntologyTermSet(datamodel.DatamodelObject):
    """
    A set of related ontology terms.
    """
    def __init__(self, ontologySource):
        self._ontologySource = ontologySource
        self._ontologyTermIdMap = dict()
        self._ontologyTermNameMap = dict()

    def add(self, ontologyTerm):
        """
        add an ontology term to the object
        """
        if ontologyTerm.id in self._ontologyTermIdMap:
            raise ValueError(
                "OntologyTerm id already exists: {}".format(str(ontologyTerm)))
        self._ontologyTermIdMap[ontologyTerm.id] = ontologyTerm
        if ontologyTerm.name is not None:
            if ontologyTerm.name in self._ontologyTermNameMap:
                raise ValueError(
                    "OntologyTerm already exists: {}".format(str(
                                                             ontologyTerm)))

    def create(self, id, name):
        """
        create a new ontology term
        """
        self.add(OntologyTerm(self._ontologySource, id, name))


class OntologyTerm(datamodel.DatamodelObject):
    """
    A specific ontology term.
    """
    def __init__(self, ontologySource, id, name):
        self._ontologySource = ontologySource
        self._id = id
        self._name = name

    def __str__(self):
        return "{}/{}/{}".format(self.ontologySource, self.id, str(self.name))

    def toProtocolElement(self):
        """
        Returns the representation of this OntologyTerm as the corresponding
        ProtocolElement.
        """
        gaOntologyTerm = protocol.OntologyTerm()
        gaOntologyTerm.ontologySource = self._ontologySource
        gaOntologyTerm.id = self._id
        gaOntologyTerm.name = self._name
        return gaOntologyTerm


class SequenceOntologyTermSet(OntologyTermSet):
    # TODO: tmp hack to encode the sequence ontology terms used by the BRCA
    # gene sets
    def __init__(self):
        super(SequenceOntologyTermSet, self).__init__(
            "http://www.sequenceontology.org/")
        self.create("SO:0000316", "CDS")
        self.create("SO:0000147", "exon")
        self.create("SO:0000704", "gene")
        self.create("SO:0000318", "start_codon")
        self.create("SO:0000319", "stop_codon")
        self.create("SO:0000710", "stop_codon_redefined_as_selenocysteine")
        self.create("SO:0000673", "transcript")
        self.create("SO:0000203", "UTR")

    _singleton = None

    @staticmethod
    def singleton():
        """
        obtain singleton instances of this class
        """
        if SequenceOntologyTermSet._singleton is None:
            SequenceOntologyTermSet._singleton = SequenceOntologyTermSet()
        return SequenceOntologyTermSet._singleton
