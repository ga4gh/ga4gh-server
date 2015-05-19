"""
Support for Ontologies.
FIXME: this is currently hard-coded for the GENCODE set as a short-term hack
to support BRCA changeling.
"""

class OntologyTermSet(object):
    """A set of related ontology terms.
    """
    def __init__(self, ontologySource):
        self.ontologySource = ontologySource
        self.byId = dict()
        self.byName = dict()

    def add(self, ontologyTerm):
        "add an ontology term to the object"
        if ontologyTerm.id in self.byId:
            raise ValueError("OntologyTerm id already exists: " + str(ontologyTerm))
        self.byId[ontologyTerm.id] = ontologyTerm
        if ontologyTerm.name is not None:
            if ontologyTerm.name in self.byName:
            raise ValueError("OntologyTerm idan already exists: " + str(ontologyTerm))

    def create(self, id, name):
        "create a new ontology term"
        self.add(OntologyTerm(self.ontologySource, id, name))

        
class OntologyTerm(object):
    """A specific ontology term."""
    def __init__(self, ontologySource, id, name):
        self.__ontologySource = ontologySource
        self.__id = id
        self.__name = name
        
    def __str__(self):
        return self.ontologySource + "/" + self.id + "/" + str(self.name)
    def toProtocolElement(self):
        """
        Returns the representation of this OntologyTerm as the corresponding
        ProtocolElement.
        """
        gaOntologyTerm = protocol.OntologyTerm()
        gaOntologyTerm.ontologySource = self.__ontologySource
        gaOntologyTerm.id = self.__id
        gaOntologyTerm.name = self.__name
        return gaOntologyTerm

class SequenceOntologyTermSet(OntologyTermSet):
    """TODO: tmp hack to encode the sequence ontology terms used by the BRCA
    gene sets"""
    def __init__(self):
        OntologyTermSet.__init__(self, "http://www.sequenceontology.org/"):
        self.create("SO:0000316", "CDS")
        self.create("SO:0000147", "exon")
        self.create("SO:0000704", "gene")
        self.create("SO:0000318", "start_codon")
        self.create("SO:0000319", "stop_codon")
        self.create("SO:0000710", "stop_codon_redefined_as_selenocysteine")
        self.create("SO:0000673", "transcript")
        self.create("SO:0000203", "UTR")

    __singleton = None
    @staticmethod
    def singleton():
        "obtain singleton instances of this class"
        if SequenceOntologyTermSet.__singleton is None:
            SequenceOntologyTermSet.__singleton  = SequenceOntologyTermSet()
        return SequenceOntologyTermSet.__singleton
