"""
Data-driven tests for ontologies.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os.path

# TODO it may be a bit circular to use obo_parser as our method of
# accessing ontology information, since this is the method we use
# in the main code. However, other libraries have very heavy dependencies.
import ga4gh.server.datamodel.obo_parser as obo_parser
import ga4gh.server.datamodel.ontologies as ontologies

import tests.datadriven as datadriven
import tests.paths as paths

from ga4gh.schemas.ga4gh.common_pb2 import OntologyTerm
import ga4gh.schemas.protocol as protocol


def testReferenceSets():
    testDataDir = os.path.join(paths.testDataDir, "ontologies")
    pattern = "*.obo"
    for test in datadriven.makeTests(
            testDataDir, OntologyTest, pattern):
        yield test


class OntologyTest(datadriven.DataDrivenTest):
    """
    Data driven test class for ontologies.
    """
    def __init__(self, localId, oboFile):
        super(OntologyTest, self).__init__(localId, oboFile)
        self._oboReader = obo_parser.OBOReader(obo_file=oboFile)

    def getDataModelInstance(self, localId, dataPath):
        ontology = ontologies.Ontology(localId)
        ontology.populateFromFile(dataPath)
        return ontology

    def testProtocolElementValid(self):
        # We don't have a protocol element here so this isn't meaningful.
        pass

    def testGoodMappings(self):
        ontology = self._gaObject
        for term in self._oboReader:
            self.assertIn(term.id, ontology.getTermIds(term.name))
            gaTerm = ontology.getGaTermByName(term.name)
            self.assertTrue(protocol.validate(protocol.toJson(gaTerm),
                                              OntologyTerm))
            self.assertEqual(gaTerm.term, term.name)
            self.assertIn(gaTerm.term_id, ontology.getTermIds(term.name))

    def testBadMappings(self):
        for badName in ["Not a term", None, 1234]:
            self.assertEqual(0, len(self._gaObject.getTermIds(badName)))
