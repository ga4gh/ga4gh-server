"""
Module responsible for translating reference sequence data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob

import ga4gh.protocol as protocol
import ga4gh.sidegraph as sidegraph


def makeLimits(start=0, end=None):
    limits = None
    if type(end) is int:
        if start >= 0 and end > start:
            limits = (start, end)
    return limits


class GraphDatabase(object):

    """
    Implements graph reference sets: Reference sequences together with
    joins, as well as some helper methods to deal with complex queries.
    Makes assumption that dataDir used for initialization contains one
    graph topology Sqlite database (.db file) together with any named
    FASTA files containing the actual sequence string data.
    """

    def __init__(self, dataDir):
        self._dataDir = dataDir
        # TODO: Throw exception if directory contains != 1 database file
        self._dbFile = glob.glob(os.path.join(self._dataDir, "*.db"))[0]

    def searchReferences(self, referenceSetId=None, start=0, end=None):
        """
        """
        limits = makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchReferencesCount()
            referenceDicts = sg.searchReferences(limits)
        references = []
        for rdict in referenceDicts:
            reference = protocol.Reference()
            reference.id = rdict['ID']
            references.append(reference)
        return (count, references)

    def searchReferenceSets(self, md5checksums=None, accessions=None,
                            assemblyId=None, start=0, end=None):
        """
        Returns a list of dictionaries holding reference set info.
        """
        # TODO: For now, just returns all reference sets, ignores arguments.
        limits = makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchReferenceSetsCount()
            referenceSetDicts = sg.searchReferenceSets(limits)
        referenceSets = []
        for rsdict in referenceSetDicts:
            referenceSet = protocol.ReferenceSet()
            referenceSet.id = rsdict['ID']
            referenceSet.includedReferenceSets = []
            referenceSets.append(referenceSet)
        return (count, referenceSets)

    def searchVariantSets(self, datasetIds=None, start=0, end=None):
        """
        Returns a list of dictionaries holding variant set info.
        """
        # TODO: For now, just returns all variant sets, ignores arguments.
        limits = makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchVariantSetsCount()
            variantSetDicts = sg.searchVariantSets(limits)
        variantSets = []
        for vsdict in variantSetDicts:
            variantSet = protocol.VariantSet()
            variantSet.id = vsdict['ID']
            variantSet.datasetId = ""  # TODO
            variantSet.referenceSetId = ""
            variantSet.metadata = []
            variantSets.append(variantSet)
        return (count, variantSets)

    def searchVariants(self):
        """
        """
        # TODO: Refactor into non-iterator SQL backed convention.
        pass

    def searchCallSets(self):
        """
        """
        # TODO: Refactor into non-iterator SQL backed convention.
        pass

    def searchSequences(self, referenceSetId=None, variantSetId=None,
                        start=0, end=None):
        """
        Returns a pair (count, sequences) where count
        is the total number of found objects for the search (without limits)
        and sequences is the requested array of protocol.Sequence objects.
        start and end are optional integer arguments that must satisfy
        start < end.

        This method, like its SQL backed siblings,
        does not follow iterator convention.
        """
        # TODO search by referenceSetId and variantSetId not yet implemented.
        limits = makeLimits(start, end)
        count = 0
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchSequencesCount()
            rawSequences = sg.searchSequences(limits)
        sequences = []
        for rawSeq in rawSequences:
            seq = protocol.Sequence()
            seq.id = rawSeq['ID']
            seq.length = rawSeq['length']
            sequences.append(seq)
        return (count, sequences)

    def searchJoins(self, referenceSetId=None, variantSetId=None,
                    sequenceId=None, start=0, end=None):
        """
        Returns a pair (count, joins) where count
        is the total number of found objects for the search (without limits)
        and joins is the requested array of protocol.Join objects.
        start and end are optional integer arguments that must satisfy
        start < end.

        This method, like its SQL backed siblings,
        does not follow iterator convention.
        """
        # TODO search by referenceSetId and variantSetId not yet implemented.
        # not to mention by sequenceId!
        limits = makeLimits(start, end)
        count = 0
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchJoinsCount()
            rawJoins = sg.searchJoins(limits)
        joins = []
        for rj in rawJoins:
            join = protocol.Join()
            # deep pile of objecty dodo
            join.side1 = protocol.Side()
            join.side1.base = protocol.Position()
            join.side1.strand = protocol.Strand.POS_STRAND\
                if rj['side1StrandIsForward'] == 'TRUE'\
                else protocol.Strand.NEG_STRAND
            join.side1.base.position = rj['side1Position']
            join.side1.base.sequenceId = rj['side1SequenceID']

            join.side2 = protocol.Side()
            join.side2.base = protocol.Position()
            join.side2.strand = protocol.Strand.POS_STRAND\
                if rj['side2StrandIsForward'] == 'TRUE'\
                else protocol.Strand.NEG_STRAND
            join.side2.base.position = rj['side2Position']
            join.side2.base.sequenceId = rj['side2SequenceID']

            joins.append(join)
        return (count, joins)

    def getSequenceBases(self, sequenceId, start=0, end=None):
        ret = protocol.GetSequenceBasesResponse()
        ret.offset = start
        ret.nextPageToken = None
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            ret.sequence = sg.getSequenceBases(sequenceId, start, end)
        return ret

    def makeSegment(self, segmentDict):
        """
        returns a properly formatted protocol.Segment object.
        The input, segmentDict, a dictionary with keys:
        {sequenceID, start, length, strandIsForward}
        """
        ret = protocol.Segment()
        ret.length = segmentDict["length"]
        ret.start = protocol.Side()
        ret.start.strand = protocol.Strand.POS_STRAND\
            if segmentDict["strandIsForward"] == 'TRUE'\
            else protocol.Strand.NEG_STRAND
        ret.start.base = protocol.Position()
        ret.start.base.sequenceId = segmentDict["sequenceID"]
        ret.start.base.position = segmentDict["start"]
        return ret

    def getAllele(self, alleleId):
        """
        Returns a GA4GH formatted allele with the requested ID.
        TODO: Crashes miserably if ID is not valid.
        """
        ret = protocol.Allele()
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            ret.id = alleleId
            ret.variantSetId = sg.getVariantSetIdForAllele(alleleId)
            ret.path = protocol.Path()
            ret.path.segments = map(self.makeSegment,
                                    sg.getAllelePathItems(alleleId))
        return ret
