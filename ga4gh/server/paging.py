"""
Iterators used by the backend for abstracting the logic of resuming
the object stream at the appropriate point
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import ga4gh.server.exceptions as exceptions


def _parsePageToken(pageToken, numValues):
    """
    Parses the specified pageToken and returns a list of the specified
    number of values. Page tokens are assumed to consist of a fixed
    number of integers seperated by colons. If the page token does
    not conform to this specification, raise a InvalidPageToken
    exception.
    """
    tokens = pageToken.split(":")
    if len(tokens) != numValues:
        msg = "Invalid number of values in page token"
        raise exceptions.BadPageTokenException(msg)
    try:
        values = map(int, tokens)
    except ValueError:
        msg = "Malformed integers in page token"
        raise exceptions.BadPageTokenException(msg)
    return values


def _parseIntegerArgument(args, key, defaultValue):
    """
    Attempts to parse the specified key in the specified argument
    dictionary into an integer. If the argument cannot be parsed,
    raises a BadRequestIntegerException. If the key is not present,
    return the specified default value.
    """
    ret = defaultValue
    try:
        if key in args:
            try:
                ret = int(args[key])
            except ValueError:
                raise exceptions.BadRequestIntegerException(key, args[key])
    except TypeError:
        raise Exception((key, args))
    return ret


class IntervalIterator(object):
    """
    Implements generator logic for types which accept a start/end
    range to search for the object. Returns an iterator over
    (object, pageToken) pairs. The pageToken is a string which allows
    us to pick up the iteration at any point, and is None for the last
    value in the iterator.
    """
    def __init__(self, request, parentContainer):
        self._request = request
        self._parentContainer = parentContainer
        self._searchIterator = None
        self._currentObject = None
        self._nextObject = None
        self._searchAnchor = None
        self._distanceFromAnchor = None
        if not request.page_token:
            self._initialiseIteration()
        else:
            # Set the search start point and the number of records to skip from
            # the page token.
            searchAnchor, objectsToSkip = _parsePageToken(
                request.page_token, 2)
            self._pickUpIteration(searchAnchor, objectsToSkip)

    def _extractProtocolObject(self, obj):
        """
        Returns the protocol object from the object passed back by iteration.
        """
        return obj

    def _initialiseIteration(self):
        """
        Starts a new iteration.
        """
        self._searchIterator = self._search(
            self._request.start,
            self._request.end if self._request.end != 0 else None)
        self._currentObject = next(self._searchIterator, None)
        if self._currentObject is not None:
            self._nextObject = next(self._searchIterator, None)
            self._searchAnchor = self._request.start
            self._distanceFromAnchor = 0
            firstObjectStart = self._getStart(self._currentObject)
            if firstObjectStart > self._request.start:
                self._searchAnchor = firstObjectStart

    def _pickUpIteration(self, searchAnchor, objectsToSkip):
        """
        Picks up iteration from a previously provided page token. There are two
        different phases here:
        1) We are iterating over the initial set of intervals in which start
        is < the search start coorindate.
        2) We are iterating over the remaining intervals in which start >= to
        the search start coordinate.
        """
        self._searchAnchor = searchAnchor
        self._distanceFromAnchor = objectsToSkip
        self._searchIterator = self._search(
            searchAnchor,
            self._request.end if self._request.end != 0 else None)
        obj = next(self._searchIterator)
        if searchAnchor == self._request.start:
            # This is the initial set of intervals, we just skip forward
            # objectsToSkip positions
            for _ in range(objectsToSkip):
                obj = next(self._searchIterator)
        else:
            # Now, we are past this initial set of intervals.
            # First, we need to skip forward over the intervals where
            # start < searchAnchor, as we've seen these already.
            while self._getStart(obj) < searchAnchor:
                obj = next(self._searchIterator)
            # Now, we skip over objectsToSkip objects such that
            # start == searchAnchor
            for _ in range(objectsToSkip):
                if self._getStart(obj) != searchAnchor:
                    raise exceptions.BadPageTokenException
                obj = next(self._searchIterator)
        self._currentObject = obj
        self._nextObject = next(self._searchIterator, None)

    def next(self):
        """
        Returns the next (object, nextPageToken) pair.
        """
        if self._currentObject is None:
            raise StopIteration()
        nextPageToken = None
        if self._nextObject is not None:
            start = self._getStart(self._nextObject)
            # If start > the search anchor, move the search anchor. Otherwise,
            # increment the distance from the anchor.
            if start > self._searchAnchor:
                self._searchAnchor = start
                self._distanceFromAnchor = 0
            else:
                self._distanceFromAnchor += 1
            nextPageToken = "{}:{}".format(
                self._searchAnchor, self._distanceFromAnchor)
        ret = self._extractProtocolObject(self._currentObject), nextPageToken
        self._currentObject = self._nextObject
        self._nextObject = next(self._searchIterator, None)
        return ret

    def __iter__(self):
        return self


class ReadsIntervalIterator(IntervalIterator):
    """
    An interval iterator for reads
    """
    def __init__(self, request, parentContainer, reference):
        self._reference = reference
        super(ReadsIntervalIterator, self).__init__(request, parentContainer)

    def _search(self, start, end):
        return self._parentContainer.getReadAlignments(
            self._reference, start, end)

    @classmethod
    def _getStart(cls, readAlignment):
        if readAlignment.alignment.position.position == 0:
            # unmapped read with mapped mate; see SAM standard 2.4.1
            return readAlignment.next_mate_position.position
        else:
            # usual case
            return readAlignment.alignment.position.position

    @classmethod
    def _getEnd(cls, readAlignment):
        return (
            cls._getStart(readAlignment) +
            len(readAlignment.aligned_sequence))


class VariantsIntervalIterator(IntervalIterator):
    """
    An interval iterator for variants
    """
    def _search(self, start, end):
        return self._parentContainer.getVariants(
            self._request.reference_name, start, end,
            self._request.call_set_ids)

    @classmethod
    def _getStart(cls, variant):
        return variant.start

    @classmethod
    def _getEnd(cls, variant):
        return variant.end


class VariantAnnotationsIntervalIterator(IntervalIterator):
    """
    An interval iterator for annotations
    """
    def __init__(self, request, parentContainer):
        super(VariantAnnotationsIntervalIterator, self).__init__(
            request, parentContainer)
        # TODO do input validation somewhere more sensible
        if self._request.effects is None:
            self._effects = []
        else:
            self._effects = self._request.effects

    def _search(self, start, end):
        return self._parentContainer.getVariantAnnotations(
            self._request.reference_name, start, end)

    def _extractProtocolObject(self, pair):
        variant, annotation = pair
        return annotation

    @classmethod
    def _getStart(cls, pair):
        variant, annotation = pair
        return variant.start

    @classmethod
    def _getEnd(cls, pair):
        variant, annotation = pair
        return variant.end

    def next(self):
        while True:
            ret = super(VariantAnnotationsIntervalIterator, self).next()
            vann = ret[0]
            if self.filterVariantAnnotation(vann):
                return self._removeNonMatchingTranscriptEffects(vann), ret[1]
        return None

    def filterVariantAnnotation(self, vann):
        """
        Returns true when an annotation should be included.
        """
        # TODO reintroduce feature ID search
        ret = False
        if len(self._effects) != 0 and not vann.transcript_effects:
            return False
        elif len(self._effects) == 0:
            return True
        for teff in vann.transcript_effects:
            if self.filterEffect(teff):
                ret = True
        return ret

    def filterEffect(self, teff):
        """
        Returns true when any of the transcript effects
        are present in the request.
        """
        ret = False
        for effect in teff.effects:
            ret = self._matchAnyEffects(effect) or ret
        return ret

    def _checkIdEquality(self, requestedEffect, effect):
        """
        Tests whether a requested effect and an effect
        present in an annotation are equal.
        """
        return self._idPresent(requestedEffect) and (
            effect.term_id == requestedEffect.term_id)

    def _idPresent(self, requestedEffect):
        return requestedEffect.term_id != ""

    def _matchAnyEffects(self, effect):
        ret = False
        for requestedEffect in self._effects:
            ret = self._checkIdEquality(requestedEffect, effect) or ret
        return ret

    def _removeNonMatchingTranscriptEffects(self, ann):
        newTxE = []
        oldTxE = ann.transcript_effects
        if len(self._effects) == 0:
            return ann
        for txe in oldTxE:
            add = False
            for effect in txe.effects:
                if self._matchAnyEffects(effect):
                    add = True
            if add:
                newTxE.append(txe)
        ann.ClearField('transcript_effects')
        ann.transcript_effects.extend(newTxE)
        return ann


class SequenceIterator(object):
    """
    Implements generator logic for types which accept a single number
    to indicate which index the iteration is currently on.
    Implements look-ahead logic for backing stores.
    """
    def __init__(self, request):
        self._request = request
        self._initialize()
        # we need to determine if another object follows the one that is
        # last returned to set nextPageToken correctly, so request an
        # additional object from the database in all cases
        if self._maxResults:
            self._maxResults += 1
        self._nextPageTokenIndex = 0
        self._objectIndex = 0
        if self._request.page_token:
            self._nextPageTokenIndex, = _parsePageToken(
                self._request.page_token, 1)
        self._numToReturn = self._request.page_size
        self._objectList = self._search()
        self._objectListLength = len(self._objectList)

    def _initialize(self):
        """
        Set _startIndex and _maxResults, and any other subclass-specific
        attributes derived from the request object
        """
        raise NotImplementedError()

    def _search(self):
        """
        Actually fetch the objects from the backing store
        """
        raise NotImplementedError()

    def _prepare(self, obj):
        """
        Perform any final transformation on the object before returning
        it to the object stream
        """
        raise NotImplementedError()

    def next(self):
        if (self._numToReturn <= 0 or self._objectIndex >=
                self._objectListLength):
            raise StopIteration()
        obj = self._objectList[self._objectIndex]
        self._nextPageTokenIndex += 1
        nextPageToken = str(self._nextPageTokenIndex)
        if self._objectIndex == self._objectListLength - 1:
            nextPageToken = None
        preparedObj = self._prepare(obj)
        self._objectIndex += 1
        self._numToReturn -= 1
        return preparedObj, nextPageToken

    def __iter__(self):
        return self


class ExpressionLevelsIterator(SequenceIterator):
    """
    Iterates through expression levels
    """
    def __init__(self, request, rnaQuant):
        self._rnaQuant = rnaQuant
        super(ExpressionLevelsIterator, self).__init__(request)

    def _initialize(self):
        self._startIndex = self._request.page_token
        self._maxResults = self._request.page_size

    def _search(self):
        iterator = list(self._rnaQuant.getExpressionLevels(
            threshold=self._request.threshold,
            names=self._request.names,
            startIndex=self._startIndex,
            maxResults=self._maxResults))
        return iterator

    def _prepare(self, obj):
        return obj.toProtocolElement()


class FeaturesIterator(SequenceIterator):
    """
    Iterates through features
    """
    def __init__(self, request, featureSet, parentId):
        self._featureSet = featureSet
        self._parentId = parentId
        super(FeaturesIterator, self).__init__(request)

    def _initialize(self):
        if self._request.start == self._request.end == 0:
            self._start = self._end = None
        else:
            self._start = self._request.start
            self._end = self._request.end
        self._startIndex = self._request.page_token
        self._maxResults = self._request.page_size

    def _search(self):
        iterator = list(self._featureSet.getFeatures(
            self._request.reference_name,
            self._start,
            self._end,
            self._startIndex,
            self._maxResults,
            self._request.feature_types,
            self._parentId,
            self._request.name,
            self._request.gene_symbol))
        return iterator

    def _prepare(self, obj):
        return obj


class ContinuousIterator(SequenceIterator):
    """
    Iterates through continuous data
    """
    def __init__(self, request, continuousSet):
        self._continuousSet = continuousSet
        super(ContinuousIterator, self).__init__(request)

    def _initialize(self):
        if self._request.start == self._request.end == 0:
            self._start = self._end = None
        else:
            self._start = self._request.start
            self._end = self._request.end
        self._startIndex = self._request.page_token
        self._maxResults = self._request.page_size

    def _search(self):
        iterator = list(self._continuousSet.getContinuous(
            self._request.reference_name,
            self._start,
            self._end))
        return iterator

    def _prepare(self, obj):
        return obj


class PeerIterator(SequenceIterator):
    """
    Iterates through peers
    """
    def __init__(self, request, dataRepo):
        self._dataRepo = dataRepo
        super(PeerIterator, self).__init__(request)

    def _initialize(self):
        if self._request.page_token != '':
            self._startIndex = int(self._request.page_token)
        else:
            self._startIndex = 0
        self._maxResults = self._request.page_size

    def _search(self):
        iterator = self._dataRepo.getPeers(
            offset=int(self._startIndex),
            limit=self._maxResults)
        return iterator

    def _prepare(self, obj):
        return obj.toProtocolElement()
