"""
Provides classes that take protocol requests, send that request to
the server, and write a particular genomics file type with the results.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


class AbstractConverter(object):
    """
    Abstract base class for converter classes
    """
    def __init__(self, request, outputStream):
        self._request = request
        self._outputStream = outputStream


class SamConverter(AbstractConverter):
    """
    Converts a request to a SAM file
    """
    def __init__(self, searchReadsRequest, outputStream):
        super(SamConverter, self).__init__(
            searchReadsRequest, outputStream)

    def convert(self):
        raise NotImplementedError()


class VcfConverter(AbstractConverter):
    """
    Converts a request to a VCF file
    """
    def __init__(self, searchVariantsRequest, outputStream):
        super(VcfConverter, self).__init__(
            searchVariantsRequest, outputStream)

    def convert(self):
        raise NotImplementedError()
