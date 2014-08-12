"""
Definitions of the GA4GH protocol types.
"""
from __future__ import print_function
from __future__ import division

import json


class ProtocolElement(object):
    
    def to_json(self):
        return json.dumps(self.__dict__)

class GASearchVariantsRequest(ProtocolElement):
    """
    Search for variants.
    """
    def __init__(self, variantSetIds, start, end):
        self.variantSetIds = variantSetIds
        self.start = start
        self.end = end

class GASearchVariantsResponse(ProtocolElement):
    def __init__(self, variants):
        self.variants = variants
        self.nextPageToken = None


