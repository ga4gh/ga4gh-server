"""
Glue to enable script access to ga4gh packages
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys


def ga4ghImportGlue():
    """
    Call this method before importing a ga4gh module in the scripts dir.
    Otherwise, you will be using the installed package instead of
    the development package.
    Assumes a certain directory structure.
    """
    path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(path)
