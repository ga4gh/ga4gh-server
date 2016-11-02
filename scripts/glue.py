"""
Glue to enable script access to ga4gh packages
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
# the following two lines are the same ones as in dev_glue.py
# they enable python to find the ga4gh.server package
import ga4gh
ga4gh.__path__.insert(0, 'ga4gh')


def ga4ghImportGlue():
    """
    Call this method before importing a ga4gh module in the scripts dir.
    Otherwise, you will be using the installed package instead of
    the development package.
    Assumes a certain directory structure.
    """
    path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(path)
