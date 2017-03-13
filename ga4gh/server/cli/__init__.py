"""
Functionality common to cli modules
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.server

import ga4gh.schemas.protocol as protocol


def addVersionArgument(parser):
    # TODO argparse strips newlines from version output
    versionString = (
        "GA4GH Server Version {}\n"
        "(Protocol Version {})".format(
            ga4gh.server.__version__, protocol.version))
    parser.add_argument(
        "--version", version=versionString, action="version")


def addDisableUrllibWarningsArgument(parser):
    parser.add_argument(
        "--disable-urllib-warnings", default=False, action="store_true",
        help="Disable urllib3 warnings")
