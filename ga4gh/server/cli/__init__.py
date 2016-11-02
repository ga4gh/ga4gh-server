"""
Functionality common to cli modules
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.server
import ga4gh.server.protocol as protocol

# the maximum value of a long type in avro = 2**63 - 1
# (64 bit signed integer)
# http://avro.apache.org/docs/1.7.7/spec.html#schema_primitive
# TODO in the meantime, this is the max value pysam can handle
# This should be removed once pysam input sanitisation has been
# implemented.
AVRO_LONG_MAX = 2**31 - 1


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


def addOutputFileArgument(parser):
    parser.add_argument(
        "--outputFile", "-o", default=None,
        help="the file to write the output to")
