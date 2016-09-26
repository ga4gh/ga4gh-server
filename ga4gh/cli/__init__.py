"""
Functionality common to cli modules
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import operator

import ga4gh
import ga4gh.protocol as protocol

# the maximum value of a long type in avro = 2**63 - 1
# (64 bit signed integer)
# http://avro.apache.org/docs/1.7.7/spec.html#schema_primitive
# TODO in the meantime, this is the max value pysam can handle
# This should be removed once pysam input sanitisation has been
# implemented.
AVRO_LONG_MAX = 2**31 - 1


class SortedHelpFormatter(argparse.HelpFormatter):
    """
    An argparse HelpFormatter that sorts the flags and subcommands
    in alphabetical order
    """
    def add_arguments(self, actions):
        """
        Sort the flags alphabetically
        """
        actions = sorted(
            actions, key=operator.attrgetter('option_strings'))
        super(SortedHelpFormatter, self).add_arguments(actions)

    def _iter_indented_subactions(self, action):
        """
        Sort the subcommands alphabetically
        """
        try:
            get_subactions = action._get_subactions
        except AttributeError:
            pass
        else:
            self._indent()
            if isinstance(action, argparse._SubParsersAction):
                for subaction in sorted(
                        get_subactions(), key=lambda x: x.dest):
                    yield subaction
            else:
                for subaction in get_subactions():
                    yield subaction
            self._dedent()


def addSubparser(subparsers, subcommand, description):
    parser = subparsers.add_parser(
        subcommand, description=description, help=description)
    return parser


def createArgumentParser(description):
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=SortedHelpFormatter)
    return parser


def addVersionArgument(parser):
    # TODO argparse strips newlines from version output
    versionString = (
        "GA4GH Server Version {}\n"
        "(Protocol Version {})".format(
            ga4gh.__version__, protocol.version))
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
