"""
Functionality common to many cli modules
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


def getRawInput(displayString):
    userResponse = raw_input(displayString)
    return userResponse


def addSubparser(subparsers, subcommand, description):
    parser = subparsers.add_parser(
        subcommand, description=description, help=description)
    return parser


def addDisableUrllibWarningsArgument(parser):
    parser.add_argument(
        "--disable-urllib-warnings", default=False, action="store_true",
        help="Disable urllib3 warnings")


def addOutputFileArgument(parser):
    parser.add_argument(
        "--outputFile", "-o", default=None,
        help="the file to write the output to")
