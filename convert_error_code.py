"""
Simple script to decode exception error codes. This translates
the error code received by clients into an exception class.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

import ga4gh.exceptions as exceptions


def parseArgs():
    parser = argparse.ArgumentParser(
        description=(
            "Converts an error code received by a clients to the  "
            "corresponding exception class."))
    parser.add_argument(
        "error_code", type=int,
        help="The error_code value in a GAException object.")
    args = parser.parse_args()
    return args


def main():
    args = parseArgs()
    exceptionClass = exceptions.getExceptionClass(args.error_code)
    print(args.error_code, exceptionClass, sep="\t")


if __name__ == '__main__':
    main()
