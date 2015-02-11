"""
End-to-end tests that run against the Google server
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import subprocess

import yaml


def getKey():
    filepath = 'scripts/auth.yml'
    with open(filepath) as stream:
        doc = yaml.load(stream)
        key = doc['google']['key']
    return key


def runTests():
    separator = "----------------------"
    path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        '..', 'client_dev.py')
    baseUrl = "https://www.googleapis.com/genomics/v1beta2"
    key = getKey()
    keyStr = "--key {0}".format(key)
    workarounds = "--workarounds=google"
    minimalOutput = "-O"
    commands = [
        "variants-search --variantSetIds 10473108253681171589 "
        "--referenceName 22 --start 51005491 --end 51005492 --pageSize 1",
        "references-list-bases --id EIaSo62VtfXT4AE --start 15000 "
        "--end 15010",
        "referencesets-get --id EMud_c37lKPXTQ",
        "references-get --id EIaSo62VtfXT4AE",
        "variantsets-search --datasetIds 10473108253681171589",
        "referencesets-search --accessions GCA_000001405.15",
        "references-search --md5checksums 1b22b98cdeb4a9304cb5d48026a85128",
        "readgroupsets-search --datasetIds 10473108253681171589 "
        "--name NA12878 --pageSize 1",
        "callsets-search --variantSetIds 10473108253681171589 "
        "--name HG00261 --pageSize 1",
        "reads-search --start 51005353 --end 51005354 --readGroupIds "
        "ChhDTXZuaHBLVEZoQ2JpT2J1enZtN25Pb0IQAA --referenceName "
        "22 --pageSize 1",
    ]
    for command in commands:
        cmdStr = """
            python {path} {keyStr} {minimalOutput} {workarounds}
            {command} {baseUrl}
            """
        cmdDict = {
            "path": path,
            "command": command,
            "workarounds": workarounds,
            "keyStr": keyStr,
            "baseUrl": baseUrl,
            "minimalOutput": minimalOutput,
        }
        cmd = cmdStr.format(**cmdDict)
        splits = cmd.split()
        cleanCmd = ' '.join(splits)
        print(separator)
        print(cleanCmd)
        print(separator)
        subprocess.check_call(splits)


if __name__ == '__main__':
    runTests()
