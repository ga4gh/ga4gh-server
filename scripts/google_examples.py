"""
End-to-end tests that run against the Google server
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import subprocess


def runTests():
    separator = "----------------------"
    path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        '..', 'client_dev.py')
    baseUrl = "https://www.googleapis.com/genomics/v1beta2"
    keyStr = "--key AIzaSyBcdyeUZSzUVE0e_HPLyCJd7bOZ6g3lr3I"
    workarounds = "--workarounds=google"
    minimalOutput = "-O"
    commands = [
        "variants-search --variantSetIds 10473108253681171589 "
        "--referenceName 22 --start 51005353 --end 51015353 --pageSize 1",
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
