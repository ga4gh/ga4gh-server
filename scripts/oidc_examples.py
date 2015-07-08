"""
End-to-end tests that run against the Google server
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import shlex

import utils


def getKey():
    doc = utils.getAuthValues()
    key = doc['google']['key']
    return key


def runTests():
    separator = "----------------------"
    path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        '..', 'client_dev.py')
    baseUrl = "https://localhost:8000/current"
    key = getKey()
    keyStr = "--key {0}".format(key)
    workarounds = ""
    minimalOutput = "-O"
    commands = [
        "readgroupsets-search "
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
        splits = shlex.split(cmd)
        cleanCmd = ' '.join(splits)
        utils.log(separator)
        utils.log(cleanCmd)
        utils.log(separator)
        utils.runCommandSplits(splits)


if __name__ == '__main__':
    runTests()
