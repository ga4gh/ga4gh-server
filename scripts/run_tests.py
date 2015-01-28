"""
Script to imitate a Travis CI run
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import subprocess

import yaml


class TravisSimulator(object):

    logStrPrefix = '***'

    yamlFileLocation = '.travis.yml'

    def parseTestCommands(self):
        yamlFile = file(self.yamlFileLocation)
        yamlData = yaml.load(yamlFile)
        return yamlData['script']

    def runTests(self):
        testCommands = self.parseTestCommands()
        for command in testCommands:
            returnCode = self.runCommand(command)
            if returnCode != 0:
                self.log('ERROR')
                return
        self.log('SUCCESS')

    def runCommand(self, command):
        self.log('Running: "{0}"'.format(command))
        splits = command.split()
        s = subprocess.Popen(splits)
        s.communicate()
        return s.returncode

    def log(self, logStr):
        print("{0} {1}".format(self.logStrPrefix, logStr))


if __name__ == '__main__':
    travisSimulator = TravisSimulator()
    travisSimulator.runTests()
