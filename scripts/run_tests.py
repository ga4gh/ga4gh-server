"""
Script to imitate a Travis CI run
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import subprocess

import yaml

import utils


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
            self.log('Running: "{}"'.format(command))
            subprocess.check_call(command, shell=True)
        self.log('SUCCESS')

    def log(self, logStr):
        utils.log("{0} {1}".format(self.logStrPrefix, logStr))


if __name__ == '__main__':
    travisSimulator = TravisSimulator()
    travisSimulator.runTests()
