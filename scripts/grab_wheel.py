"""
Helps with installing ga4gh* packages coming from places other than
pypi into one's environment
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import contextlib
import glob
import os
import pip
import shlex
import subprocess
import urlparse

import ga4gh_common.utils as utils
import ga4gh_common.cli as cli


# TODO take this out when we have a similar method in ga4gh-common
@contextlib.contextmanager
def performInDirectory(dirPath):
    """
    Change the current working directory to dirPath before performing
    an operation, then restore the original working directory after
    """
    originalDirectoryPath = os.getcwd()
    try:
        os.chdir(dirPath)
        yield
    finally:
        os.chdir(originalDirectoryPath)


# TODO take this out when we have a similar method in ga4gh-common
def runCommandReturnOutput(cmd):
    """
    Runs a shell command and returns the stdout
    """
    splits = shlex.split(cmd)
    output = subprocess.check_output(splits)
    return output


# TODO take this out when we have a similar method in ga4gh-common
def chomp(line):
    """
    Returns a string stripped of its trailing newline character
    """
    assert line[-1] == '\n'
    return line[:-1]


def uninstallPackage(packageName):
    """
    Uninstalls the specified package.
    If the package is not installed, a warning message is printed,
    which can be ignored.
    """
    pip.main(['uninstall', '-y', packageName])


def getWheelPath(directoryPath):
    """
    Get the file path of an installable wheel under directoryPath
    """
    distDirPath = os.path.join(directoryPath, 'dist')
    wheelFileNames = glob.glob(os.path.join(distDirPath, '*.whl'))
    if len(wheelFileNames) > 1:
        msg = ("More than one wheel found\n"
               "Clear out {} before running again".format(distDirPath))
        utils.log(msg)
        raise Exception(msg)
    if len(wheelFileNames) < 1:
        msg = "No wheel found"
        utils.log(msg)
        raise Exception(msg)
    return wheelFileNames[0]


def installFromRemote(repoUrl):
    """
    Install a package given a git repository's url
    """
    # For now, a workaround is to download the git repositories,
    # set them at the commits that you want to built the packages
    # at, and then use the local functionality of this script
    raise NotImplementedError()


def installFromLocal(repoPath):
    """
    Installs a package from a local git repository
    """
    packageName = None
    with performInDirectory(repoPath):
        # get the package name
        packageName = runCommandReturnOutput("python setup.py --name")
        packageName = chomp(packageName)
        # create the wheel
        cmd = 'python setup.py bdist_wheel --universal'
        utils.runCommand(cmd)
    uninstallPackage(packageName)
    wheelPath = getWheelPath(repoPath)
    # here we do not use the pip library directly, since we want
    # this script to fail if the install fails
    cmd = 'pip install {}'.format(wheelPath)
    utils.runCommand(cmd)


def isUrl(url):
    """
    Determines if the argument is a valid url
    """
    return urlparse.urlparse(url).scheme != ""


def main():
    parser = cli.createArgumentParser("Installs wheels of ga4gh packages")
    parser.add_argument(
        "packageLocation",
        help=("The location of the package.  "
              "Can be either a git URL or a filepath of the root directory "
              "of a git repository"))
    args = parser.parse_args()
    if isUrl(args.packageLocation):
        repoUrl = args.packageLocation
        installFromRemote(repoUrl)
    else:
        repoPath = args.packageLocation
        installFromLocal(repoPath)


if __name__ == '__main__':
    main()
