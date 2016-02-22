"""
The repo manager command line interface
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

import ga4gh.repo_manager as repo_manager
import ga4gh.exceptions as exceptions
import ga4gh.cli_common as cli_common


class AbstractRepoCommandRunner(object):

    def __init__(self, args):
        self.args = args
        self.repoPath = args.repoPath
        self.repoManager = repo_manager.RepoManager(self.repoPath)
        if 'force' in args:
            self.force = args.force

    def confirmRun(self, func, deleteString):
        if self.force:
            func()
        else:
            displayString = (
                "Are you sure you want to delete data in {}? "
                "[y|N] ".format(deleteString))
            userResponse = cli_common.getRawInput(displayString)
            if userResponse.strip() == 'y':
                func()
            else:
                print("Aborted")


class AbstractRepoAddCommandRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(AbstractRepoAddCommandRunner, self).__init__(args)
        self.filePath = args.filePath


class AbstractRepoAddMoveCommandRunner(AbstractRepoAddCommandRunner):

    def __init__(self, args):
        super(AbstractRepoAddMoveCommandRunner, self).__init__(args)
        self.moveMode = args.moveMode


class AbstractRepoDatasetCommandRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(AbstractRepoDatasetCommandRunner, self).__init__(args)
        self.datasetName = args.datasetName


class AbstractRepoDatasetFilepathCommandRunner(
        AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(AbstractRepoDatasetFilepathCommandRunner, self).__init__(args)
        self.filePath = args.filePath
        self.moveMode = args.moveMode


class CheckRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(CheckRunner, self).__init__(args)

    def run(self):
        self.repoManager.check()


class ListRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(ListRunner, self).__init__(args)

    def run(self):
        self.repoManager.list()


class DestroyRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(DestroyRunner, self).__init__(args)

    def run(self):
        def func(): self.repoManager.destroy()
        self.confirmRun(func, 'the Data Repository')


class InitRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(InitRunner, self).__init__(args)

    def run(self):
        self.repoManager.init()


class AddDatasetRunner(AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(AddDatasetRunner, self).__init__(args)

    def run(self):
        self.repoManager.addDataset(self.datasetName)


class RemoveDatasetRunner(AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(RemoveDatasetRunner, self).__init__(args)

    def run(self):
        def func(): self.repoManager.removeDataset(self.datasetName)
        self.confirmRun(func, 'dataset {}'.format(self.datasetName))


class AddReferenceSetRunner(AbstractRepoAddMoveCommandRunner):

    def __init__(self, args):
        super(AddReferenceSetRunner, self).__init__(args)
        self.metadata = {
            'description': args.description,
        }

    def run(self):
        self.repoManager.addReferenceSet(
            self.filePath, self.moveMode, self.metadata)


class RemoveReferenceSetRunner(AbstractRepoCommandRunner):

    def __init__(self, args):
        super(RemoveReferenceSetRunner, self).__init__(args)
        self.referenceSetName = args.referenceSetName

    def run(self):
        def func():
            self.repoManager.removeReferenceSet(self.referenceSetName)
        self.confirmRun(
            func, 'reference set {}'.format(self.referenceSetName))


class AddReadGroupSetRunner(AbstractRepoDatasetFilepathCommandRunner):

    def __init__(self, args):
        super(AddReadGroupSetRunner, self).__init__(args)

    def run(self):
        self.repoManager.addReadGroupSet(
            self.datasetName, self.filePath, self.moveMode)


class RemoveReadGroupSetRunner(AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(RemoveReadGroupSetRunner, self).__init__(args)
        self.readGroupSetName = args.readGroupSetName

    def run(self):
        def func():
            self.repoManager.removeReadGroupSet(
                self.datasetName, self.readGroupSetName)
        self.confirmRun(
            func, 'read group set {}'.format(self.readGroupSetName))


class AddVariantSetRunner(AbstractRepoDatasetFilepathCommandRunner):

    def __init__(self, args):
        super(AddVariantSetRunner, self).__init__(args)

    def run(self):
        self.repoManager.addVariantSet(
            self.datasetName, self.filePath, self.moveMode)


class RemoveVariantSetRunner(AbstractRepoDatasetCommandRunner):

    def __init__(self, args):
        super(RemoveVariantSetRunner, self).__init__(args)
        self.variantSetName = args.variantSetName

    def run(self):
        def func():
            self.repoManager.removeVariantSet(
                self.datasetName, self.variantSetName)
        self.confirmRun(func, 'variant set {}'.format(self.variantSetName))


def addRepoArgument(subparser):
    subparser.add_argument(
        "repoPath", help="the file path of the data repository")


def addForceArgument(subparser):
    subparser.add_argument(
        "-f", "--force", action='store_true',
        default=False, help="do not prompt for confirmation")


def addDatasetNameArgument(subparser):
    subparser.add_argument(
        "datasetName", help="the name of the dataset to create/modify")


def addReadGroupSetNameArgument(subparser):
    subparser.add_argument(
        "readGroupSetName",
        help="the name of the read group set")


def addVariantSetNameArgument(subparser):
    subparser.add_argument(
        "variantSetName",
        help="the name of the variant set")


def addFilePathArgument(subparser):
    subparser.add_argument(
        "filePath", help="the path of the file to be moved into the repo")


def addMoveModeArgument(subparser):
    subparser.add_argument(
        "--moveMode",
        help="move, copy or link the target file (default link)",
        choices=['move', 'copy', 'link'],
        default='link')


def addReferenceSetMetadataArguments(subparser):
    subparser.add_argument(
        "--description",
        help="description of the reference set",
        default="TODO")


def getRepoParser():
    parser = argparse.ArgumentParser(
        description="GA4GH data repository management tool")
    subparsers = parser.add_subparsers(title='subcommands',)
    parser.add_argument(
        "--loud", default=False, action="store_true",
        help="propagate exceptions from RepoManager")

    initParser = cli_common.addSubparser(
        subparsers, "init", "Initialize a data repository")
    initParser.set_defaults(runner=InitRunner)
    addRepoArgument(initParser)

    checkParser = cli_common.addSubparser(
        subparsers, "check", "Check to see if repo is well-formed")
    checkParser.set_defaults(runner=CheckRunner)
    addRepoArgument(checkParser)

    listParser = cli_common.addSubparser(
        subparsers, "list", "List the contents of the repo")
    listParser.set_defaults(runner=ListRunner)
    addRepoArgument(listParser)

    destroyParser = cli_common.addSubparser(
        subparsers, "destroy", "Destroy the repo")
    destroyParser.set_defaults(runner=DestroyRunner)
    addRepoArgument(destroyParser)
    addForceArgument(destroyParser)

    addDatasetParser = cli_common.addSubparser(
        subparsers, "add-dataset", "Add a dataset to the data repo")
    addDatasetParser.set_defaults(runner=AddDatasetRunner)
    addRepoArgument(addDatasetParser)
    addDatasetNameArgument(addDatasetParser)

    removeDatasetParser = cli_common.addSubparser(
        subparsers, "remove-dataset",
        "Remove a dataset from the data repo")
    removeDatasetParser.set_defaults(runner=RemoveDatasetRunner)
    addRepoArgument(removeDatasetParser)
    addDatasetNameArgument(removeDatasetParser)
    addForceArgument(removeDatasetParser)

    addReferenceSetParser = cli_common.addSubparser(
        subparsers, "add-referenceset",
        "Add a reference set to the data repo")
    addReferenceSetParser.set_defaults(runner=AddReferenceSetRunner)
    addRepoArgument(addReferenceSetParser)
    addFilePathArgument(addReferenceSetParser)
    addMoveModeArgument(addReferenceSetParser)
    addReferenceSetMetadataArguments(addReferenceSetParser)

    removeReferenceSetParser = cli_common.addSubparser(
        subparsers, "remove-referenceset",
        "Remove a reference set from the repo")
    removeReferenceSetParser.set_defaults(runner=RemoveReferenceSetRunner)
    addRepoArgument(removeReferenceSetParser)
    removeReferenceSetParser.add_argument(
        "referenceSetName",
        help="the name of the reference set")
    addForceArgument(removeReferenceSetParser)

    addReadGroupSetParser = cli_common.addSubparser(
        subparsers, "add-readgroupset",
        "Add a read group set to the data repo")
    addReadGroupSetParser.set_defaults(runner=AddReadGroupSetRunner)
    addRepoArgument(addReadGroupSetParser)
    addDatasetNameArgument(addReadGroupSetParser)
    addFilePathArgument(addReadGroupSetParser)
    addMoveModeArgument(addReadGroupSetParser)

    removeReadGroupSetParser = cli_common.addSubparser(
        subparsers, "remove-readgroupset",
        "Remove a read group set from the repo")
    removeReadGroupSetParser.set_defaults(runner=RemoveReadGroupSetRunner)
    addRepoArgument(removeReadGroupSetParser)
    addDatasetNameArgument(removeReadGroupSetParser)
    addReadGroupSetNameArgument(removeReadGroupSetParser)
    addForceArgument(removeReadGroupSetParser)

    addVariantSetParser = cli_common.addSubparser(
        subparsers, "add-variantset", "Add a variant set to the data repo")
    addVariantSetParser.set_defaults(runner=AddVariantSetRunner)
    addRepoArgument(addVariantSetParser)
    addDatasetNameArgument(addVariantSetParser)
    addFilePathArgument(addVariantSetParser)
    addMoveModeArgument(addVariantSetParser)

    removeVariantSetParser = cli_common.addSubparser(
        subparsers, "remove-variantset",
        "Remove a variant set from the repo")
    removeVariantSetParser.set_defaults(runner=RemoveVariantSetRunner)
    addRepoArgument(removeVariantSetParser)
    addDatasetNameArgument(removeVariantSetParser)
    addVariantSetNameArgument(removeVariantSetParser)
    addForceArgument(removeVariantSetParser)

    return parser


def repo_main(args=None):
    parser = getRepoParser()
    parsedArgs = parser.parse_args(args)
    if "runner" not in parsedArgs:
        parser.print_help()
    else:
        runner = parsedArgs.runner(parsedArgs)
        try:
            runner.run()
        except exceptions.RepoManagerException as exception:
            print(exception.message)
            if parsedArgs.loud:
                raise
