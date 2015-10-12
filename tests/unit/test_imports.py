"""
Tests that the project's module import graph conforms to certain policies
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy
import fnmatch
import itertools
import logging
import operator
import os
import pprint
import unittest
import sys

from snakefood.util import iter_pyfiles, setup_logging, is_python
from snakefood.find import find_dependencies
from snakefood.find import ERROR_IMPORT, ERROR_SYMBOL, ERROR_UNUSED
from snakefood.fallback.collections import defaultdict
from snakefood.roots import find_roots, relfile

import tests.utils as utils


class TestImports(unittest.TestCase):
    """
    Tests that the import graph:
    - doesn't contain any cycles
    - doesn't violate layering constraints
    """
    @classmethod
    def setUpClass(cls):
        snakefoodScanner = SnakefoodScanner()
        cls.graph = snakefoodScanner.scan()

    def testNoCycles(self):
        checker = ImportGraphCycleChecker(self.graph)
        checker.checkNoCycles()

    def testLayering(self):
        checker = ImportGraphLayerChecker(self.graph)
        checker.checkLayeringEnforced()


##############
# Exceptions #
##############


class ConfigurationException(Exception):
    """
    The configuration of a policy checker is invalid
    """


class PolicyException(Exception):
    """
    The code violates some enforced policy
    """


class SnakefoodScannerException(Exception):
    """
    Something went wrong in the snakefood wrapper
    """


###############
# ImportGraph #
###############


class ImportGraphNodeColor(object):
    """
    Node color constants for cycle detection
    """
    WHITE = "WHITE"  # unvisited
    BLACK = "BLACK"  # visited
    GREY = "GREY"  # currently visiting


class ImportGraphNode(object):
    """
    A node in the import graph
    """
    def __init__(self, entry):
        self.name = entry.from_filename
        self.deps = set([entry.to_filename])
        # below fields are for cycle detection
        self.color = ImportGraphNodeColor.WHITE

    def __repr__(self):
        return "ImportGraphNode: {} -> {}".format(
            self.name, repr(list(self.deps)))


class ImportGraph(object):
    """
    A directed graph of import relationships.
    Nodes are files/modules and edges are dependencies.
    """
    def __init__(self):
        self.graph = {}

    def iterNodes(self):
        return self.graph.items()

    def getNodeFor(self, name):
        return self.graph[name]

    def addEntry(self, entry):
        if entry.to_filename is None:
            return
        if entry.from_filename in self.graph:
            self.graph[entry.from_filename].deps.add(entry.to_filename)
        else:
            node = ImportGraphNode(entry)
            self.graph[entry.from_filename] = node

    def hasAnyDependencies(self, name):
        return name in self.graph and len(self.graph[name].deps) != 0

    def hasDependencyOn(self, name, dependency):
        if name not in self.graph:
            return False
        return dependency in self.graph[name].deps

    def removeDependency(self, name, dependency):
        self.graph[name].deps.remove(dependency)

    def printGraph(self):
        pprint.pprint(self.graph)


############
# Checkers #
############


class ImportGraphLayerChecker(object):
    """
    Checks the import graph layering policy

    TODO this class could be more efficient with some optimizations,
    but as it stands the time for these tests is dominated by file
    operations and and parsing the ASTs
    """
    excludedPythonFilenames = set([
        '__init__.py',
    ])

    # each file/module is in one and only one moduleGroup
    moduleGroupNames = {
        'cli': ['ga4gh/cli.py'],
        'client': ['ga4gh/client.py'],
        'frontend': ['ga4gh/frontend.py'],
        'backend': ['ga4gh/backend.py'],
        'exceptions': ['ga4gh/exceptions.py'],
        'datamodel': ['ga4gh/datamodel/reads.py',
                      'ga4gh/datamodel/references.py',
                      'ga4gh/datamodel/variants.py',
                      'ga4gh/datamodel/datasets.py',
                      'ga4gh/datamodel/genotype_phenotype.py'],
        'libraries': ['ga4gh/converters.py',
                      'ga4gh/configtest.py'],
        'protocol': ['ga4gh/protocol.py',
                     'ga4gh/_protocol_definitions.py'],
        'config': ['ga4gh/serverconfig.py'],
        'avrotools': ['ga4gh/avrotools.py'],
    }

    # each moduleGroupName has one and only one entry here
    layers = [
        ['cli'],
        ['client'],
        ['frontend'],
        ['backend'],
        ['libraries'],
        ['datamodel'],
        ['exceptions'],
        ['avrotools'],
        ['config'],
        ['protocol'],
    ]

    def __init__(self, graph):
        self._checkConfiguration()
        self.graph = graph
        self.moduleGroupToOrderIndex = {}
        for i, layerRow in enumerate(self.layers):
            for moduleGroup in layerRow:
                self.moduleGroupToOrderIndex[moduleGroup] = i
        self.moduleToModuleGroup = {}
        for moduleGroup, modules in self.moduleGroupNames.items():
            for module in modules:
                self.moduleToModuleGroup[module] = moduleGroup

    def checkLayeringEnforced(self):
        # rules:
        # - no module can import from modules in layers above it
        # - no module can import from modules in moduleGroups in
        #       same layer as it
        for layer in self.layers:
            for moduleGroup in layer:
                modulesInGroup = self.moduleGroupNames[moduleGroup]
                self._sameLayerCheck(layer, moduleGroup, modulesInGroup)
                self._aboveLayerCheck(layer, moduleGroup, modulesInGroup)

    def _allModules(self):
        modules = list(itertools.chain(*self.moduleGroupNames.values()))
        return modules

    def _checkConfiguration(self):
        # each module that exists in the file tree appears in moduleGroupNames
        pythonFiles = []
        for root, dirnames, filenames in os.walk(utils.getGa4ghFilePath()):
            for filename in fnmatch.filter(filenames, '*.py'):
                pythonFilename = os.path.relpath(
                    os.path.join(root, filename))
                if (pythonFilename not in self.excludedPythonFilenames and
                        filename not in self.excludedPythonFilenames):
                    pythonFiles.append(pythonFilename)
        modules = self._allModules()
        moduleSet = set(modules)
        for pythonFile in pythonFiles:
            if pythonFile not in moduleSet:
                message = "file {} is not listed in moduleGroupNames".format(
                    pythonFile)
                raise ConfigurationException(message)

        # each module should only appear once in moduleGroupNames
        modules = self._allModules()
        moduleSet = set(modules)
        if len(modules) != len(moduleSet):
            for module in moduleSet:
                modules.remove(module)
            message = "duplicate module names in moduleGroupNames: {}"
            raise ConfigurationException(message.format(', '.join(modules)))

        # each moduleGroup should only appear once in layers
        # every defined moduleGroup appears in layers
        moduleGroups = self.moduleGroupNames.keys()
        layersModuleGroups = list(itertools.chain(*self.layers))
        if set(moduleGroups) != set(layersModuleGroups):
            message = "moduleGroupNames and layer moduleGroups not equal"
            raise ConfigurationException(message)

    def _layerIndex(self, layerName):
        return self.moduleGroupToOrderIndex[layerName]

    def _moduleGroupNamesAtSameLayerAs(self, moduleGroup):
        layerIndex = self._layerIndex(moduleGroup)
        layerCopy = self.layers[layerIndex][::]
        layerCopy.remove(moduleGroup)
        return layerCopy

    def _modulesInModuleGroup(self, moduleGroup):
        return self.moduleGroupNames[moduleGroup]

    def _modulesAtSameLayerAs(self, moduleGroup):
        moduleGroupNamesAtSameLayer = self._moduleGroupNamesAtSameLayerAs(
            moduleGroup)
        modules = []
        for moduleGroupName in moduleGroupNamesAtSameLayer:
            layerModules = self._modulesInModuleGroup(moduleGroupName)
            modules.extend(layerModules)
        return modules

    def _modulesAtLayerIndex(self, layerIndex):
        modules = []
        for moduleGroup in self.layers[layerIndex]:
            modules.extend(self._modulesInModuleGroup(moduleGroup))
        return modules

    def _modulesInLayersAbove(self, moduleGroup):
        layerIndex = self._layerIndex(moduleGroup)
        layersAbove = self.layers[:layerIndex]
        modules = []
        for i, layer in enumerate(layersAbove):
            layerModules = self._modulesAtLayerIndex(i)
            modules.extend(layerModules)
        return modules

    def _sameLayerCheck(self, layer, moduleGroup, modulesInGroup):
        modulesAtSameLayer = self._modulesAtSameLayerAs(moduleGroup)
        for module in modulesInGroup:
            for sameLayerModule in modulesAtSameLayer:
                if self.graph.hasDependencyOn(module, sameLayerModule):
                    message = "module '{}' in moduleGroup '{}' " \
                        "has dependency on module '{}' in same layer '{}'"
                    exceptionString = message.format(
                        module, moduleGroup, sameLayerModule, layer)
                    raise PolicyException(exceptionString)

    def _aboveLayerCheck(self, layer, moduleGroup, modulesInGroup):
        modulesAboveLayer = self._modulesInLayersAbove(moduleGroup)
        for module in modulesInGroup:
            for aboveLayerModule in modulesAboveLayer:
                if self.graph.hasDependencyOn(module, aboveLayerModule):
                    group = self.moduleToModuleGroup[aboveLayerModule]
                    message = "module '{}' in moduleGroup '{}' " \
                        "has dependency on module '{}' in moduleGroup '{}'"
                    exceptionString = message.format(
                        module, moduleGroup,
                        aboveLayerModule, group)
                    raise PolicyException(exceptionString)


class ImportGraphCycleChecker(object):
    """
    Checks that there are no cycles in the import graph
    (except those that are explicitly allowed)
    """

    # cyclic dependencies that we want to exclude from validation;
    # essentially, an entry here removes an edge from the dependency
    # graph as far as cycle detection is concerned
    cycleExclusions = [
        ['ga4gh/_protocol_definitions.py', 'ga4gh/protocol.py']
    ]

    def __init__(self, graph):
        self.graph = graph
        self.visitStack = []

    def checkNoCycles(self):
        graph = self._getPreprocessedGraph()
        for name, node in graph.iterNodes():
            if node.color == ImportGraphNodeColor.WHITE:
                self._visitNode(graph, node)

    def _getPreprocessedGraph(self):
        graph = copy.deepcopy(self.graph)
        for name, dependency in self.cycleExclusions:
            graph.removeDependency(name, dependency)
        return graph

    def _visitNode(self, graph, node):
        self.visitStack.append(node)
        node.color = ImportGraphNodeColor.GREY
        for dependency in node.deps:
            if not graph.hasAnyDependencies(dependency):
                continue
            dependencyNode = graph.getNodeFor(dependency)
            if dependencyNode.color == ImportGraphNodeColor.GREY:
                self.visitStack.append(dependencyNode)
                pathString = ' --> '.join(
                    [visited.name for visited in self.visitStack])
                exceptionStr = "Circular import reference: {}".format(
                    pathString)
                raise PolicyException(exceptionStr)
            elif dependencyNode.color == ImportGraphNodeColor.WHITE:
                self._visitNode(graph, dependencyNode)
        node.color = ImportGraphNodeColor.BLACK
        self.visitStack.pop()


#############
# Snakefood #
#############


class SnakefoodEntries(object):
    """
    A list of import entries that snakefood generates
    """
    def __init__(self):
        self.entries = []

    def append(self, entry):
        self.entries.append(entry)

    def printEntries(self):
        pprint.pprint(self.entries)

    def iterEntries(self):
        return iter(self.entries)


class SnakefoodEntry(object):
    """
    An import record that snakefood generates
    """
    def __init__(self, from_root, from_filename, to_root, to_filename):
        self.from_root = from_root
        self.from_filename = from_filename
        self.to_root = to_root
        self.to_filename = to_filename

    def __repr__(self):
        return "SnakefoodEntry: {} -> {}".format(
            self.from_filename, self.to_filename)


class SnakefoodScanner(object):
    """
    Scans for imports within modules in the project.
    Mostly taken from here:
    https://bitbucket.org/blais/snakefood/src/
    e0a74fa6260dcd44716d40b4eb404ca024323eac/
    lib/python/snakefood/gendeps.py?at=default
    """
    def __init__(self):
        self.optsIgnoreUnused = None
        self.optsVerbose = 0
        self.optsDoPragmas = True
        self.optsQuiet = 1
        self.optsInternal = 1
        self.optsExternal = None
        self.optsIgnores = ['.svn', 'CVS', 'build', '.hg', '.git']
        self.optsPrintRoots = None
        self.optsFollow = True
        self.args = [utils.packageName]

    def scan(self):
        """
        Returns an ImportGraph
        """
        self.optsVerbose -= self.optsQuiet
        setup_logging(self.optsVerbose)
        info = logging.info
        warning = logging.warning
        debug = logging.debug
        if self.optsInternal and self.optsExternal:
            message = "Using --internal and --external at the same time " \
                "does not make sense."
            raise SnakefoodScannerException(message)
        if self.optsPrintRoots:
            inroots = find_roots(self.args, self.optsIgnores)
            for dn in sorted(inroots):
                print(dn)
            return
        info("")
        info("Input paths:")
        for arg in self.args:
            fn = os.path.realpath(arg)
            info('  {}'.format(fn))
            if not os.path.exists(fn):
                message = "Filename '{}' does not exist.".format(fn)
                raise SnakefoodScannerException(message)
        # Get the list of package roots for our input files and prepend
        # them to the module search path to insure localized imports.
        inroots = find_roots(self.args, self.optsIgnores)
        if (self.optsInternal or self.optsExternal) and not inroots:
            message = "No package roots found from the given files or " \
                "directories. Using --internal with these roots will  " \
                "generate no dependencies."
            raise SnakefoodScannerException(message)
        info("")
        info("Roots of the input files:")
        for root in inroots:
            info('  {}'.format(root))
        info("")
        info("Using the following import path to search for modules:")
        sys.path = inroots + sys.path
        for dn in sys.path:
            info("  {}".format(dn))
        inroots = frozenset(inroots)
        # Find all the dependencies.
        info("")
        info("Processing files:")
        info("")
        allfiles = defaultdict(set)
        allerrors = []
        processed_files = set()
        fiter = iter_pyfiles(self.args, self.optsIgnores, False)
        while 1:
            newfiles = set()
            for fn in fiter:
                if fn in processed_files:
                    continue  # Make sure we process each file only once.
                info("  {}".format(fn))
                processed_files.add(fn)
                if is_python(fn):
                    files, errors = find_dependencies(
                        fn, self.optsVerbose,
                        self.optsDoPragmas, self.optsVerbose)
                    allerrors.extend(errors)
                else:
                    # If the file is not a source file, we don't know how
                    # to get the dependencies of that (without importing,
                    # which we want to avoid).
                    files = []
                # When packages are the source of dependencies, remove the
                # __init__ file.  This is important because the targets
                # also do not include the __init__ (i.e. when "from
                # <package> import <subpackage>" is seen).
                if os.path.basename(fn) == '__init__.py':
                    fn = os.path.dirname(fn)
                # Make sure all the files at least appear in the output,
                # even if it has no dependency.
                from_ = relfile(fn, self.optsIgnores)
                if from_ is None:
                    continue
                infrom = from_[0] in inroots
                if self.optsInternal and not infrom:
                    continue
                if not self.optsExternal:
                    allfiles[from_].add((None, None))
                # Add the dependencies.
                for dfn in files:
                    xfn = dfn
                    if os.path.basename(xfn) == '__init__.py':
                        xfn = os.path.dirname(xfn)
                    to_ = relfile(xfn, self.optsIgnores)
                    into = to_[0] in inroots
                    if (self.optsInternal and not into) or \
                            (self.optsExternal and into):
                        continue
                    allfiles[from_].add(to_)
                    newfiles.add(dfn)
            if not (self.optsFollow and newfiles):
                break
            else:
                fiter = iter(sorted(newfiles))
        # If internal is used twice, we filter down  further the
        # dependencies to the set of files that were processed only,
        # not just to the files that live in the same roots.
        if self.optsInternal >= 2:
            filtfiles = type(allfiles)()
            for from_, tolist in allfiles.iteritems():
                filtfiles[from_] = set(
                    x for x in tolist if x in allfiles or x == (None, None))
            allfiles = filtfiles
        info("")
        info("SUMMARY")
        info("=======")
        # Output a list of the symbols that could not
        # be imported as modules.
        reports = [
            ("Modules that were ignored because not used:",
                ERROR_UNUSED, info),
            ("Modules that could not be imported:",
                ERROR_IMPORT, warning),
            ]
        if self.optsVerbose >= 2:
            reports.append(
                ("Symbols that could not be imported as modules:",
                    ERROR_SYMBOL, debug))
        for msg, errtype, efun in reports:
            names = set(name for (err, name) in allerrors if err is errtype)
            if names:
                efun("")
                efun(msg)
                for name in sorted(names):
                    efun("  {}".format(name))
        # Output the list of roots found.
        info("")
        info("Found roots:")
        foundRoots = set()
        for key, files in allfiles.iteritems():
            foundRoots.add(key[0])
            foundRoots.update(map(operator.itemgetter(0), files))
        if None in foundRoots:
            foundRoots.remove(None)
        for root in sorted(foundRoots):
            info("  {}".format(root))
        # Output the dependencies.
        entries = SnakefoodEntries()
        info("")
        for (from_root, from_), targets in sorted(
                allfiles.iteritems(), key=operator.itemgetter(0)):
            for to_root, to_ in sorted(targets):
                entry = SnakefoodEntry(from_root, from_, to_root, to_)
                entries.append(entry)
        graph = ImportGraph()
        for entry in entries.iterEntries():
            graph.addEntry(entry)
        return graph
