"""
Functions and utility classes for testing.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import StringIO
import functools
import humanize
import itertools
import os
import random
import signal
import sys
import time

import ga4gh.registry as registry
import ga4gh.datasource.simulator as simulator

packageName = 'ga4gh'


def captureOutput(func, *args, **kwargs):
    """
    Runs the specified function and arguments, and returns the
    tuple (stdout, stderr) as strings.
    """
    stdout = sys.stdout
    sys.stdout = StringIO.StringIO()
    stderr = sys.stderr
    sys.stderr = StringIO.StringIO()
    try:
        func(*args, **kwargs)
        stdoutOutput = sys.stdout.getvalue()
        stderrOutput = sys.stderr.getvalue()
    finally:
        sys.stdout.close()
        sys.stdout = stdout
        sys.stderr.close()
        sys.stderr = stderr
    return stdoutOutput, stderrOutput


def zipLists(*lists):
    """
    Checks to see if all of the lists are the same length, and throws
    an AssertionError otherwise.  Returns the zipped lists.
    """
    length = len(lists[0])
    for list_ in lists[1:]:
        assert len(list_) == length
    return zip(*lists)


def getLinesFromLogFile(stream):
    stream.flush()
    stream.seek(0)
    lines = stream.readlines()
    return lines


def getProjectRootFilePath():
    # assumes we're in a directory one level below the project root
    return os.path.dirname(os.path.dirname(__file__))


def getGa4ghFilePath():
    return os.path.join(getProjectRootFilePath(), packageName)


def powerset(iterable, maxSets=None):
    """
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)

    See https://docs.python.org/2/library/itertools.html#recipes
    """
    s = list(iterable)
    return itertools.islice(itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(len(s) + 1)),
        0, maxSets)


# ---------------- Decorators ----------------


class TimeoutException(Exception):
    """
    A process has taken too long to execute
    """


class Timed(object):
    """
    Decorator that times a method, reporting runtime at finish
    """
    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            self.start = time.time()
            result = func(*args, **kwargs)
            self.end = time.time()
            self._report()
            return result
        return wrapper

    def _report(self):
        delta = self.end - self.start
        timeString = humanize.time.naturaldelta(delta)
        print("Finished in {} ({} seconds)".format(timeString, delta))


class Repeat(object):
    """
    A decorator to use for repeating a tagged function.
    The tagged function should return true if it wants to run again,
    and false if it wants to stop repeating.
    """
    defaultSleepSeconds = 0.1

    def __init__(self, sleepSeconds=defaultSleepSeconds):
        self.sleepSeconds = sleepSeconds

    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            while func(*args, **kwargs):
                time.sleep(self.sleepSeconds)
        return wrapper


class Timeout(object):
    """
    A decorator to use for only allowing a function to run
    for a limited amount of time
    """
    defaultTimeoutSeconds = 60

    def __init__(self, timeoutSeconds=defaultTimeoutSeconds):
        self.timeoutSeconds = timeoutSeconds

    def __call__(self, func):

        def _handle_timeout(signum, frame):
            raise TimeoutException()

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                # set the alarm and execute func
                signal.signal(signal.SIGALRM, _handle_timeout)
                signal.alarm(self.timeoutSeconds)
                result = func(*args, **kwargs)
            finally:
                # clear the alarm
                signal.alarm(0)
            return result
        return wrapper


def create_simulated_registry_db(
        db_url="sqlite:///:memory:", random_seed=1, num_datasets=3,
        num_reference_sets=3, num_references_per_reference_set=3,
        num_variant_sets=3, num_calls=3, num_variant_annotation_sets=3,
        num_read_group_sets=3, num_read_groups_per_read_group_set=3,
        num_feature_sets=3):
    """
    Creates an in-memory registry DB, and populates it with random data
    according to the specified parameters.
    """
    registry_db = registry.RegistryDb(db_url)
    registry_db.open()
    registry_db.initialise()
    rng = random.Random()
    rng.seed(random_seed)
    reference_sets = []
    for j in range(num_reference_sets):
        reference_set = simulator.SimulatedReferenceSet(
            "sim_ref_set_{}".format(j), random_seed=random.randint(0, 2**31),
            num_references=num_references_per_reference_set)
        registry_db.add_reference_set(reference_set)
        reference_sets.append(reference_set)
    for j in range(num_datasets):
        dataset = registry.Dataset("sim_ds_{}".format(j))
        registry_db.add_dataset(dataset)
        for k in range(num_calls):
            individual = registry.Individual("sim_ind_{}".format(k))
            individual.dataset = dataset
            registry_db.add(individual)
            bio_sample = registry.BioSample("sim_bio_sample_{}".format(k))
            bio_sample.dataset = dataset
            bio_sample.individual = individual
            registry_db.add(bio_sample)
        bio_samples = list(dataset.bio_samples)
        for k in range(num_variant_sets):
            variant_set = simulator.SimulatedVariantSet(
                "sim_var_set_{}".format(k),
                random_seed=random.randint(0, 2**31),
                bio_samples=bio_samples)
            variant_set.dataset = dataset
            variant_set.reference_set = random.choice(reference_sets)
            for l in range(num_variant_annotation_sets):
                annotation_set = simulator.SimulatedVariantAnnotationSet(
                    "sim_vas_{}".format(l),
                    random_seed=random.randint(0, 2**31))
                variant_set.variant_annotation_sets.append(annotation_set)
            registry_db.add_variant_set(variant_set)
        for k in range(num_read_group_sets):
            read_group_set = simulator.SimulatedReadGroupSet(
                "sim_rg_set_{}".format(k), bio_samples,
                random_seed=random.randint(0, 2**31),
                num_read_groups=num_read_groups_per_read_group_set)
            read_group_set.dataset = dataset
            read_group_set.reference_set = random.choice(reference_sets)
            registry_db.add_read_group_set(read_group_set)
        for k in range(num_feature_sets):
            feature_set = simulator.SimulatedFeatureSet(
                "sim_feature_set_{}".format(k),
                random_seed=random.randint(0, 2**31))
            feature_set.dataset = dataset
            feature_set.reference_set = random.choice(reference_sets)
            registry_db.add_feature_set(feature_set)
    return registry_db
