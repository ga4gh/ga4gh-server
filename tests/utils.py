"""
Functions and utility classes for testing.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import functools
import humanize
import itertools
import os
import random
import signal
import string
import time

import avro.schema

import ga4gh.protocol as protocol
import ga4gh.client as client


packageName = 'ga4gh'


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


def getProjectRootFilePath():
    # assumes we're in a directory one level below the project root
    return os.path.dirname(os.path.dirname(__file__))


def getGa4ghFilePath():
    return os.path.join(getProjectRootFilePath(), packageName)


def makeHttpClient():
    url = "http://example.com"
    debugLevel = 0
    workarounds = set()
    key = "KEY"
    httpClient = client.HttpClient(url, debugLevel, workarounds, key)
    return httpClient


class TimeoutException(Exception):
    """
    A process has taken too long to execute
    """


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


def applyVersion(route):
    return "/{0}{1}".format(protocol.version, route)


def powerset(iterable, maxSets=None):
    """
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)

    See https://docs.python.org/2/library/itertools.html#recipes
    """
    s = list(iterable)
    return itertools.islice(itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(len(s) + 1)),
        0, maxSets)


def valueFromSetGenerator(possibilities=[]):
    """
    Creates a generator which returns a random selection from the supplied list
    of possiblities.
    """
    while True:
        yield random.choice(possibilities)


def noneGenerator():
    """
    Creates a generator which always returns None.
    """
    while True:
        yield None


def booleanGenerator():
    """
    Creates a generator which returns only True and False.
    """
    gen = valueFromSetGenerator([True, False])
    while True:
        yield gen.next()


def integerGenerator(min=-2 ** 31, max=2 ** 31):
    """
    Creates a generator which returns integers between min and max.
    """
    while True:
        yield random.randint(min, max)


def longGenerator(min=-2 ** 63, max=2 ** 63):
    """
    Creates a generator which returns integers between min and max.
    """
    while True:
        yield random.randint(min, max)


def floatGenerator(min=-2 ** 31, max=2 ** 31):
    """
    Creates a generator which returns a floating point number between min and
    max.
    """
    while True:
        yield random.uniform(min, max)


def doubleGenerator(min=-2 ** 63, max=2 ** 63):
    """
    Creates a generator which returns a floating point number between min and
    max.
    """
    while True:
        yield random.uniform(min, max)


def stringGenerator(length=10, characters=string.printable):
    """
    Creates a generator which returns a string with the supplied character set
    (defaults to printable characters).
    """
    gen = valueFromSetGenerator(characters)
    while True:
        yield "".join(itertools.islice(gen, length))


def unionGenerator(generators):
    """
    Creates a generator that returns a value from a randomly chosen generator
    from the supplied list.
    """
    while True:
        yield random.choice(generators).next()


def arrayGenerator(generator, length=10):
    """
    Creates a generator that returns an an array of values taken from the
    supplied generator.
    """
    while True:
        yield list(itertools.islice(generator, length))


def mapGenerator(valueGen, size=10, keyGen=stringGenerator()):
    """
    Creates a generator that returns an a dictionary with keys generated by
    keyGen and values generated by valueGen.
    """
    while True:
        keys = list(itertools.islice(keyGen, size))
        values = list(itertools.islice(valueGen, size))
        yield dict(zip(keys, values))


def recordGenerator(recordSchema):
    """
    Creates a generator which returns a dictionary where keys and values are
    determined by the recordSchema. The recordSchema is a map of the records
    keys to a generator of valid values for that key.
    """
    while True:
        record = {}
        for key, gen in recordSchema.iteritems():
            record[key] = gen.next()
        yield record


class InstanceGenerator:
    """
    A tool for generating instances of avro schemas.
    """
    typeToGenerator = {
        "null": noneGenerator(),
        "boolean": booleanGenerator(),
        "int": integerGenerator(),
        "long": longGenerator(),
        "string": stringGenerator(),
        "float": floatGenerator(),
        "double": doubleGenerator()
    }

    def generateInstance(self, cls):
        """
        When supplied with a subclass of ProtocolElement (not an instance),
        generates a valid instance.
        """
        generator = self.__generatorForSchema(cls.schema)
        rawInstance = generator.next()
        return cls.fromJsonDict(rawInstance)

    def generateInvalidateTypeValue(self, *typeStrings):
        """
        Returns a value that is invalid for the supplied type.
        """
        if any(t in self.typeToGenerator for t in typeStrings):
            remaining = dict((typeString, gen) for typeString, gen
                             in self.typeToGenerator.iteritems()
                             if typeString not in typeStrings)
            return random.choice(remaining.values())
        else:  # Not a primitive type, can't be a bool, so just return boolean
            return self.typeToGenerator["boolean"]

    def __generatorForSchema(self, schema):
        if isinstance(schema, avro.schema.RecordSchema):
            generators = {}
            for field in schema.fields:
                generators[field.name] = self.__generatorForSchema(field.type)
            return recordGenerator(generators)
        elif isinstance(schema, avro.schema.UnionSchema):
            generators = list(self.__generatorForSchema(s)
                              for s in schema.schemas)
            return unionGenerator(generators)
        elif isinstance(schema, avro.schema.MapSchema):
            return mapGenerator(self.__generatorForSchema(schema.values))
        elif isinstance(schema, avro.schema.ArraySchema):
            return arrayGenerator(self.__generatorForSchema(schema.items))
        elif isinstance(schema, avro.schema.PrimitiveSchema):
            return self.__generatorForType(schema.type)
        elif isinstance(schema, avro.schema.EnumSchema):
            return valueFromSetGenerator(schema.symbols)
        else:
            raise Exception("unknown schema type %s" % schema)

    def __generatorForType(self, typeString):
        try:
            return self.typeToGenerator[typeString]
        except KeyError:
            raise Exception("unknown type %s" % typeString)
