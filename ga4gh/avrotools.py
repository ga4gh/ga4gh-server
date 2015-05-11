"""
Provides additional functionality based on avro
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import itertools
import random
import string

import avro

import ga4gh.protocol as protocol


class SchemaToolException(Exception):
    """
    Something that went wrong in the SchemaTool
    """


class SchemaTool(object):
    """
    Provides additional validation methods based on avro

    Some code adapted from
    https://github.com/apache/avro/blob/trunk/lang/py/src/avro/io.py#L101
    """
    sinkValue = 'aRandomStringBecauseWeCanNotUseNone'
    INT_MIN_VALUE = -(1 << 31)
    INT_MAX_VALUE = (1 << 31) - 1
    LONG_MIN_VALUE = -(1 << 63)
    LONG_MAX_VALUE = (1 << 63) - 1

    def __init__(self, class_):
        self.class_ = class_
        self._ensureProtocolSubclass()
        self.schema = class_.schema
        self.typeToGenerator = None

    def getInvalidFields(self, jsonDict):
        """
        Takes a jsonDict instantiation of class class_.  If the jsonDict is
        valid for that class, return an empty dictionary.  Otherwise, return a
        multilevel dictionary of the invalid data members and only the
        invalid data members.

        This is an expensive function to run and shouldn't be on the main
        data path.
        """
        self._ensureNoExtraFields(jsonDict)
        fields = self._getInvalidFields(self.schema, jsonDict)
        if fields == self.sinkValue:
            return {}
        else:
            return fields

    def getTypicalInstance(self):
        """
        Return a typical instance of the class
        """
        jsonDict = self._getTypicalInstance(self.schema)
        instance = self.class_.fromJsonDict(jsonDict)
        return instance

    def getDefaultInstance(self):
        """
        Return an instance of the class with required values set
        """
        # TODO this could be more efficient
        defaultInstance = self.class_()
        typicalInstance = self.getTypicalInstance()
        for key in defaultInstance.__slots__:
            if key in defaultInstance.requiredFields:
                typicalValue = getattr(typicalInstance, key)
                setattr(defaultInstance, key, typicalValue)
        return defaultInstance

    def getRandomInstance(self):
        randomTool = RandomInstanceGenerator(self.class_)
        instance = randomTool.getInstance()
        return instance

    def _getTypicalInstance(self, schema):
        schemaType = schema.type
        if schemaType == 'null':
            raise SchemaToolException()
        elif schemaType == 'boolean':
            return True
        elif schemaType == 'string':
            return 'aString'
        elif schemaType == 'bytes':
            return 'someBytes'
        elif schemaType == 'int':
            return 5
        elif schemaType == 'long':
            return 6
        elif schemaType in ['float', 'double']:
            return 7.0
        elif schemaType == 'fixed':
            return 'x' * schema.size
        elif schemaType == 'enum':
            return schema.symbols[0]
        elif schemaType == 'array':
            val1 = self._getTypicalInstance(schema.items)
            val2 = self._getTypicalInstance(schema.items)
            arr = [val1, val2]
            return arr
        elif schemaType == 'map':
            value = self._getTypicalInstance(schema.values)
            mapping = {'key': value}
            return mapping
        elif schemaType in ['union', 'error_union']:
            for memberSchema in schema.schemas:
                if memberSchema.type == 'null':
                    continue
                return self._getTypicalInstance(memberSchema)
        elif schemaType in ['record', 'error', 'request']:
            record = {}
            for field in schema.fields:
                record[field.name] = self._getTypicalInstance(field.type)
            return record
        raise SchemaToolException()

    def _ensureProtocolSubclass(self):
        if not issubclass(self.class_, protocol.ProtocolElement):
            message = "class '{}' is not subclass of ProtocolElement".format(
                self.class_.__name__)
            raise SchemaToolException(message)

    def _ensureNoExtraFields(self, jsonDict):
        # TODO we could probably move this inside _getInvalidFields
        # when the underlying object is a python dict
        schemaJson = self.schema.to_json()
        schemaJsonFields = schemaJson['fields']
        fieldNamesSet = set()
        for field in schemaJsonFields:
            fieldName = field['name']
            fieldNamesSet.add(fieldName)
        for fieldName in jsonDict.keys():
            try:
                fieldNamesSet.remove(fieldName)
            except KeyError:
                message = "jsonDict has extra field '{}'".format(fieldName)
                raise SchemaToolException(message)
        if len(fieldNamesSet) != 0:
            message = "jsonDict does not contain field '{}' for " \
                "class '{}'".format(fieldNamesSet.pop(), self.class_.__name__)
            raise SchemaToolException(message)

    def _getInvalidFields(self, expectedSchema, datum):
        schemaType = expectedSchema.type
        if schemaType == 'null':
            if datum is not None:
                return datum
        elif schemaType == 'boolean':
            if not isinstance(datum, bool):
                return datum
        elif schemaType == 'string':
            if not isinstance(datum, basestring):
                return datum
        elif schemaType == 'bytes':
            if not isinstance(datum, str):
                return datum
        elif schemaType == 'int':
            if not ((isinstance(datum, int) or isinstance(datum, long)) and
                    self.INT_MIN_VALUE <= datum <= self.INT_MAX_VALUE):
                return datum
        elif schemaType == 'long':
            if not ((isinstance(datum, int) or isinstance(datum, long)) and
                    self.LONG_MIN_VALUE <= datum <= self.LONG_MAX_VALUE):
                return datum
        elif schemaType in ['float', 'double']:
            if not (isinstance(datum, int) or isinstance(datum, long) or
                    isinstance(datum, float)):
                return datum
        elif schemaType == 'fixed':
            if not (isinstance(datum, str) and
                    len(datum) == expectedSchema.size):
                return datum
        elif schemaType == 'enum':
            if datum not in expectedSchema.symbols:
                return datum
        elif schemaType == 'array':
            if not isinstance(datum, list):
                return datum
            arr = []
            for dat in datum:
                result = self._getInvalidFields(expectedSchema.items, dat)
                if result != self.sinkValue:
                    arr.append(result)
            if len(arr):
                return arr
        elif schemaType == 'map':
            if not isinstance(datum, dict):
                return datum
            dic = {}
            for key, value in datum.items():
                if not isinstance(key, basestring):
                    dic[key] = value
                result = self._getInvalidFields(expectedSchema.values, value)
                if result != self.sinkValue:
                    dic[key] = result
            if len(dic):
                return dic
        elif schemaType in ['union', 'error_union']:
            # If a union type is invalid, we can't tell what was the
            # intended type was to be evaluated.  Therefore, we just
            # take an arbitrary non-null type.
            returnVal = datum
            for schema in expectedSchema.schemas:
                result = self._getInvalidFields(schema, datum)
                if result == self.sinkValue:
                    return self.sinkValue
                elif schema.type != 'null':
                    returnVal = result
            return returnVal
        elif schemaType in ['record', 'error', 'request']:
            if not isinstance(datum, dict):
                return datum
            dic = {}
            for field in expectedSchema.fields:
                key = field.name
                value = datum.get(key)
                result = self._getInvalidFields(field.type, value)
                if result != self.sinkValue:
                    dic[key] = result
            if len(dic):
                return dic
        return self.sinkValue


class RandomInstanceGenerator(object):
    """
    Generates random instances and values
    """
    INT_MIN_VALUE = -(1 << 31)
    INT_MAX_VALUE = (1 << 31) - 1
    LONG_MIN_VALUE = -(1 << 63)
    LONG_MAX_VALUE = (1 << 63) - 1

    @classmethod
    def generateInvalidTypeValue(cls, *typeStrings):
        """
        Returns a value that is invalid for the supplied type.
        """
        if any(t in cls._getTypeToGenerator() for t in typeStrings):
            remaining = dict((typeString, gen) for typeString, gen
                             in cls._getTypeToGenerator().iteritems()
                             if typeString not in typeStrings)
            return random.choice(remaining.values())
        else:  # Not a primitive type, can't be a bool, so just return boolean
            return cls._getTypeToGenerator()["boolean"]

    @classmethod
    def _generatorForType(cls, typeString):
        try:
            return cls._getTypeToGenerator()[typeString]
        except KeyError:
            raise SchemaToolException("unknown type %s" % typeString)

    @classmethod
    def _getTypeToGenerator(cls):
        typeToGenerator = {
            "null": cls._noneGenerator(),
            "boolean": cls._booleanGenerator(),
            "int": cls._integerGenerator(),
            "long": cls._longGenerator(),
            "string": cls._stringGenerator(),
            "float": cls._floatGenerator(),
            "double": cls._doubleGenerator()
        }
        return typeToGenerator

    @classmethod
    def _valueFromSetGenerator(cls, possibilities):
        """
        Creates a generator which returns a random selection from the
        supplied list of possiblities.
        """
        while True:
            yield random.choice(possibilities)

    @classmethod
    def _noneGenerator(cls):
        """
        Creates a generator which always returns None.
        """
        while True:
            yield None

    @classmethod
    def _booleanGenerator(cls):
        """
        Creates a generator which returns only True and False.
        """
        gen = cls._valueFromSetGenerator([True, False])
        while True:
            yield gen.next()

    @classmethod
    def _integerGenerator(cls, min_=None, max_=None):
        """
        Creates a generator which returns integers between min and max.
        """
        if min_ is None:
            min_ = cls.INT_MIN_VALUE
        if max_ is None:
            max_ = cls.INT_MAX_VALUE
        while True:
            yield random.randint(min_, max_)

    @classmethod
    def _longGenerator(cls, min_=None, max_=None):
        """
        Creates a generator which returns integers between min and max.
        """
        if min_ is None:
            min_ = cls.LONG_MIN_VALUE
        if max_ is None:
            max_ = cls.LONG_MAX_VALUE
        while True:
            yield random.randint(min_, max_)

    @classmethod
    def _floatGenerator(cls, min_=None, max_=None):
        """
        Creates a generator which returns a floating point number between min
        and max.
        """
        if min_ is None:
            min_ = cls.INT_MIN_VALUE
        if max_ is None:
            max_ = cls.INT_MAX_VALUE
        while True:
            yield random.uniform(min_, max_)

    @classmethod
    def _doubleGenerator(cls, min_=None, max_=None):
        """
        Creates a generator which returns a floating point number between min
        and max.
        """
        if min_ is None:
            min_ = cls.LONG_MIN_VALUE
        if max_ is None:
            max_ = cls.LONG_MAX_VALUE
        while True:
            yield random.uniform(min_, max_)

    @classmethod
    def _stringGenerator(cls, length=10, characters=string.printable):
        """
        Creates a generator which returns a string with the supplied
        character set (defaults to printable characters).
        """
        gen = cls._valueFromSetGenerator(characters)
        while True:
            yield "".join(itertools.islice(gen, length))

    @classmethod
    def _unionGenerator(cls, generators):
        """
        Creates a generator that returns a value from a randomly chosen
        generator from the supplied list.
        """
        while True:
            yield random.choice(generators).next()

    @classmethod
    def _arrayGenerator(cls, generator, length=10):
        """
        Creates a generator that returns an an array of values taken from the
        supplied generator.
        """
        while True:
            yield list(itertools.islice(generator, length))

    @classmethod
    def _mapGenerator(cls, valueGen, size=10):
        """
        Creates a generator that returns an a dictionary with keys generated
        by keyGen and values generated by valueGen.
        """
        keyGen = cls._stringGenerator()
        while True:
            keys = list(itertools.islice(keyGen, size))
            values = list(itertools.islice(valueGen, size))
            yield dict(zip(keys, values))

    @classmethod
    def _recordGenerator(self, recordSchema):
        """
        Creates a generator which returns a dictionary where keys and values
        are determined by the recordSchema. The recordSchema is a map of the
        records keys to a generator of valid values for that key.
        """
        while True:
            record = {}
            for key, gen in recordSchema.iteritems():
                record[key] = gen.next()
            yield record

    def __init__(self, class_):
        self.class_ = class_
        self.schema = class_.schema

    def getInstance(self):
        generator = self._generateRandomInstance(self.schema)
        jsonDict = generator.next()
        instance = self.class_.fromJsonDict(jsonDict)
        return instance

    def _generateRandomInstance(self, schema):
        if isinstance(schema, avro.schema.RecordSchema):
            generators = {}
            for field in schema.fields:
                generators[field.name] = self._generateRandomInstance(
                    field.type)
            return self._recordGenerator(generators)
        elif isinstance(schema, avro.schema.UnionSchema):
            generators = list(self._generateRandomInstance(s)
                              for s in schema.schemas)
            return self._unionGenerator(generators)
        elif isinstance(schema, avro.schema.MapSchema):
            return self._mapGenerator(
                self._generateRandomInstance(schema.values))
        elif isinstance(schema, avro.schema.ArraySchema):
            return self._arrayGenerator(
                self._generateRandomInstance(schema.items))
        elif isinstance(schema, avro.schema.PrimitiveSchema):
            return self._generatorForType(schema.type)
        elif isinstance(schema, avro.schema.EnumSchema):
            return self._valueFromSetGenerator(schema.symbols)
        else:
            raise SchemaToolException("unknown schema type %s" % schema)
