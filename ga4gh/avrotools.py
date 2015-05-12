"""
Provides additional functionality based on avro
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import random
import string

import ga4gh.protocol as protocol


INT_MIN_VALUE = -(1 << 31)
INT_MAX_VALUE = (1 << 31) - 1
LONG_MIN_VALUE = -(1 << 63)
LONG_MAX_VALUE = (1 << 63) - 1


class AvrotoolsException(Exception):
    """
    Something that went wrong in the avrotools classes
    """


class AvroTool(object):
    """
    Base class for avrotool classes
    """
    def __init__(self, class_):
        self.class_ = class_
        self.assertProtocolSubclass()
        self.schema = class_.schema

    def assertProtocolSubclass(self):
        if not issubclass(self.class_, protocol.ProtocolElement):
            message = "class '{}' is not subclass of ProtocolElement".format(
                self.class_.__name__)
            raise AvrotoolsException(message)


class Validator(AvroTool):
    """
    Provides methods for validating schemas
    """
    def getInvalidFields(self, jsonDict):
        """
        Takes a jsonDict instantiation of class class_.  If the jsonDict is
        valid for that class, return an empty dictionary.  Otherwise,
        return a multilevel dictionary of the invalid data members and only the
        invalid data members.

        This is an expensive function to run and shouldn't be on the main
        data path.
        """
        validator = SchemaValidator(self.class_)
        invalidFields = validator.getInvalidFields(jsonDict)
        return invalidFields


class Creator(AvroTool):
    """
    Provides methods for creating instances of protocol classes
    according to different rules
    """
    def __init__(self, class_):
        self.class_ = class_
        self.assertProtocolSubclass()
        self.schema = class_.schema

    def getTypicalInstance(self):
        """
        Return a typical instance of the class
        """
        return self._getInstance(TypicalInstanceCreator)

    def getRandomInstance(self):
        """
        Return a random instance of the class
        """
        return self._getInstance(RandomInstanceCreator)

    def getInvalidInstance(self):
        """
        Return an invalid instance of the class
        """
        return self._getInstance(InvalidInstanceCreator)

    def getDefaultInstance(self):
        """
        Return an instance of the class with required values set
        """
        return self._getInstance(DefaultInstanceCreator)

    def _getInstance(self, creatorClass):
        creator = creatorClass(self.class_)
        instance = creator.getInstance()
        return instance

    def getInvalidField(self, fieldName):
        """
        Return an invalid value of the attribute with name fieldName
        """
        return self._getField(InvalidInstanceCreator, fieldName)

    def _getField(self, creatorClass, fieldName):
        creator = creatorClass(self.class_)
        value = creator.getFieldValue(fieldName)
        return value


class AvroTypeSwitch(object):
    """
    Provides reusable logic for branching based on avro types
    """
    schemaRequired = set(
        ['record', 'enum', 'map', 'array', 'union',
         'request', 'error', 'error_union', 'fixed'])

    def __init__(self, class_):
        self.class_ = class_
        self.schema = class_.schema

    def getInstance(self):
        """
        Return an instance of the class instantiated according to
        the rules defined by the subclass
        """
        jsonDict = self.handleSchema(self.schema)
        instance = self.class_.fromJsonDict(jsonDict)
        return instance

    def getFieldValue(self, fieldName):
        """
        Return the value of a field of the class instantiated according
        to the rules defined by the subclass
        """
        fieldType = None
        for field in self.schema.fields:
            if field.name == fieldName:
                fieldType = field.type
                break
        value = self.handleSchema(fieldType)
        return value

    def handleSchema(self, schema, extra=None):
        """
        Return a jsonDict representation of an instance of an object
        defined by schema instantiated according to the rules defined by the
        subclass
        """
        switch = {
            'null': self.handleNull,
            'boolean': self.handleBoolean,
            'string': self.handleString,
            'bytes': self.handleBytes,
            'int': self.handleInt,
            'long': self.handleLong,
            'float': self.handleFloat,
            'double': self.handleDouble,
            'fixed': self.handleFixed,
            'enum': self.handleEnum,
            'array': self.handleArray,
            'map': self.handleMap,
            'union': self.handleUnion,
            'error_union': self.handleErrorUnion,
            'record': self.handleRecord,
            'error': self.handleError,
            'request': self.handleRequest,
        }
        handler = switch[schema.type]
        return self.handleSchemaDispatch(schema, handler, extra)

    def handleSchemaDispatch(self, schema, handler, extra):
        """
        Invoke the handler provided by the handler switch in handleSchema.
        Override if there is a need for extra information to be passed
        into the individual handlers to invoke them correctly.
        """
        if schema.type in self.schemaRequired:
            value = handler(schema)
        else:
            value = handler()
        return value

    # the below methods just delegate to other methods that the subclasses
    # implement; there's no need to re-implement them in (most) subclasses

    def handleErrorUnion(self, schema):
        return self.handleUnion(schema)

    def handleError(self, schema):
        return self.handleRecord(schema)

    def handleRequest(self, schema):
        return self.handleRecord(schema)


class SchemaValidator(AvroTypeSwitch):
    """
    Provides methods for schema validation
    """
    sinkValue = 'aRandomStringBecauseWeCanNotUseNone'

    def handleSchemaDispatch(self, schema, handler, extra):
        if schema.type in self.schemaRequired:
            value = handler(schema, extra)
        else:
            value = handler(extra)
        return value

    def getInvalidFields(self, jsonDict):
        fields = self.handleSchema(self.schema, jsonDict)
        if fields == self.sinkValue:
            return {}
        else:
            return fields

    def handleNull(self, datum):
        if datum is not None:
            return datum
        else:
            return self.sinkValue

    def handleBoolean(self, datum):
        if not isinstance(datum, bool):
            return datum
        else:
            return self.sinkValue

    def handleString(self, datum):
        if not isinstance(datum, basestring):
            return datum
        else:
            return self.sinkValue

    def handleBytes(self, datum):
        if not isinstance(datum, str):
            return datum
        else:
            return self.sinkValue

    def handleInt(self, datum):
        if not ((isinstance(datum, int) or isinstance(datum, long)) and
                INT_MIN_VALUE <= datum <= INT_MAX_VALUE):
            return datum
        else:
            return self.sinkValue

    def handleLong(self, datum):
        if not ((isinstance(datum, int) or isinstance(datum, long)) and
                LONG_MIN_VALUE <= datum <= LONG_MAX_VALUE):
            return datum
        else:
            return self.sinkValue

    def handleFloat(self, datum):
        if not (isinstance(datum, int) or isinstance(datum, long) or
                isinstance(datum, float)):
            return datum
        else:
            return self.sinkValue

    def handleDouble(self, datum):
        return self.handleFloat(datum)

    def handleFixed(self, schema, datum):
        if not (isinstance(datum, str) and len(datum) == schema.size):
            return datum
        else:
            return self.sinkValue

    def handleEnum(self, schema, datum):
        if datum not in schema.symbols:
            return datum
        else:
            return self.sinkValue

    def handleArray(self, schema, datum):
        if not isinstance(datum, list):
            return datum
        arr = []
        for dat in datum:
            result = self.handleSchema(schema.items, dat)
            if result != self.sinkValue:
                arr.append(result)
        if len(arr):
            return arr
        else:
            return self.sinkValue

    def handleMap(self, schema, datum):
        if not isinstance(datum, dict):
            return datum
        dic = {}
        for key, value in datum.items():
            if not isinstance(key, basestring):
                dic[key] = value
            result = self.handleSchema(schema.values, value)
            if result != self.sinkValue:
                dic[key] = result
        if len(dic):
            return dic
        else:
            return self.sinkValue

    def handleUnion(self, schema, datum):
        # If a union type is invalid, we can't tell what was the
        # intended type was to be evaluated.  Therefore, we just
        # take an arbitrary non-null type.
        returnVal = datum
        for unionSchema in schema.schemas:
            result = self.handleSchema(unionSchema, datum)
            if result == self.sinkValue:
                return self.sinkValue
            elif unionSchema.type != 'null':
                returnVal = result
        return returnVal

    def handleRecord(self, schema, datum):
        if not isinstance(datum, dict):
            return datum
        dic = {}
        datumKeys = set(datum.keys())
        for field in schema.fields:
            key = field.name
            value = datum.get(key)
            result = self.handleSchema(field.type, value)
            if result != self.sinkValue:
                dic[key] = result
            try:
                datumKeys.remove(key)
            except KeyError:
                # field that the schema defines as in the record is
                # missing from the jsonDict
                raise AvrotoolsException(key)
        if len(datumKeys):
            # field that the schema does not define is present
            # in the jsonDict
            raise AvrotoolsException(datumKeys)
        if len(dic):
            return dic
        else:
            return self.sinkValue

    def handleErrorUnion(self, schema, datum):
        return self.handleUnion(schema, datum)

    def handleError(self, schema, datum):
        return self.handleRecord(schema, datum)

    def handleRequest(self, schema, datum):
        return self.handleRecord(schema, datum)


class RandomInstanceCreator(AvroTypeSwitch):
    """
    Generates random instances and values
    """
    def handleNull(self):
        return None

    def handleBoolean(self):
        return random.choice([True, False])

    def handleString(self, length=10, characters=string.printable):
        return ''.join(
            [random.choice(characters) for _ in range(length)])

    def handleBytes(self):
        return self.handleString()

    def handleInt(self, min_=INT_MIN_VALUE, max_=INT_MAX_VALUE):
        return random.randint(min_, max_)

    def handleLong(self, min_=LONG_MIN_VALUE, max_=LONG_MAX_VALUE):
        return random.randint(min_, max_)

    def handleFloat(self, min_=INT_MIN_VALUE, max_=INT_MAX_VALUE):
        return random.uniform(min_, max_)

    def handleDouble(self, min_=LONG_MIN_VALUE, max_=LONG_MAX_VALUE):
        return random.uniform(min_, max_)

    def handleFixed(self, schema):
        return self.handleString(schema.size)

    def handleEnum(self, schema):
        return random.choice(schema.symbols)

    def handleArray(self, schema, length=10):
        return [self.handleSchema(schema.items) for _ in range(length)]

    def handleMap(self, schema, size=10):
        return dict(
            (self.handleString(), self.handleSchema(schema.values))
            for _ in range(size))

    def handleUnion(self, schema):
        chosenSchema = random.choice(schema.schemas)
        return self.handleSchema(chosenSchema)

    def handleRecord(self, schema):
        return dict(
            (field.name, self.handleSchema(field.type))
            for field in schema.fields)


class TypicalInstanceCreator(AvroTypeSwitch):
    """
    Generates typical instances and values
    """
    def handleNull(self):
        raise AvrotoolsException()

    def handleBoolean(self):
        return True

    def handleString(self):
        return 'aString'

    def handleBytes(self):
        return b'someBytes'

    def handleInt(self):
        return 5

    def handleLong(self):
        return 6

    def handleFloat(self):
        return 7.0

    def handleDouble(self):
        return 8.0

    def handleFixed(self, schema):
        return 'x' * schema.size

    def handleEnum(self, schema):
        return schema.symbols[0]

    def handleArray(self, schema):
        return [self.handleSchema(schema.items) for _ in range(2)]

    def handleMap(self, schema):
        return {'key': self.handleSchema(schema.values)}

    def handleUnion(self, schema):
        # return an instance of the first non-null schema in the union
        for memberSchema in schema.schemas:
            if memberSchema.type != 'null':
                return self.handleSchema(memberSchema)
        raise AvrotoolsException()  # should never happen

    def handleRecord(self, schema):
        record = {}
        for field in schema.fields:
            record[field.name] = self.handleSchema(field.type)
        return record


class InvalidInstanceCreator(AvroTypeSwitch):
    """
    Generates invalid instances and values
    """
    def handleNull(self):
        return 1

    def handleBoolean(self):
        return "invalidBoolean"

    def handleString(self):
        return 2

    def handleBytes(self):
        return 3

    def handleInt(self):
        return "invalidInt"

    def handleLong(self):
        return "invalidLong"

    def handleFloat(self):
        return "invalidFloat"

    def handleDouble(self):
        return "invalidDouble"

    def handleFixed(self, schema):
        return 4

    def handleEnum(self, schema):
        return "invalidEnum"

    def handleArray(self, schema):
        return "invalidArray"

    def handleMap(self, schema):
        return "invalidMap"

    def handleUnion(self, schema):
        validTypes = set()
        for subSchema in schema.schemas:
            validTypes.add(subSchema.type)
        if 'string' not in validTypes:
            return 'invalidString'
        elif 'map' not in validTypes:
            return {}
        elif 'array' not in validTypes:
            return []
        # no union should have this many types not in it...
        raise AvrotoolsException()

    def handleRecord(self, schema):
        return "invalidRecord"


class DefaultInstanceCreator(TypicalInstanceCreator):
    """
    Generates typical instances with only the default fields set
    """
    def getInstance(self):
        defaultInstance = self.class_()
        jsonDict = defaultInstance.toJsonDict()
        for fieldName in defaultInstance.__slots__:
            if fieldName in defaultInstance.requiredFields:
                typicalValue = self.getFieldValue(fieldName)
                jsonDict[fieldName] = typicalValue
        instance = self.class_.fromJsonDict(jsonDict)
        return instance
