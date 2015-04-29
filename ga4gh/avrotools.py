"""
Provides additional functionality based on avro
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.protocol as protocol


# constants from avro
INT_MIN_VALUE = -(1 << 31)
INT_MAX_VALUE = (1 << 31) - 1
LONG_MIN_VALUE = -(1 << 63)
LONG_MAX_VALUE = (1 << 63) - 1


class ValidationToolException(Exception):
    """
    Something that went wrong in the ValidationTool
    """


class ValidationTool(object):
    """
    Provides additional validation methods based on avro
    """
    sinkValue = 'aRandomStringBecauseWeCanNotUseNone'

    @classmethod
    def getInvalidFields(cls, class_, jsonDict):
        """
        Takes a ProtocolElement subclass class_ and a jsonDict instantiation
        of that class.  If the jsonDict is valid for that class, return
        an empty dictionary.  Otherwise, return a multilevel dictionary
        of the invalid data members and only the invalid data members.

        This is an expensive function to run and shouldn't be on the main
        data path.
        """
        cls._ensureProtocolSubclass(class_)
        cls._ensureNoExtraFields(class_, jsonDict)
        fields = cls._getInvalidFields(class_.schema, jsonDict)
        if fields == cls.sinkValue:
            return {}
        else:
            return fields

    @classmethod
    def _ensureProtocolSubclass(cls, class_):
        if not issubclass(class_, protocol.ProtocolElement):
            message = "class '{}' is not subclass of ProtocolElement".format(
                class_.__name__)
            raise ValidationToolException(message)

    @classmethod
    def _ensureNoExtraFields(cls, class_, jsonDict):
        # TODO we could probably move this inside _getInvalidFields
        # when the underlying object is a python dict
        schema = class_.schema
        schemaJson = schema.to_json()
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
                raise ValidationToolException(message)
        if len(fieldNamesSet) != 0:
            message = "jsonDict does not contain field '{}' for " \
                "class '{}'".format(fieldNamesSet.pop(), class_.__name__)
            raise ValidationToolException(message)

    # Adapted from
    # https://github.com/apache/avro/blob/trunk/lang/py/src/avro/io.py#L101
    @classmethod
    def _getInvalidFields(cls, expectedSchema, datum):
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
                    INT_MIN_VALUE <= datum <= INT_MAX_VALUE):
                return datum
        elif schemaType == 'long':
            if not ((isinstance(datum, int) or isinstance(datum, long)) and
                    LONG_MIN_VALUE <= datum <= LONG_MAX_VALUE):
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
                result = cls._getInvalidFields(expectedSchema.items, dat)
                if result != cls.sinkValue:
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
                result = cls._getInvalidFields(expectedSchema.values, value)
                if result != cls.sinkValue:
                    dic[key] = result
            if len(dic):
                return dic
        elif schemaType in ['union', 'error_union']:
            # If a union type is invalid, we can't tell what was the
            # intended type was to be evaluated.  Therefore, we just
            # take an arbitrary non-null type.
            returnVal = datum
            for schema in expectedSchema.schemas:
                result = cls._getInvalidFields(schema, datum)
                if result == cls.sinkValue:
                    return cls.sinkValue
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
                result = cls._getInvalidFields(field.type, value)
                if result != cls.sinkValue:
                    dic[key] = result
            if len(dic):
                return dic
        return cls.sinkValue
