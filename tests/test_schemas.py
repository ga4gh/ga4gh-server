"""
Tests for auto generated schemas and conversion to and from JSON.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import string
import random
import inspect
import unittest

import avro.schema

import tests
import ga4gh.protocol as protocol

# TODO We should add some tests here - see comment from @pashields on
# https://github.com/ga4gh/server/pull/39


def randomString():
    """
    Returns a randomly generated short string.
    """
    randInt = random.randint(0, 10)
    randStr = ''.join(random.choice(
        string.ascii_letters) for _ in range(randInt))
    return randStr


class SchemaTest(unittest.TestCase):
    """
    Superclass of schema tests.
    """
    typicalValueMap = {
        "string": "string value",
        "int": 1000,
        "long": 10000,
        "boolean": True,
        "double": 0.125,
        "float": 0.25
    }
    badValueMap = {
        "string": 42,
        "int": "string",
        "long": "string",
        "boolean": {},
        "double": {},
        "float": {},
    }
    randomValueMap = {
        "string": randomString,
        "int": lambda: random.randint(-42, 42),
        "long": lambda: random.randint(-42, 42),
        "boolean": lambda: bool(random.randint(0, 1)),
        "double": lambda: random.uniform(-42, 42),
        "float": lambda: random.uniform(-42, 42),
    }

    def getProtocolClasses(self):
        """
        Returns all the protocol classes defined by the schema.
        """
        for name, obj in inspect.getmembers(protocol):
            if inspect.isclass(obj):
                # We're only interested in sublasses of ProtocolElement
                pe = protocol.ProtocolElement
                if issubclass(obj, pe) and obj is not pe:
                    yield obj

    def getAvroSchema(self, cls, fieldName):
        """
        Returns the avro schema for the specified field.
        """
        field = None
        for fld in cls.schema.fields:
            if fld.name == fieldName:
                field = fld
        return field

    def getInvalidValue(self, cls, fieldName):
        """
        Returns a value that should trigger a schema validation failure.
        """
        field = self.getAvroSchema(cls, fieldName)
        typ = field.type
        if isinstance(typ, avro.schema.UnionSchema):
            typ = typ.schemas[1]
        if isinstance(typ, avro.schema.MapSchema):
            ret = ["string"]
        elif isinstance(typ, avro.schema.ArraySchema):
            ret = 42.0
        elif isinstance(typ, avro.schema.EnumSchema):
            ret = "NO ENUM COULD HAVE THIS NAME"
        elif isinstance(typ, avro.schema.RecordSchema):
            ret = 42
        elif typ.type in self.badValueMap:
            ret = self.badValueMap[typ.type]
        else:
            raise Exception("schema assumptions violated")
        return ret

    def getRandomValue(self, cls, fieldName):
        """
        Returns a value randomly generated conforming to the schema.
        """
        maxListSize = 10
        field = self.getAvroSchema(cls, fieldName)
        typ = field.type
        if isinstance(typ, avro.schema.UnionSchema):
            # TODO put in some probability that return None if
            # we have a union schema
            typ = typ.schemas[1]
        if isinstance(typ, avro.schema.MapSchema):
            ret = {}
            for j in range(random.randint(0, maxListSize)):
                randStr = [randomString() for _ in range(
                    random.randint(0, maxListSize))]
                ret["key{0}".format(j)] = randStr
        elif isinstance(typ, avro.schema.ArraySchema):
            randInt = random.randint(0, maxListSize)
            if cls.isEmbeddedType(field.name):
                embeddedClass = cls.getEmbeddedType(field.name)

                def getRandInst():
                    return self.getRandomInstance(embeddedClass)
            else:
                getRandInst = self.randomValueMap[typ.items.type]
            ret = [getRandInst() for j in range(randInt)]
        elif isinstance(typ, avro.schema.EnumSchema):
            ret = random.choice(typ.symbols)
        elif isinstance(typ, avro.schema.RecordSchema):
            self.assertTrue(cls.isEmbeddedType(fieldName))
            embeddedClass = cls.getEmbeddedType(fieldName)
            ret = self.getRandomInstance(embeddedClass)
        elif typ.type in self.randomValueMap:
            ret = self.randomValueMap[typ.type]()
        else:
            raise Exception("schema assumptions violated")
        return ret

    def getTypicalValue(self, cls, fieldName):
        """
        Returns a typical value for the specified field on the specified
        Protocol class.
        """
        # We make some simplifying assumptions about how the schema is
        # structured which fits the way the GA4GH protocol is currently
        # designed but may break in the future. We try to at least flag
        # this fact here.
        err = "Schema structure assumptions violated"
        field = self.getAvroSchema(cls, fieldName)
        typ = field.type
        if isinstance(typ, avro.schema.UnionSchema):
            t0 = typ.schemas[0]
            if (isinstance(t0, avro.schema.PrimitiveSchema)
                    and t0.type == "null"):
                typ = typ.schemas[1]
            else:
                raise Exception(err)
        ret = None
        if isinstance(typ, avro.schema.MapSchema):
            ret = {"key": ["value1", "value2"]}
            if not isinstance(typ.values, avro.schema.ArraySchema):
                raise Exception(err)
            if not isinstance(typ.values.items, avro.schema.PrimitiveSchema):
                raise Exception(err)
            if typ.values.items.type != "string":
                raise Exception(err)
        elif isinstance(typ, avro.schema.ArraySchema):
            if cls.isEmbeddedType(field.name):
                embeddedClass = cls.getEmbeddedType(field.name)
                ret = [self.getTypicalInstance(embeddedClass)]
            else:
                if not isinstance(typ.items, avro.schema.PrimitiveSchema):
                    raise Exception(err)
                ret = [self.typicalValueMap[typ.items.type]]
        elif isinstance(typ, avro.schema.EnumSchema):
            ret = typ.symbols[0]
        elif isinstance(typ, avro.schema.RecordSchema):
            self.assertTrue(cls.isEmbeddedType(fieldName))
            embeddedClass = cls.getEmbeddedType(fieldName)
            ret = self.getTypicalInstance(embeddedClass)
        elif typ.type in self.typicalValueMap:
            ret = self.typicalValueMap[typ.type]
        else:
            raise Exception("schema assumptions violated")
        return ret

    def getTypicalInstance(self, cls):
        """
        Returns a typical instance of the specified protocol class.
        """
        instance = cls()
        for field in cls.schema.fields:
            setattr(instance, field.name,
                    self.getTypicalValue(cls, field.name))
        return instance

    def getRandomInstance(self, cls):
        """
        Returns an instance of the specified class with randomly generated
        values conforming to the schema.
        """
        instance = cls()
        for field in cls.schema.fields:
            setattr(instance, field.name,
                    self.getRandomValue(cls, field.name))
        return instance

    def setRequiredValues(self, instance):
        """
        Sets the required values in the specified instance to typical values.
        """
        for key in instance.__slots__:
            if key in instance.requiredFields:
                value = self.getTypicalValue(type(instance), key)
                setattr(instance, key, value)

    def getDefaultInstance(self, cls):
        """
        Returns a new instance with the required values set.
        """
        instance = cls()
        self.setRequiredValues(instance)
        return instance


class EqualityTest(SchemaTest):
    """
    Tests equality is correctly calculated for different protocol elements.
    """
    def verifyEqualityOperations(self, i1, i2):
        self.assertEqual(i1, i1)
        self.assertTrue(i1 == i1)
        self.assertFalse(i1 != i1)
        self.assertEqual(i1, i2)
        self.assertTrue(i1 == i2)
        self.assertFalse(i1 != i2)
        for val in [None, {}, [], object, ""]:
            self.assertFalse(i1 == val)
            self.assertTrue(i1 != val)
            self.assertFalse(val == i1)
            self.assertTrue(val != i1)
        # Now change an attribute on one and check if equality fails.
        for field in i1.schema.fields:
            setattr(i1, field.name, "value unique to i1")
            self.assertFalse(i1 == i2)
            self.assertTrue(i1 != i2)

    def testSameClasses(self):
        factories = [self.getDefaultInstance, self.getTypicalInstance,
                     self.getRandomInstance]
        for cls in self.getProtocolClasses():
            for factory in factories:
                i1 = factory(cls)
                i2 = cls.fromJsonDict(i1.toJsonDict())
                self.verifyEqualityOperations(i1, i2)

    def testDifferentValues(self):
        def factory(cls):
            return cls()
        factories = [factory, self.getTypicalInstance, self.getDefaultInstance,
                     self.getRandomInstance]
        classes = list(self.getProtocolClasses())
        c1 = classes[0]
        for c2 in classes[1:]:
            for factory in factories:
                i1 = factory(c1)
                i2 = factory(c2)
                self.assertFalse(i1 == i2)
                self.assertTrue(i1 != i2)

    def testDifferentLengthArrays(self):
        i1 = self.getTypicalInstance(protocol.GACallSet)
        i2 = protocol.GACallSet.fromJsonDict(i1.toJsonDict())
        i2.variantSetIds.append("extra")
        self.assertFalse(i1 == i2)


class SerialisationTest(SchemaTest):
    """
    Tests the serialisation and deserialisation code for the schema classes
    """
    def validateClasses(self, factory):
        for cls in self.getProtocolClasses():
            instance = factory(cls)
            jsonStr = instance.toJsonString()
            otherInstance = cls.fromJsonString(jsonStr)
            self.assertEqual(instance, otherInstance)

            jsonDict = instance.toJsonDict()
            otherInstance = cls.fromJsonDict(jsonDict)
            self.assertEqual(instance, otherInstance)

    def testSerialiseDefaultValues(self):
        self.validateClasses(self.getDefaultInstance)

    def testSerialiseTypicalValues(self):
        self.validateClasses(self.getTypicalInstance)

    def testSerialiseRandomValues(self):
        self.validateClasses(self.getRandomInstance)


class ValidatorTest(SchemaTest):
    """
    Tests the validator to see if it will correctly identify instances
    that do not match the schema and also that it correctly identifies
    instances that do match the schema
    """
    def validateClasses(self, factory):
        for cls in self.getProtocolClasses():
            instance = factory(cls)
            jsonDict = instance.toJsonDict()
            self.assertTrue(cls.validate(jsonDict))

    def testValidateDefaultValues(self):
        self.validateClasses(self.getDefaultInstance)

    def testValidateTypicalValues(self):
        self.validateClasses(self.getTypicalInstance)

    def testValidateRandomValues(self):
        self.validateClasses(self.getRandomInstance)

    def testValidateBadValues(self):
        for cls in self.getProtocolClasses():
            instance = self.getTypicalInstance(cls)
            jsonDict = instance.toJsonDict()
            self.assertFalse(cls.validate(None))
            self.assertFalse(cls.validate([]))
            self.assertFalse(cls.validate(1))
            # setting values to bad values should be invalid
            for key in jsonDict.keys():
                dct = dict(jsonDict)
                dct[key] = self.getInvalidValue(cls, key)
                self.assertFalse(cls.validate(dct))
            for c in tests.powerset(jsonDict.keys(), 10):
                if len(c) > 0:
                    dct = dict(jsonDict)
                    for f in c:
                        dct[f] = self.getInvalidValue(cls, f)
                    self.assertFalse(cls.validate(dct))
