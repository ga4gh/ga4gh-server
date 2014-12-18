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
    n = random.randint(0, 10)
    s = ''.join(random.choice(string.ascii_letters) for x in range(n))
    return s


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
        for f in cls.schema.fields:
            if f.name == fieldName:
                field = f
        return field

    def getInvalidValue(self, cls, fieldName):
        """
        Returns a value that should trigger a schema validation failure.
        """
        field = self.getAvroSchema(cls, fieldName)
        t = field.type
        if isinstance(t, avro.schema.UnionSchema):
            t = t.schemas[1]
        if isinstance(t, avro.schema.MapSchema):
            ret = ["string"]
        elif isinstance(t, avro.schema.ArraySchema):
            ret = 42.0
        elif isinstance(t, avro.schema.EnumSchema):
            ret = "NO ENUM COULD HAVE THIS NAME"
        elif isinstance(t, avro.schema.RecordSchema):
            ret = 42
        elif t.type in self.badValueMap:
            ret = self.badValueMap[t.type]
        else:
            raise Exception("schema assumptions violated")
        return ret

    def getRandomValue(self, cls, fieldName):
        """
        Returns a value randomly generated conforming to the schema.
        """
        maxListSize = 10
        field = self.getAvroSchema(cls, fieldName)
        t = field.type
        if isinstance(t, avro.schema.UnionSchema):
            # TODO put in some probability that return None if
            # we have a union schema
            t = t.schemas[1]
        if isinstance(t, avro.schema.MapSchema):
            ret = {}
            for j in range(random.randint(0, maxListSize)):
                l = [randomString() for j in range(
                    random.randint(0, maxListSize))]
                ret["key{0}".format(j)] = l
        elif isinstance(t, avro.schema.ArraySchema):
            n = random.randint(0, maxListSize)
            if cls.isEmbeddedType(field.name):
                embeddedClass = cls.getEmbeddedType(field.name)

                def f():
                    return self.getRandomInstance(embeddedClass)
            else:
                f = self.randomValueMap[t.items.type]
            ret = [f() for j in range(n)]
        elif isinstance(t, avro.schema.EnumSchema):
            ret = random.choice(t.symbols)
        elif isinstance(t, avro.schema.RecordSchema):
            assert cls.isEmbeddedType(fieldName)
            embeddedClass = cls.getEmbeddedType(fieldName)
            ret = self.getRandomInstance(embeddedClass)
        elif t.type in self.randomValueMap:
            ret = self.randomValueMap[t.type]()
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
        t = field.type
        if isinstance(t, avro.schema.UnionSchema):
            t0 = t.schemas[0]
            if (isinstance(t0, avro.schema.PrimitiveSchema)
                    and t0.type == "null"):
                t = t.schemas[1]
            else:
                raise Exception(err)
        ret = None
        if isinstance(t, avro.schema.MapSchema):
            ret = {"key": ["value1", "value2"]}
            if not isinstance(t.values, avro.schema.ArraySchema):
                raise Exception(err)
            if not isinstance(t.values.items, avro.schema.PrimitiveSchema):
                raise Exception(err)
            if t.values.items.type != "string":
                raise Exception(err)
        elif isinstance(t, avro.schema.ArraySchema):
            if cls.isEmbeddedType(field.name):
                embeddedClass = cls.getEmbeddedType(field.name)
                ret = [self.getTypicalInstance(embeddedClass)]
            else:
                if not isinstance(t.items, avro.schema.PrimitiveSchema):
                    raise Exception(err)
                ret = [self.typicalValueMap[t.items.type]]
        elif isinstance(t, avro.schema.EnumSchema):
            ret = t.symbols[0]
        elif isinstance(t, avro.schema.RecordSchema):
            assert cls.isEmbeddedType(fieldName)
            embeddedClass = cls.getEmbeddedType(fieldName)
            ret = self.getTypicalInstance(embeddedClass)
        elif t.type in self.typicalValueMap:
            ret = self.typicalValueMap[t.type]
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
        for k in instance.__dict__.keys():
            if k in instance.requiredFields:
                v = self.getTypicalValue(type(instance), k)
                setattr(instance, k, v)

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
        for v in [None, {}, [], object, ""]:
            self.assertFalse(i1 == v)
            self.assertTrue(i1 != v)
            self.assertFalse(v == i1)
            self.assertTrue(v != i1)
        # Now change an attribute on one and check if equality fails.
        for f in i1.schema.fields:
            setattr(i1, f.name, "value unique to i1")
            self.assertFalse(i1 == i2)
            self.assertTrue(i1 != i2)

    def testSameClasses(self):
        factories = [self.getDefaultInstance, self.getTypicalInstance,
                     self.getRandomInstance]
        for cls in self.getProtocolClasses():
            for f in factories:
                i1 = f(cls)
                i2 = cls.fromJSON(i1.toJSON())
                self.verifyEqualityOperations(i1, i2)

    def testDifferentValues(self):
        f = lambda cls: cls()
        factories = [f, self.getTypicalInstance, self.getDefaultInstance,
                     self.getRandomInstance]
        classes = list(self.getProtocolClasses())
        c1 = classes[0]
        for c2 in classes[1:]:
            for f in factories:
                i1 = f(c1)
                i2 = f(c2)
                self.assertFalse(i1 == i2)
                self.assertTrue(i1 != i2)


class SerialisationTest(SchemaTest):
    """
    Tests the serialisation and deserialisation code for the schema classes
    """
    def validateClasses(self, factory):
        for cls in self.getProtocolClasses():
            instance = factory(cls)
            jsonStr = instance.toJSON()
            otherInstance = cls.fromJSON(jsonStr)
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
            jsonDict = instance.toJSONDict()
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
            jsonDict = instance.toJSONDict()
            self.assertFalse(cls.validate(None))
            self.assertFalse(cls.validate([]))
            self.assertFalse(cls.validate(1))
            # setting values to bad values should be invalid
            for f in jsonDict.keys():
                d = dict(jsonDict)
                d[f] = self.getInvalidValue(cls, f)
                self.assertFalse(cls.validate(d))
            for c in tests.powerset(jsonDict.keys(), 10):
                if len(c) > 0:
                    d = dict(jsonDict)
                    for f in c:
                        d[f] = self.getInvalidValue(cls, f)
                    self.assertFalse(cls.validate(d))
