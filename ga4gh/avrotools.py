"""
Provides additional functionality based on avro
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import avro.io

import ga4gh.protocol as protocol


class ValidationToolException(Exception):
    """
    Something that went wrong in the ValidationTool
    """


class ValidationTool(object):
    """
    Provides additional validation methods based on avro
    """
    @classmethod
    def getInvalidFields(cls, class_, jsonDict):
        """
        Takes a ProtocolElement subclass class_ and a jsonDict instantiation
        of that class.  Returns a dictionary whose keys are the names of
        fields in the jsonDict that failed validation.  The values of the
        returned dictionary are the jsonDict's values for those fields.

        This is an expensive function to run and shouldn't be on the main
        data path.

        TODO is there an easier/faster way to do this than generating
        an on-demand schema for every field we want to test?
        """
        if not issubclass(class_, protocol.ProtocolElement):
            message = "class '{}' is not subclass of ProtocolElement".format(
                class_.__name__)
            raise ValidationToolException(message)
        invalidFields = {}
        schema = class_.schema
        schemaJson = schema.to_json()
        schemaJsonFields = schemaJson['fields']
        cls._ensureNoExtraFields(schemaJsonFields, jsonDict)
        for field in schemaJsonFields:
            fieldName = field['name']
            try:
                fieldValue = jsonDict[fieldName]
            except KeyError:
                message = "jsonDict does not contain field '{}' for " \
                    "class '{}'".format(fieldName, class_.__name__)
                raise ValidationToolException(message)
            fieldSchema = cls._createFieldSchema(schema, fieldName)
            fieldDict = {fieldName: fieldValue}
            isFieldValid = avro.io.validate(fieldSchema, fieldDict)
            if not isFieldValid:
                invalidFields[fieldName] = fieldValue
        return invalidFields

    @classmethod
    def _ensureNoExtraFields(cls, schemaJsonFields, jsonDict):
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

    @classmethod
    def _createFieldSchema(cls, schema, fieldName):
        schemaJson = schema.to_json()
        # TODO this is an ugly hack
        # I think GAReadAlignment.nextMatePosition's schema is not getting
        # written correctly.
        # In the schema its type attribute is:
        #   "type": ["null", "GAPosition"]
        # which looks wrong.  I think it should be:
        #   "type": ["null", {...schema record of GAPosition...}]
        # which looks more similar to other attributes with embedded types,
        # like GAReadAlignment.alignment
        if fieldName == "nextMatePosition":
            new = cls._createNewFieldsValue(schemaJson, fieldName)
            gaPositionSchemaSource = \
                json.loads(protocol.GAPosition._schemaSource)
            new['type'] = [u'null', gaPositionSchemaSource]
            schemaJson['fields'] = [new]
        else:
            schemaJson['fields'] = [cls._createNewFieldsValue(
                schemaJson, fieldName)]
        fieldSchemaSource = json.dumps(schemaJson)
        fieldSchema = avro.schema.parse(fieldSchemaSource)
        return fieldSchema

    @classmethod
    def _createNewFieldsValue(cls, schemaJson, fieldName):
        for fieldSchemaDict in schemaJson['fields']:
            if fieldSchemaDict['name'] == fieldName:
                return fieldSchemaDict
