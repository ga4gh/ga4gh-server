import sqlite3
import uuid
import json

import ga4gh.exceptions as exceptions

import datetime

class SQLBackend():
    # TODO generalize for the rest of the API?

    def __init__(self, objClass, dbFile, nestedFields=["info"]):
        self.nestedFields = nestedFields
        self.objClass = objClass
        self.dbFile = dbFile
        blankObj = objClass()
        self.columns = [str(k) for k in blankObj.toJsonDict()]
        self.tableName = type(blankObj).__name__
        # TODO protobuf will have different types so
        # convert them to SQLite ones here
        self.coltypes = ("text" for k in blankObj.toJsonDict())
        self.safeMakeTable()

    def insertString(self):
        """
        Creates an insert statement with length equal to
        the number of columns for the object being stored.

        :return: sql insert statement
        """
        sql = "INSERT INTO " + self.tableName + " VALUES("
        for c in self.columns:
            sql += "?, "
        sql = sql.rstrip(", ")
        sql += ");"
        return sql

    def orderedValuesFromDict(self, obj):
        """
        Receives a protocol object and converts its values
        into an ordered list in accordance with the columns
        for writing insert statements.

        Nested fields are coerced into a JSON string.
        :param obj:
        :return: ordered list of values
        """
        # TODO can coerce types here
        ret = []
        for k in self.columns:
            if k not in self.nestedFields:
                ret.append(str(obj.__getattribute__(k)))
            else:
                ret.append(str(json.dumps(obj.__getattribute__(k))))
        return ret

    def objToClassObj(self, obj):
        """
        Converts a plain old python dict and coerces it into
        the class object in self.objClass().
        :param obj:
        :return: class object of type self.objClass()
        """
        classObj = self.objClass()
        for k in classObj.toJsonDict():
            if k in obj:
                classObj.__setattr__(k, obj[k])
        return classObj

    def safeMakeTable(self):
        """
        Creates a table in the self.dbFile according to the
        class object stored.
        :return:
        """
        # TODO fail when there is a table with different cols
        self._dbconn = sqlite3.connect(self.dbFile)
        dbconn = self._dbconn
        dbcur = dbconn.cursor()
        dbcur.execute(self.tableString())
        dbcur.close()
        return

    def tableString(self):
        """
        Assembles the `create table` statement out of
        self.tableName, self.columns, and self.coltypes.
        :return: sql string for table creation
        """
        string = "create table if not exists " + self.tableName
        string += "("
        for c in zip(self.columns, self.coltypes):
            string += c[0] + " " + c[1] + ", "
        string = string.rstrip(", ")
        string += ")"
        return string

    def hydrate(self, row):
        """
        Assembles a row into a class object and returns
        that class instance.
        :param row:
        :return: object of class self.objClass
        """
        hydrated = dict(zip(self.columns, row))
        for k in hydrated:
            if k in self.nestedFields:
                hydrated[k] = json.loads(hydrated[k])
        return self.objToClassObj(hydrated)

    def getById(self, id):
        """
        Gets a single row by ID and returns a class
        instance for the given row.
        :param id:
        :return: object of class self.objClass
        """
        # TODO sanitize inputs
        sql = "SELECT * FROM " + self.tableName + " where id=='" + id + "';"
        query = self._dbconn.execute(sql)
        results = query.fetchone()
        if results:
            return self.hydrate(results)
        else:
            return None

    def delete(self, obj):
        # TODO
        pass

    def insert(self, obj):
        """
        Performs SQL insertion of a given object and returns
        the newly created object with ID field populated.
        :param obj:
        :return: obj representing new row
        """
        # generate ID
        obj = self.objToClassObj(obj)
        obj.id = str(uuid.uuid4())
        # TODO make sure it actually has this field (created/updated)
        obj.created = datetime.datetime.now().isoformat()
        obj.updated = obj.created
        sql = self.insertString()
        query = self._dbconn.execute(sql, self.orderedValuesFromDict(obj))
        # TODO add message using cursor dereferencing
        self._dbconn.commit()
        return self.getById(obj.id)

    def updateString(self, obj):
        """
        Assembles an SQL statement for performing an update given
        an object with new fields.
        :param obj:
        :return: sql string that performs update
        """
        # TODO sanitize inputs and type coercion
        sql = "UPDATE " + self.tableName + " SET "
        for c in self.columns:
            if c in self.nestedFields:
                sql += c + " = '" + str(obj.__getattribute__(c)) + "', "
            elif obj.__getattribute__(c):
                sql += c + " = '" + obj.__getattribute__(c) + "', "
        sql = sql.rstrip(", ")
        sql += " WHERE id=='" + obj.id + "'"
        return sql

    def update(self, obj):
        """
        Performs SQL update of an object using its ID.
        :param obj:
        :return: Newly updated object.
        """
        # get the record we want to update
        obj = self.objToClassObj(obj)
        old = self.getById(obj.id)
        if not old:
            # there should be an existing record to update
            raise exceptions.BadRequestException
        # overwrite dict values with new
        for k in old.toJsonDict():
            if k not in obj.toJsonDict():
                obj.__setattr__(k, old.__getattribute__(k))
        sql = self.updateString(obj)
        print(sql)
        query = self._dbconn.execute(sql)
        print(query)
        # TODO add message using cursor dereferencing
        self._dbconn.commit()
        # update SQL
        return self.getById(obj.id)