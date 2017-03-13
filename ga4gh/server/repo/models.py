"""
peewee is a lightweight ORM with SQLite, postgresql,
and MySQL support. This file presents models for the
registry database.

Partially auto-generated using pwiz.

    python -m pwiz -e sqlite ga4gh-example-data/registry.db > models.py

For more on the peewee model API see:

https://peewee.readthedocs.io/en/latest/peewee/models.html

"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import datetime
import peewee as pw

# The databaseProxy is used to dynamically changed the
# backing database and needs to be set to an actual
# database instance to use these models.
databaseProxy = pw.Proxy()


class SqliteDatabase(pw.SqliteDatabase):
    def __init__(self, *_, **__):
        super(SqliteDatabase, self).__init__(*_, **__)


class BaseModel(pw.Model):
    attributes = pw.TextField(null=True)

    class Meta:
        database = databaseProxy


class Peer(BaseModel):
    url = pw.TextField(primary_key=True, null=False, unique=True)


class Announcement(BaseModel):
    # Provides the storage layer for AnnouncePeerRequest
    # Primary key for the record autoincrement
    id = pw.PrimaryKeyField()
    # URL submitted as a possible peer
    url = pw.TextField(null=False)
    # Other information about the request stored as JSON.
    remote_addr = pw.TextField(null=True)
    user_agent = pw.TextField(null=True)
    # The time at which this record was created.
    created = pw.DateTimeField()

    def save(self, *args, **kwargs):
        self.created = datetime.datetime.now()
        return super(Announcement, self).save(*args, **kwargs)


class Dataset(BaseModel):
    description = pw.TextField(null=True)
    id = pw.TextField(primary_key=True)
    info = pw.TextField(null=True)
    name = pw.TextField(unique=True)


class Biosample(BaseModel):
    created = pw.TextField(null=True)
    datasetid = pw.ForeignKeyField(
        db_column='datasetId', rel_model=Dataset, to_field='id')
    description = pw.TextField(null=True)
    disease = pw.TextField(null=True)
    id = pw.TextField(primary_key=True)
    individualid = pw.TextField(db_column='individualId', null=True)
    info = pw.TextField(null=True)
    name = pw.TextField()
    updated = pw.TextField(null=True)
    individualAgeAtCollection = pw.TextField(null=True)

    class Meta:
        indexes = (
            (('datasetid', 'name'), True),
        )


class Referenceset(BaseModel):
    assemblyid = pw.TextField(db_column='assemblyId', null=True)
    dataurl = pw.TextField(db_column='dataUrl')
    description = pw.TextField(null=True)
    id = pw.TextField(primary_key=True)
    isderived = pw.IntegerField(db_column='isDerived', null=True)
    md5checksum = pw.TextField(null=True)
    name = pw.TextField(unique=True)
    species = pw.TextField(db_column='species', null=True)
    sourceaccessions = pw.TextField(db_column='sourceAccessions', null=True)
    sourceuri = pw.TextField(db_column='sourceUri', null=True)


class Variantset(BaseModel):
    created = pw.TextField(null=True)
    dataurlindexmap = pw.TextField(db_column='dataUrlIndexMap')
    datasetid = pw.ForeignKeyField(
        db_column='datasetId', rel_model=Dataset, to_field='id')
    id = pw.TextField(primary_key=True)
    metadata = pw.TextField(null=True)
    name = pw.TextField()
    referencesetid = pw.ForeignKeyField(
        db_column='referenceSetId', rel_model=Referenceset, to_field='id')
    updated = pw.TextField(null=True)

    class Meta:
        indexes = (
            (('datasetid', 'name'), True),
        )


class Callset(BaseModel):
    biosampleid = pw.TextField(db_column='biosampleId', null=True)
    id = pw.TextField(primary_key=True)
    name = pw.TextField()
    variantsetid = pw.ForeignKeyField(
        db_column='variantSetId', rel_model=Variantset, to_field='id')

    class Meta:
        indexes = (
            (('variantsetid', 'name'), True),
        )


class Ontology(BaseModel):
    dataurl = pw.TextField(db_column='dataUrl')
    id = pw.TextField(primary_key=True)
    name = pw.TextField(unique=True)
    ontologyprefix = pw.TextField(db_column='ontologyPrefix')


class Featureset(BaseModel):
    dataurl = pw.TextField(db_column='dataUrl')
    datasetid = pw.ForeignKeyField(
        db_column='datasetId', rel_model=Dataset, to_field='id')
    id = pw.TextField(primary_key=True)
    info = pw.TextField(null=True)
    name = pw.TextField()
    ontologyid = pw.ForeignKeyField(
        db_column='ontologyId', rel_model=Ontology, to_field='id')
    referencesetid = pw.ForeignKeyField(
        db_column='referenceSetId', rel_model=Referenceset, to_field='id')
    sourceuri = pw.TextField(
        db_column='sourceUri', null=True)

    class Meta:
        indexes = (
            (('datasetid', 'name'), True),
        )


class ContinuousSet(BaseModel):
    dataurl = pw.TextField(db_column='dataUrl')
    datasetid = pw.ForeignKeyField(
        db_column='datasetId', rel_model=Dataset, to_field='id')
    id = pw.TextField(primary_key=True)
    info = pw.TextField(null=True)
    name = pw.TextField()
    referencesetid = pw.ForeignKeyField(
        db_column='referenceSetId', rel_model=Referenceset, to_field='id')
    sourceuri = pw.TextField(
        db_column='sourceUri', null=True)

    class Meta:
        indexes = (
            (('datasetid', 'name'), True),
        )


class Individual(BaseModel):
    created = pw.TextField()
    datasetid = pw.ForeignKeyField(
        db_column='datasetId', rel_model=Dataset, to_field='id')
    description = pw.TextField(null=True)
    id = pw.TextField(primary_key=True)
    info = pw.TextField(null=True)
    name = pw.TextField(null=True)
    sex = pw.TextField(null=True)
    species = pw.TextField(null=True)
    updated = pw.TextField(null=True)

    class Meta:
        indexes = (
            (('datasetid', 'name'), True),
        )


class Phenotypeassociationset(BaseModel):
    dataurl = pw.TextField(db_column='dataUrl')
    datasetid = pw.ForeignKeyField(
        db_column='datasetId', rel_model=Dataset, to_field='id')
    id = pw.TextField(primary_key=True)
    name = pw.TextField(null=True)

    class Meta:
        indexes = (
            (('datasetid', 'name'), True),
        )


class Readgroupset(BaseModel):
    dataurl = pw.TextField(db_column='dataUrl')
    datasetid = pw.ForeignKeyField(
        db_column='datasetId', rel_model=Dataset, to_field='id')
    id = pw.TextField(primary_key=True)
    indexfile = pw.TextField(db_column='indexFile')
    name = pw.TextField()
    programs = pw.TextField(null=True)
    referencesetid = pw.ForeignKeyField(
        db_column='referenceSetId', rel_model=Referenceset, to_field='id')
    stats = pw.TextField()

    class Meta:
        indexes = (
            (('datasetid', 'name'), True),
        )


class Readgroup(BaseModel):
    biosampleid = pw.TextField(db_column='biosampleId', null=True)
    created = pw.TextField(null=True)
    description = pw.TextField(null=True)
    experiment = pw.TextField()
    id = pw.TextField(primary_key=True)
    name = pw.TextField()
    predictedinsertsize = pw.IntegerField(
        db_column='predictedInsertSize', null=True)
    readgroupsetid = pw.ForeignKeyField(
        db_column='readGroupSetId', rel_model=Readgroupset, to_field='id')
    samplename = pw.TextField(db_column='sampleName', null=True)
    stats = pw.TextField()
    updated = pw.TextField(null=True)

    class Meta:
        indexes = (
            (('readgroupsetid', 'name'), True),
        )


class Reference(BaseModel):
    id = pw.TextField(null=True, primary_key=True)
    isderived = pw.IntegerField(db_column='isDerived', null=True)
    length = pw.IntegerField(null=True)
    md5checksum = pw.TextField(null=True)
    name = pw.TextField()
    species = pw.TextField(db_column='species', null=True)
    referencesetid = pw.ForeignKeyField(
        db_column='referenceSetId', rel_model=Referenceset, to_field='id')
    sourceaccessions = pw.TextField(db_column='sourceAccessions', null=True)
    sourcedivergence = pw.FloatField(db_column='sourceDivergence', null=True)
    sourceuri = pw.TextField(db_column='sourceUri', null=True)

    class Meta:
        indexes = (
            (('referencesetid', 'name'), True),
        )


class Rnaquantificationset(BaseModel):
    dataurl = pw.TextField(db_column='dataUrl')
    datasetid = pw.ForeignKeyField(
        db_column='datasetId', rel_model=Dataset, to_field='id')
    id = pw.TextField(primary_key=True)
    info = pw.TextField(null=True)
    name = pw.TextField()
    referencesetid = pw.ForeignKeyField(
        db_column='referenceSetId', rel_model=Referenceset, to_field='id')

    class Meta:
        indexes = (
            (('datasetid', 'name'), True),
        )


class System(BaseModel):
    key = pw.TextField(primary_key=True)
    value = pw.TextField()


class Variantannotationset(BaseModel):
    analysis = pw.TextField(null=True)
    annotationtype = pw.TextField(db_column='annotationType', null=True)
    created = pw.TextField(null=True)
    id = pw.TextField(primary_key=True)
    name = pw.TextField()
    ontologyid = pw.ForeignKeyField(
        db_column='ontologyId', rel_model=Ontology, to_field='id')
    updated = pw.TextField(null=True)
    variantsetid = pw.ForeignKeyField(
        db_column='variantSetId', rel_model=Variantset, to_field='id')

    class Meta:
        indexes = (
            (('variantsetid', 'name'), True),
        )
