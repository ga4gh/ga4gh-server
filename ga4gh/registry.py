"""
The registry database for all the top-level objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import random
import hashlib
import datetime
import json
import base64

import sqlalchemy
import sqlalchemy.exc
import sqlalchemy.orm as orm
import sqlalchemy.ext.declarative

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


def get_external_id():
    # TODO make this use SystemRandom and make a table specific prefix.
    return random.randint(0, 2**31)


def create_id_column():
    return sqlalchemy.Column(
        sqlalchemy.Integer, primary_key=True, autoincrement=False,
        default=get_external_id)


def create_timestamp_column():
    # TODO we should add a update default here as well.
    return sqlalchemy.Column(
        sqlalchemy.DateTime, nullable=False,
        default=datetime.datetime.now())


SqlAlchemyBase = sqlalchemy.ext.declarative.declarative_base()


#####################################################################
#
# Datasets and tables to represent schema wide objects.
#
#####################################################################


class Dataset(SqlAlchemyBase):
    __tablename__ = 'Dataset'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False, unique=True)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    # We may have a large number of child collections per dataset,
    # so we don't want to populate these lists unless we need them.
    variant_sets = orm.relationship(
        "VariantSet", back_populates="dataset", lazy="dynamic",
        cascade="all, delete, delete-orphan")
    read_group_sets = orm.relationship(
        "ReadGroupSet", back_populates="dataset", lazy="dynamic",
        cascade="all, delete, delete-orphan")
    feature_sets = orm.relationship(
        "FeatureSet", back_populates="dataset", lazy="dynamic",
        cascade="all, delete, delete-orphan")
    individuals = orm.relationship(
        "Individual", back_populates="dataset", lazy="dynamic",
        cascade="all, delete, delete-orphan")
    bio_samples = orm.relationship(
        "BioSample", back_populates="dataset", lazy="dynamic",
        cascade="all, delete, delete-orphan")

    def __init__(self, name):
        self.name = name
        self.description = ""

    def get_protobuf(self):
        ret = protocol.Dataset()
        ret.id = str(self.id)
        ret.name = self.name
        ret.description = self.description
        return ret


class BioSample(SqlAlchemyBase):
    """
    A BioSample refers to a unit of biological material from which the
    substrate molecules (e.g. genomic DNA, RNA, proteins) for molecular
    analyses (e.g. sequencing, array hybridisation, mass-spectrometry) are
    extracted. Examples would be a tissue biopsy, a single cell from a culture
    for single cell genome sequencing or a protein fraction from a gradient
    centrifugation. Several instances (e.g. technical replicates) or types of
    experiments (e.g. genomic array as well as RNA-seq experiments) may refer
    to the same BioSample. In the context of the GA4GH metadata schema,
    BioSample constitutes the central reference object.
    """
    __tablename__ = "BioSample"

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    disease = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()

    individual_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Individual.id"),
        nullable=False)
    individual = orm.relationship(
        "Individual", back_populates="bio_samples", single_parent=True)
    dataset_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Dataset.id"),
        nullable=False)
    dataset = orm.relationship(
        "Dataset", back_populates="bio_samples", single_parent=True)

    # TODO Info

    __table_args__ = (
        # BioSample names must be unique within a dataset
        sqlalchemy.UniqueConstraint("dataset_id", "name"),
    )

    def __init__(self, name):
        self.name = name
        self.description = ""
        self.disease = ""

    def get_protobuf(self):
        ret = protocol.BioSample()
        ret.id = str(self.id)
        ret.name = self.name
        ret.description = self.description
        ret.individual_id = str(self.individual_id)
        ret.dataset_id = str(self.dataset_id)
        # .message_create_time = protocol.datetime_to_iso8601(
        #     self.creation_timestamp)
        # experiment.message_update_time = protocol.datetime_to_iso8601(
        #     self.update_timestamp)
        # TODO disease
        return ret


class Individual(SqlAlchemyBase):
    """
    An individual (or subject) typically corresponds to an individual
    human or another organism.
    """
    __tablename__ = "Individual"

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    sex = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()

    bio_samples = orm.relationship(
        "BioSample", back_populates="individual",
        cascade="all, delete, delete-orphan")
    dataset_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Dataset.id"),
        nullable=False)
    dataset = orm.relationship(
        "Dataset", back_populates="individuals", single_parent=True)
    # TODO Info

    __table_args__ = (
        # Individual names must be unique within a dataset
        sqlalchemy.UniqueConstraint("dataset_id", "name"),
    )

    def __init__(self, name):
        self.name = name
        self.description = ""
        self.sex = ""

    def get_protobuf(self):
        ret = protocol.Individual()
        ret.id = str(self.id)
        ret.name = self.name
        ret.description = self.description
        ret.dataset_id = str(self.dataset_id)
        # .message_create_time = protocol.datetime_to_iso8601(
        #     self.creation_timestamp)
        # experiment.message_update_time = protocol.datetime_to_iso8601(
        #     self.update_timestamp)
        # TODO sex
        return ret


class Info(SqlAlchemyBase):
    """
    A generic key-value table used to store info maps used in various
    places in the schema.
    """
    # TODO complete this. We probably need a more complicated
    # relationship because we can multiple values per key.
    __tablename__ = 'Info'

    id = sqlalchemy.Column(
        sqlalchemy.Integer, primary_key=True, autoincrement=True)
    key = sqlalchemy.Column(sqlalchemy.String, nullable=False, unique=True)
    value = sqlalchemy.Column(sqlalchemy.String, nullable=False)


#####################################################################
#
# Ontologies
#
#####################################################################

class Ontology(SqlAlchemyBase):
    __tablename__ = 'Ontology'

    id = create_id_column()
    """
    The server ID for this Ontology.
    """

    name = sqlalchemy.Column(sqlalchemy.String, nullable=False, unique=True)
    """
    The name of this ontology. Unique within the repository.
    """

    ontology_prefix = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    The ontology prefix string, i.e. "SO" for a sequence ontology
    and "GO" for gene a ontology.
    """

    # TODO Should we add more fields here to (hopefully) avoid having to
    # upgrade the tables in the future??

    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity': 'Ontology',
         'polymorphic_on': type
    }


#####################################################################
#
# References
#
#####################################################################


class Accession(SqlAlchemyBase):
    """
    Accession idenfiers for sequences.
    """
    __tablename__ = "Accession"

    id = sqlalchemy.Column(
        sqlalchemy.Integer, primary_key=True, autoincrement=True)

    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    The accession ID string.
    """

    def __init__(self, name):
        self.name = name

    __table_args__ = (
        # Accession names must be unique.
        sqlalchemy.UniqueConstraint("name"),
    )


# There is a many-to-many association between References and Accessions,
# so we need an association table.
_reference_accessions_table = sqlalchemy.Table(
    'reference_accession_association', SqlAlchemyBase.metadata,
    sqlalchemy.Column(
        'reference_id', sqlalchemy.Integer,
        sqlalchemy.ForeignKey('Reference.id')),
    sqlalchemy.Column(
        'accession_id', sqlalchemy.Integer,
        sqlalchemy.ForeignKey('Accession.id')))


class Reference(SqlAlchemyBase):
    """
    A Reference is a canonical assembled contig, intended to act as a
    reference coordinate space for other genomic annotations. A single
    Reference might represent the human chromosome 1, for instance.
    """

    __tablename__ = 'Reference'

    id = create_id_column()

    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    The name of this reference, eg. "22". This must be unique
    within a ReferenceSet.
    """

    md5checksum = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    The MD5 checksum uniquely representing this `Reference` as a
    lower-case hexadecimal string, calculated as the MD5 of the upper-case
    sequence excluding all whitespace characters.
    """

    length = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    """
    The length of this reference's sequence string.
    """

    ncbi_taxon_id = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    """
    The NCBI Taxon ID for this reference. This is the
    ID from http://www.ncbi.nlm.nih.gov/taxonomy (e.g. 9606->human)
    indicating the species which this reference is intended to model.
    Note that `Reference`s within a ReferenceSet may specify a different
    `ncbiTaxonId`, as assemblies may contain reference sequences
    which do not belong to the modeled species, e.g.  EBV in a
    human reference genome.
    """

    source_uri = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    The URI from which the sequence was obtained. Specifies a FASTA format
    file/string with one name, sequence pair.
    """

    is_derived = sqlalchemy.Column(sqlalchemy.Boolean, nullable=False)
    """
    True if this Reference is derived. A sequence X is said to be
    derived from source sequence Y, if X and Y are of the same length and
    the per-base sequence divergence at A/C/G/T bases is sufficiently
    small. Two sequences derived from the same official sequence share the
    same coordinates and annotations, and can be replaced with the official
    sequence for certain use cases.
    """

    source_divergence = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    """
    The source divergence reference of a reference is the fraction of non-indel
    bases that do not match the reference this record was derived from.
    """

    source_accessions = orm.relationship(
        "Accession", secondary=_reference_accessions_table)
    """
    The list of source accessions for this Reference. These are all known
    corresponding accession IDs in INSDC (GenBank/ENA/DDBJ) ideally
    with a version number, e.g. `NC_000001.11`.
    """

    reference_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReferenceSet.id"),
        nullable=False)
    reference_set = orm.relationship(
        "ReferenceSet", back_populates="references")

    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity': 'Reference',
         'polymorphic_on': type
    }
    __table_args__ = (
        # Reference names must be unique within a reference set.
        sqlalchemy.UniqueConstraint("reference_set_id", "name"),
    )

    def __init__(self, name):
        self.name = name
        self.source_uri = ""
        self.ncbi_taxon_id = 0
        self.is_derived = False
        self.source_divergence = 0

    def get_protobuf(self):
        ret = protocol.Reference()
        ret.id = str(self.id)
        ret.name = self.name
        ret.md5checksum = self.md5checksum
        ret.length = self.length
        ret.ncbi_taxon_id = self.ncbi_taxon_id
        ret.source_accessions.extend(
            acc.name for acc in self.source_accessions)
        ret.is_derived = self.is_derived
        ret.source_divergence = self.source_divergence
        ret.source_uri = self.source_uri
        return ret

    def check_query_range(self, start, end):
        """
        Checks to ensure that the query range is valid within this reference.
        If not, raise ReferenceRangeErrorException.
        """
        if start < 0 or end > self.length or start > end:
            raise exceptions.ReferenceRangeErrorException(self.id, start, end)


# There is a many-to-many association between ReferenceSets and Accessions,
# so we need an association table.
_reference_set_accessions_table = sqlalchemy.Table(
    'reference_set_accession_association', SqlAlchemyBase.metadata,
    sqlalchemy.Column(
        'reference_set_id', sqlalchemy.Integer,
        sqlalchemy.ForeignKey('ReferenceSet.id')),
    sqlalchemy.Column(
        'accession_id', sqlalchemy.Integer,
        sqlalchemy.ForeignKey('Accession.id')))


class ReferenceSet(SqlAlchemyBase):
    """
    A ReferenceSet is a set of References which typically comprise a
    reference assembly, such as GRCh38.
    """
    __tablename__ = 'ReferenceSet'

    id = create_id_column()

    name = sqlalchemy.Column(sqlalchemy.String, nullable=False, unique=True)
    """
    The name of this ReferenceSet. This must be unique within the repository.
    """

    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    Returns the free text description of this reference set.
    """

    assembly_id = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    This is the public id of this reference set, such as `GRCh37`
    """

    md5checksum = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    The MD5 checksum for this reference set. This checksum is
    calculated by making a list of `Reference.md5checksum` for all
    `Reference`s in this set. We then sort this list, and take the
    MD5 hash of all the strings concatenated together.
    """

    is_derived = sqlalchemy.Column(sqlalchemy.Boolean, nullable=False)
    """
    True if this ReferenceSet is derived. A ReferenceSet
    may be derived from a source if it contains additional sequences,
    or some of the sequences within it are derived.
    """

    ncbi_taxon_id = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    """
    The NCBI Taxon ID for this reference set. This is the
    ID from http://www.ncbi.nlm.nih.gov/taxonomy (e.g. 9606->human)
    indicating the species which this assembly is intended to model.
    Note that contained `Reference`s may specify a different
    `ncbiTaxonId`, as assemblies may contain reference sequences
    which do not belong to the modeled species, e.g.  EBV in a
    human reference genome.
    """

    source_uri = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    The sourceURI for this ReferenceSet.

    TODO: clarify the purpose and meaning of this field.
    """

    source_accessions = orm.relationship(
        "Accession", secondary=_reference_set_accessions_table)
    """
    The list of source accessions for this ReferenceSet. These are all known
    corresponding accession IDs in INSDC (GenBank/ENA/DDBJ) ideally
    with a version number, e.g. `NC_000001.11`.
    """

    references = orm.relationship(
        "Reference", back_populates="reference_set",
        cascade="all, delete, delete-orphan")
    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity': 'ReferenceSet',
         'polymorphic_on': type
    }

    def __init__(self, name):
        self.name = name
        self.assembly_id = ""
        self.description = ""
        self.is_derived = False
        self.ncbi_taxon_id = 0
        self.source_uri = ""

    def _compute_md5checksum(self):
        """
        Computes the MD5 checksum for this reference set. See the
        documentation for the md5checksum for how this is defined.
        """
        checksums = sorted([ref.md5checksum for ref in self.references])
        self.md5checksum = hashlib.md5("".join(checksums)).hexdigest()

    def get_protobuf(self):
        ret = protocol.ReferenceSet()
        ret.id = str(self.id)
        ret.name = self.name
        ret.assembly_id = self.assembly_id
        ret.description = self.description
        ret.is_derived = self.is_derived
        ret.md5checksum = self.md5checksum
        ret.ncbi_taxon_id = self.ncbi_taxon_id
        ret.source_accessions.extend(
            acc.name for acc in self.source_accessions)
        ret.source_uri = self.source_uri
        return ret

#####################################################################
#
# Variants
#
#####################################################################


class VariantSetMetadata(SqlAlchemyBase):
    __tablename__ = "VariantSetMetadata"

    id = create_id_column()
    key = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    value = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    type = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    number = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    variant_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("VariantSet.id"),
        nullable=False)
    variant_set = orm.relationship(
        "VariantSet", back_populates="variant_set_metadata")

    def __init__(self, key=None, value="", type="", number="", description=""):
        self.key = key
        self.value = value
        self.type = type
        self.number = number
        self.description = description

    def get_protobuf(self):
        ret = protocol.VariantSetMetadata()
        ret.id = str(self.id)
        ret.key = self.key
        ret.value = self.value
        ret.type = self.type
        ret.number = self.number
        ret.description = self.description
        return ret


# There is a many-to-many association between VariantSets and CallSets,
# so we need an association table.
_call_set_association_table = sqlalchemy.Table(
    'call_set_association', SqlAlchemyBase.metadata,
    sqlalchemy.Column(
        'variant_set_id', sqlalchemy.Integer,
        sqlalchemy.ForeignKey('VariantSet.id')),
    sqlalchemy.Column(
        'call_set_id', sqlalchemy.Integer,
        sqlalchemy.ForeignKey('CallSet.id')))


class CallSet(SqlAlchemyBase):
    __tablename__ = 'CallSet'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()

    bio_sample_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("BioSample.id"),
        nullable=True)
    bio_sample = orm.relationship("BioSample")
    variant_sets = orm.relationship(
        "VariantSet", secondary=_call_set_association_table,
        back_populates="call_sets")

    def get_protobuf(self):
        ret = protocol.CallSet()
        ret.id = str(self.id)
        ret.name = self.name
        if self.bio_sample_id is not None:
            ret.bio_sample_id = str(self.bio_sample_id)
        ret.created = protocol.datetime_to_milliseconds(
            self.creation_timestamp)
        ret.updated = protocol.datetime_to_milliseconds(
            self.update_timestamp)
        ret.variant_set_ids.extend(str(vs.id) for vs in self.variant_sets)
        # TODO create a generic Info table that can hold these values
        # for a variety of different tables. Or it may be simpler to create
        # a per-class table.
        # for key in self._info:
        #     gaCallSet.info[key].values.extend(_encodeValue(self._info[key]))
        return ret

    def __init__(self, name=None):
        self.name = name


class VariantSet(SqlAlchemyBase):
    __tablename__ = 'VariantSet'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()

    dataset_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Dataset.id"),
        nullable=False)
    dataset = orm.relationship(
        "Dataset", back_populates="variant_sets", single_parent=True)
    reference_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReferenceSet.id"),
        nullable=False)
    reference_set = orm.relationship("ReferenceSet", single_parent=True)

    variant_set_metadata = orm.relationship(
        "VariantSetMetadata", cascade="all, delete, delete-orphan")
    call_sets = orm.relationship(
        "CallSet", secondary=_call_set_association_table,
        back_populates="variant_sets", lazy="dynamic")
    variant_annotation_sets = orm.relationship(
        "VariantAnnotationSet", back_populates="variant_set",
        cascade="all, delete, delete-orphan")

    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity': 'VariantSet',
         'polymorphic_on': type
    }
    __table_args__ = (
        # VariantSet names must be unique within a dataset
        sqlalchemy.UniqueConstraint("dataset_id", "name"),
    )

    def __init__(self, name):
        self.name = name

    def get_protobuf(self):
        ret = protocol.VariantSet()
        ret.id = str(self.id)
        ret.name = self.name
        ret.dataset_id = str(self.dataset_id)
        ret.reference_set_id = str(self.reference_set_id)
        metadata = [m.get_protobuf() for m in self.variant_set_metadata]
        ret.metadata.extend(metadata)
        return ret

    def is_annotated(self):
        """
        Returns True if there are any variant annotation sets associated with
        this VariantSet.
        """
        return len(self.variant_annotation_sets) > 0

    def get_variants(self, reference_name, start=0, end=None):
        """
        Convenience method for testing. Returns an iterator over the
        variants in the specified region. If end is not specified
        all variants from start are returned.
        """
        raise NotImplementedError()


class Analysis(SqlAlchemyBase):

    __tablename__ = "Analysis"

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    # TODO This should really be a program like in the ReadGroup
    # so I'm not going to bother making a proper list.
    software = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()

    def __init__(
            self, name=None, description=None, software=None,
            creation_timestamp=None, update_timestamp=None):
        self.name = name
        self.description = description
        self.software = software
        # TODO add support for INFO
        if creation_timestamp is not None:
            self.creation_timestamp = creation_timestamp
        if update_timestamp is not None:
            self.update_timestamp = update_timestamp

    def get_protobuf(self):
        ret = protocol.Analysis()
        ret.id = str(self.id)
        ret.name = self.name
        ret.description = self.description
        # We're making software a single string for now. This should be
        # updated to use Program when we harmonise the data model
        ret.software.append(self.software)
        ret.created = protocol.datetime_to_iso8601(self.creation_timestamp)
        ret.updated = protocol.datetime_to_iso8601(self.update_timestamp)
        return ret


class VariantAnnotationSet(SqlAlchemyBase):

    __tablename__ = "VariantAnnotationSet"

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()
    variant_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("VariantSet.id"),
        nullable=False)
    variant_set = orm.relationship(
        "VariantSet", back_populates="variant_annotation_sets",
        single_parent=True)
    analysis_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Analysis.id"),
        nullable=False)
    analysis = orm.relationship("Analysis", single_parent=True)

    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity': 'VariantAnnotationSet',
         'polymorphic_on': type
    }
    __table_args__ = (
        # Annotation set names must be unique within a variant set
        sqlalchemy.UniqueConstraint("variant_set_id", "name"),
    )

    def __init__(self, name):
        self.name = name

    def get_protobuf(self):
        ret = protocol.VariantAnnotationSet()
        ret.id = str(self.id)
        ret.name = self.name
        # Casing needs to be fixed at the schema level
        ret.variantSet_id = str(self.variant_set_id)
        ret.analysis.CopyFrom(self.analysis.get_protobuf())
        return ret

    def getTranscriptEffectId(self, gaTranscriptEffect):
        effs = [eff.term for eff in gaTranscriptEffect.effects]
        return hashlib.md5(
            "{}\t{}\t{}\t{}".format(
                gaTranscriptEffect.alternate_bases,
                gaTranscriptEffect.feature_id,
                effs, gaTranscriptEffect.hgvs_annotation)
            ).hexdigest()

    def hashVariantAnnotation(cls, gaVariant, gaVariantAnnotation):
        """
        Produces an MD5 hash of the gaVariant and gaVariantAnnotation objects
        """
        treffs = [treff.id for treff in gaVariantAnnotation.transcript_effects]
        return hashlib.md5(
            "{}\t{}\t{}\t".format(
                gaVariant.reference_bases, tuple(gaVariant.alternate_bases),
                treffs)
            ).hexdigest()

    def accept_annotation(self, annotation, request):
        """
        Returns True if the specified annotation should be included in a
        response to the specified request.
        """
        # TODO This is not thoroughly tested and is probably naive!
        effect_ids = set(effect.id for effect in request.effects)
        if len(effect_ids) == 0:
            # If no effect IDs are provided, we accept everything
            return True
        for treff in annotation.transcript_effects:
            for effect in treff.effects:
                if effect.id in effect_ids:
                    return True
        return False


#####################################################################
#
# Reads
#
#####################################################################


class ReadStats(SqlAlchemyBase):
    __tablename__ = 'ReadStats'

    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    aligned_read_count = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    unaligned_read_count = sqlalchemy.Column(
        sqlalchemy.Integer, nullable=False)
    base_count = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)

    def __init__(
            self, aligned_read_count=0, unaligned_read_count=0, base_count=0):
        self.aligned_read_count = aligned_read_count
        self.unaligned_read_count = unaligned_read_count
        self.base_count = base_count

    def update_protobuf(self, stats):
        stats.aligned_read_count = self.aligned_read_count
        stats.unaligned_read_count = self.unaligned_read_count
        stats.base_count = self.base_count


class Program(SqlAlchemyBase):
    __tablename__ = 'Program'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    version = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    command_line = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    prev_program_id = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    index = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    """
    The position of this program within a given chain. Used to sort
    programs to ensure they are returned in the most useful order.
    """
    read_group_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReadGroupSet.id"),
        nullable=False)

    def __init__(
            self, name="", version="", command_line="", prev_program_id="",
            index=None):
        self.name = name
        self.index = index
        self.version = version
        self.command_line = command_line
        self.prev_program_id = prev_program_id

    def get_protobuf(self):
        ret = protocol.ReadGroup.Program()
        ret.id = str(self.id)
        ret.name = self.name
        ret.version = self.version
        ret.command_line = self.command_line
        ret.prev_program_id = self.prev_program_id
        return ret


class Experiment(SqlAlchemyBase):
    __tablename__ = 'Experiment'

    id = create_id_column()
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()
    instrument_model = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    sequencing_center = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    library = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    platform_unit = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    run_time = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    def __init__(
            self, instrument_model="", sequencing_center="", description="",
            library="", platform_unit="", run_time=""):
        self.instrument_model = instrument_model
        self.sequencing_center = sequencing_center
        self.description = description
        self.library = library
        self.platform_unit = platform_unit
        self.run_time = run_time
        self.creation_timestamp = datetime.datetime.now()
        self.update_timestamp = datetime.datetime.now()

    def update_protobuf(self, experiment):
        experiment.id = str(self.id)
        experiment.instrument_model = self.instrument_model
        experiment.sequencing_center = self.sequencing_center
        experiment.description = self.description
        experiment.library = self.library
        experiment.platform_unit = self.platform_unit
        experiment.run_time = self.run_time
        experiment.message_create_time = protocol.datetime_to_iso8601(
            self.creation_timestamp)
        experiment.message_update_time = protocol.datetime_to_iso8601(
            self.update_timestamp)


class ReadGroup(SqlAlchemyBase):
    """
    Class representing a ReadGroup. A ReadGroup is all the data that's
    processed the same way by the sequencer.  There are typically 1-10
    ReadGroups in a ReadGroupSet.
    """
    __tablename__ = 'ReadGroup'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()
    predicted_insert_size = sqlalchemy.Column(
        sqlalchemy.Integer, nullable=False)
    sample_name = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    bio_sample_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("BioSample.id"),
        nullable=True)
    bio_sample = orm.relationship("BioSample")

    read_group_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReadGroupSet.id"),
        nullable=False)
    read_group_set = orm.relationship(
        "ReadGroupSet", back_populates="read_groups")
    read_stats_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReadStats.id"),
        nullable=False)
    read_stats = orm.relationship(
        "ReadStats",  cascade="all, delete, delete-orphan", single_parent=True)
    experiment_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Experiment.id"),
        nullable=True)
    experiment = orm.relationship(
        "Experiment",  cascade="all, delete, delete-orphan",
        single_parent=True)

    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity': 'ReadGroup',
         'polymorphic_on': type
    }
    __table_args__ = (
        # ReadGroup names must be unique within a read group set.
        sqlalchemy.UniqueConstraint("read_group_set_id", "name"),
    )

    def __init__(self, name=None, sample_name=None):
        self.name = name
        self.sample_name = sample_name
        self.description = ""
        self.predicted_insert_size = 0
        self.read_stats = ReadStats()
        self.creation_timestamp = datetime.datetime.now()
        self.update_timestamp = datetime.datetime.now()

    def get_protobuf(self):
        ret = protocol.ReadGroup()
        ret.id = str(self.id)
        ret.name = self.name
        ret.description = self.description
        if self.sample_name is not None:
            ret.sample_name = self.sample_name
        if self.bio_sample_id is not None:
            ret.bio_sample_id = str(self.bio_sample_id)
        ret.created = protocol.datetime_to_milliseconds(
            self.creation_timestamp)
        ret.updated = protocol.datetime_to_milliseconds(
            self.update_timestamp)
        ret.dataset_id = str(self.read_group_set.dataset_id)
        ret.reference_set_id = str(self.read_group_set.reference_set_id)
        ret.predicted_insert_size = self.predicted_insert_size
        self.read_stats.update_protobuf(ret.stats)
        ret.programs.extend(
            p.get_protobuf() for p in self.read_group_set.programs)
        if self.experiment is not None:
            self.experiment.update_protobuf(ret.experiment)
        return ret


class ReadGroupSet(SqlAlchemyBase):

    __tablename__ = 'ReadGroupSet'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    dataset_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Dataset.id"),
        nullable=False)
    dataset = orm.relationship(
        "Dataset", back_populates="read_group_sets", single_parent=True)
    reference_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReferenceSet.id"),
        nullable=False)
    reference_set = orm.relationship("ReferenceSet", single_parent=True)
    read_stats_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReadStats.id"),
        nullable=False)
    read_stats = orm.relationship(
        "ReadStats",  cascade="all, delete, delete-orphan", single_parent=True)
    # For simplicity we associate the programs with the ReadGroupSet rather
    # than the read group, as the protocol does.
    programs = orm.relationship(
        "Program",  cascade="all, delete, delete-orphan",
        single_parent=True, order_by="Program.index")
    read_groups = orm.relationship(
        "ReadGroup", back_populates="read_group_set",
        cascade="all, delete, delete-orphan")

    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity': 'ReadGroupSet',
         'polymorphic_on': type
    }
    __table_args__ = (
        # ReadGroupSet names must be unique within a dataset
        sqlalchemy.UniqueConstraint("dataset_id", "name"),
    )

    def __init__(self, name):
        self.name = name
        self.read_stats = ReadStats()
        self.experiment = Experiment()

    def get_protobuf(self):
        ret = protocol.ReadGroupSet()
        ret.id = str(self.id)
        ret.name = self.name
        ret.dataset_id = str(self.dataset_id)
        ret.read_groups.extend(
            rg.get_protobuf() for rg in self.read_groups)
        self.read_stats.update_protobuf(ret.stats)
        return ret


#####################################################################
#
# Features
#
#####################################################################


class FeatureSet(SqlAlchemyBase):
    """
    A set of sequence features annotations
    """
    __tablename__ = 'FeatureSet'

    id = create_id_column()
    """
    The FeatureSet id. Unique within the repository.
    """

    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    """
    The name of this FeatureSet. Unique within a dataset.
    """

    source_uri = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    dataset_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Dataset.id"),
        nullable=False)
    dataset = orm.relationship(
        "Dataset", back_populates="feature_sets", single_parent=True)
    reference_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReferenceSet.id"),
        nullable=False)
    reference_set = orm.relationship("ReferenceSet", single_parent=True)

    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity': 'FeatureSet',
         'polymorphic_on': type
    }
    __table_args__ = (
        # FeatureSet names must be unique within a dataset
        sqlalchemy.UniqueConstraint("dataset_id", "name"),
    )

    def __init__(self, name, source_uri=""):
        self.name = name
        self.source_uri = source_uri

    def get_protobuf(self):
        ret = protocol.FeatureSet()
        ret.id = str(self.id)
        ret.name = self.name
        ret.dataset_id = str(self.dataset_id)
        ret.source_uri = self.source_uri
        # TODO Add info values
        return ret


#####################################################################
#
# The registry DB.
#
#####################################################################

class RegistryDb(object):
    """
    The database representing all of the top-level objects in a given
    data repository.
    """
    def __init__(self, db_url):
        self._db_url = db_url
        self._engine = None
        self._session = None

    def get_new_id(self):
        # TODO use an instance of SystemRandom to do this and check against the
        # global IDs table. We should not have the possibility of having two
        # objects with the same ID. We could avoid the problem by using the
        # first 4 bits to identify the type.
        return random.randint(0, 2**31)

    # Miscellaneous house-keeping methods for the Registry DB.

    def get_db_url(self):
        """
        Returns the SqlAlchemy URL for the underlying database.
        """
        return self._db_url

    def exists(self):
        """
        Returns True if there is a fully initialised registry in the DB url.
        """
        return True

    def open(self):
        """
        Opens this registry DB.
        """
        self._engine = sqlalchemy.create_engine(self._db_url, echo=False)
        Session = orm.sessionmaker()
        Session.configure(bind=self._engine)
        self._session = Session()

    def close(self):
        """
        Closes this repo.
        """
        self._session.close()
        self._session = None
        self._engine = None

    def delete(self):
        """
        Drops all tables in the DB.
        """
        self.open()
        metadata = SqlAlchemyBase.metadata
        metadata.drop_all(self._engine)
        self.close()

    def initialise(self):
        """
        Initialise this data repostitory by creating the registry DB.
        """
        metadata = SqlAlchemyBase.metadata
        metadata.drop_all(self._engine)
        metadata.create_all(self._engine)

    def verify(self):
        print("NOT IMPLEMENTED YET!!!")

    def print_summary(self):
        """
        Prints out a summary of the registry.
        """
        print("Datasets = ", self._session.query(Dataset).count())

    # Top level API to add objects to the data repository.

    def remove(self, obj):
        """
        Removes the specified object from the database.
        """
        self._session.delete(obj)
        self._session.commit()

    def _is_unique_name_violation(self, integrity_error):
        """
        Returns True if the specified exception corresponds to a
        unique name violation.
        """
        # This seems to work for Postgres and Sqlite...
        return "name" in str(integrity_error.orig)

    def __add(self, db_object):
        try:
            self._session.add(db_object)
            self._session.commit()
        except sqlalchemy.exc.IntegrityError as ie:
            if self._is_unique_name_violation(ie):
                raise exceptions.DuplicateNameException(db_object.name)
            else:
                raise ie

    def add(self, db_object):
        # TODO remove the other specific methods, and updated the signature
        # so that we can add a bunch of objects in a single transaction.
        self.__add(db_object)

    def add_ontology(self, ontology):
        """
        Adds the specified ontology to this data repository.
        """
        self.__add(ontology)

    def add_reference_set(self, reference_set):
        """
        Adds the specified reference set to this data repository.
        """
        self.__add(reference_set)

    def add_dataset(self, dataset):
        """
        Adds the specified dataset to this data repository.
        """
        self.__add(dataset)

    def add_variant_set(self, variant_set):
        """
        Adds the specified variant_set to this data repository.
        """
        self.__add(variant_set)

    def add_read_group_set(self, read_group_set):
        """
        Adds the specified variant_set to this data repository.
        """
        self.__add(read_group_set)

    def add_feature_set(self, feature_set):
        """
        Adds the specified feature_set to this data repository.
        """
        self.__add(feature_set)

    # Object accessors by name.

    def get_dataset_by_name(self, name):
        """
        Returns the dataset with the specified name.
        """
        result = self._session.query(Dataset).filter(
            Dataset.name == name).first()
        if result is None:
            raise exceptions.DatasetNameNotFoundException(name)
        return result

    def get_ontology_by_name(self, name):
        """
        Returns the ontology with the specified name.
        """
        result = self._session.query(Ontology).filter(
            Ontology.name == name).first()
        if result is None:
            raise exceptions.OntologyNameNotFoundException(name)
        return result

    def get_reference_set_by_name(self, name):
        """
        Returns the reference_set with the specified name.
        """
        result = self._session.query(ReferenceSet).filter(
            ReferenceSet.name == name).first()
        if result is None:
            raise exceptions.ReferenceSetNameNotFoundException(name)
        return result

    def get_variant_set_by_name(self, dataset_name, name):
        """
        Returns the variant_set with the specified name from the dataset
        with the specified name.
        """
        query = self._session.query(VariantSet).join(VariantSet.dataset)
        query = query.filter(
            VariantSet.name == name, Dataset.name == dataset_name)
        result = query.first()
        if result is None:
            raise exceptions.VariantSetNameNotFoundException(name)
        return result

    def get_feature_set_by_name(self, dataset_name, name):
        """
        Returns the feature_set with the specified name from the dataset
        with the specified name.
        """
        query = self._session.query(FeatureSet).join(FeatureSet.dataset)
        query = query.filter(
            FeatureSet.name == name, Dataset.name == dataset_name)
        result = query.first()
        if result is None:
            raise exceptions.FeatureSetNameNotFoundException(name)
        return result

    def get_read_group_set_by_name(self, dataset_name, name):
        """
        Returns the read_group_set with the specified name from the dataset
        with the specified name.
        """
        query = self._session.query(ReadGroupSet).join(ReadGroupSet.dataset)
        query = query.filter(
            ReadGroupSet.name == name, Dataset.name == dataset_name)
        result = query.first()
        if result is None:
            raise exceptions.ReadGroupSetNameNotFoundException(name)
        return result

    # Object accessors by ID. These are run when the corresponding GET request
    # is received.

    def get_dataset(self, id_):
        """
        Retuns the Dataset with the specified ID, or raises a
        DatasetNotFoundException if it does not exist.
        """
        result = self._session.query(Dataset).filter(
            Dataset.id == id_).first()
        if result is None:
            raise exceptions.DatasetNotFoundException(id_)
        return result

    def get_read_group(self, id_):
        """
        Retuns the ReadGroup with the specified ID, or raises a
        ReadGroupNotFoundException if it does not exist.
        """
        result = self._session.query(ReadGroup).filter(
            ReadGroup.id == id_).first()
        if result is None:
            raise exceptions.ReadGroupNotFoundException(id_)
        return result

    def get_read_group_set(self, id_):
        """
        Retuns the ReadGroupSet with the specified ID, or raises a
        ReadGroupSetNotFoundException if it does not exist.
        """
        result = self._session.query(ReadGroupSet).filter(
            ReadGroupSet.id == id_).first()
        if result is None:
            raise exceptions.ReadGroupSetNotFoundException(id_)
        return result

    def get_call_set(self, id_):
        """
        Retuns the CallSet with the specified ID, or raises a
        CallSetNotFoundException if it does not exist.
        """
        result = self._session.query(CallSet).filter(CallSet.id == id_).first()
        if result is None:
            raise exceptions.CallSetNotFoundException(id_)
        return result

    def get_variant_set(self, id_):
        """
        Retuns the VariantSet with the specified ID, or raises a
        VariantSetNotFoundException if it does not exist.
        """
        result = self._session.query(VariantSet).filter(
            VariantSet.id == id_).first()
        if result is None:
            raise exceptions.VariantSetNotFoundException(id_)
        return result

    def get_variant_annotation_set(self, id_):
        """
        Retuns the VariantAnnotationSet with the specified ID, or raises a
        VariantAnnotationSetNotFoundException if it does not exist.
        """
        result = self._session.query(VariantAnnotationSet).filter(
            VariantAnnotationSet.id == id_).first()
        if result is None:
            raise exceptions.VariantAnnotationSetNotFoundException(id_)
        return result

    def get_feature_set(self, id_):
        """
        Retuns the FeatureSet with the specified ID, or raises a
        FeatureSetNotFoundException if it does not exist.
        """
        result = self._session.query(FeatureSet).filter(
            FeatureSet.id == id_).first()
        if result is None:
            raise exceptions.FeatureSetNotFoundException(id_)
        return result

    def get_reference_set(self, id_):
        """
        Retuns the ReferenceSet with the specified ID, or raises a
        ReferenceSetNotFoundException if it does not exist.
        """
        result = self._session.query(ReferenceSet).filter(
            ReferenceSet.id == id_).first()
        if result is None:
            raise exceptions.ReferenceSetNotFoundException(id_)
        return result

    def get_reference(self, id_):
        """
        Returns the Reference with the specified ID or raises a
        ReferenceNotFoundException if it does not exist.
        """
        result = self._session.query(Reference).filter(
            Reference.id == id_).first()
        if result is None:
            raise exceptions.ReferenceNotFoundException(id_)
        return result

    def get_bio_sample(self, id_):
        """
        Returns the BioSample with the specified ID or raises a
        BioSampleNotFoundException if it does not exist.
        """
        result = self._session.query(BioSample).filter(
            BioSample.id == id_).first()
        if result is None:
            raise exceptions.BioSampleNotFoundException(id_)
        return result

    def get_individual(self, id_):
        """
        Returns the Individual with the specified ID or raises a
        IndividualNotFoundException if it does not exist.
        """
        result = self._session.query(Individual).filter(
            Individual.id == id_).first()
        if result is None:
            raise exceptions.IndividualNotFoundException(id_)
        return result

    # General getters to iterate over all objects in a class.

    def get_reference_sets(self):
        return self._session.query(ReferenceSet).all()

    def get_datasets(self):
        return self._session.query(Dataset).all()

    def get_variant_sets(self):
        return self._session.query(VariantSet).all()

    def get_variant_annotation_sets(self):
        return self._session.query(VariantAnnotationSet).all()

    # Getters to provide queries corresponding to the external search requests.

    def get_references_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchReferencesRequest.
        """
        query = self._session.query(Reference).filter(
            Reference.reference_set_id == request.reference_set_id)
        if request.accession:
            query = query.filter(
                Reference.source_accessions.any(name=request.accession))
        if request.md5checksum:
            query = query.filter(Reference.md5checksum == request.md5checksum)
        return query

    def get_reference_sets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchReferenceSetsRequest.
        """
        query = self._session.query(ReferenceSet)
        if request.accession:
            query = query.filter(
                ReferenceSet.source_accessions.any(name=request.accession))
        if request.md5checksum:
            query = query.filter(
                ReferenceSet.md5checksum == request.md5checksum)
        if request.assembly_id:
            query = query.filter(
                ReferenceSet.assembly_id == request.assembly_id)
        return query

    def get_datasets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchDatasetsRequest.
        """
        query = self._session.query(Dataset)
        return query

    def get_variant_sets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchVariantSets request.
        """
        query = self._session.query(VariantSet).filter(
            VariantSet.dataset_id == request.dataset_id)
        return query

    def get_variant_annotation_sets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchVariantAnnotationSets request.
        """
        query = self._session.query(VariantAnnotationSet).filter(
            VariantAnnotationSet.variant_set_id == request.variant_set_id)
        return query

    def get_feature_sets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchFeatureSets request.
        """
        query = self._session.query(FeatureSet).filter(
            FeatureSet.dataset_id == request.dataset_id)
        return query

    def get_call_sets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchCallSets request.
        """
        query = self._session.query(CallSet)
        query = query.filter(
            CallSet.variant_sets.any(id=request.variant_set_id))
        if request.name:
            query = query.filter(CallSet.name == request.name)
        if request.bio_sample_id:
            query = query.filter(
                CallSet.bio_sample_id == request.bio_sample_id)
        return query

    def get_read_group_sets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchReadGroupSetsRequest.
        """
        query = self._session.query(ReadGroupSet).filter(
            ReadGroupSet.dataset_id == request.dataset_id)
        if request.name:
            query = query.filter(ReadGroupSet.name == request.name)
        if request.bio_sample_id:
            raise exceptions.NotImplementedException(
                "search by bio sample ID not implemented")
        return query

    def get_individuals_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchIndividualsRequest.
        """
        query = self._session.query(Individual).filter(
            Individual.dataset_id == request.dataset_id)
        if request.name:
            query = query.filter(Individual.name == request.name)
        return query

    def get_bio_samples_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchBioSamplesRequest.
        """
        query = self._session.query(BioSample).filter(
            BioSample.dataset_id == request.dataset_id)
        if request.name:
            query = query.filter(BioSample.name == request.name)
        if request.individual_id:
            query = query.filter(
                BioSample.individual_id == request.individual_id)
        return query


#################
# TODO Compound ID class should be greatly simplified.
#################

class CompoundId(object):
    """
    Base class for an id composed of several different parts.  Each
    compound ID consists of a set of fields, each of which corresponds to a
    local ID in the data hierarchy. For example, we might have fields like
    ["dataset", "variantSet"] for a variantSet.  These are available as
    cid.dataset, and cid.variantSet.  The actual IDs of the containing
    objects can be obtained using the corresponding attributes, e.g.
    cid.datasetId and cid.variantSetId.
    """
    fields = []
    """
    The fields that the compound ID is composed of. These are parsed and
    made available as attributes on the object.
    """
    containerIds = []
    """
    The fields of the ID form a breadcrumb trail through the data
    hierarchy, and successive prefixes provide the IDs for objects
    further up the tree. This list is a set of tuples giving the
    name and length of a given prefix forming an identifier.
    """
    differentiator = None
    """
    A string used to guarantee unique ids for objects.  A value of None
    indicates no string is used.  Otherwise, this string will be spliced
    into the object's id.
    """
    differentiatorFieldName = 'differentiator'
    """
    The name of the differentiator field in the fields array for CompoundId
    subclasses.
    """

    def __init__(self, parentCompoundId, *localIds):
        """
        Allocates a new CompoundId for the specified parentCompoundId and
        local identifiers. This compoundId inherits all of the fields and
        values from the parent compound ID, and must have localIds
        corresponding to its fields. If no parent id is present,
        parentCompoundId should be set to None.
        """
        index = 0
        if parentCompoundId is not None:
            for field in parentCompoundId.fields:
                setattr(self, field, getattr(parentCompoundId, field))
                index += 1
        if (self.differentiator is not None and
                self.differentiatorFieldName in self.fields[index:]):
            # insert a differentiator into the localIds if appropriate
            # for this class and we haven't advanced beyond it already
            differentiatorIndex = self.fields[index:].index(
                self.differentiatorFieldName)
            localIds = localIds[:differentiatorIndex] + tuple([
                self.differentiator]) + localIds[differentiatorIndex:]
        for field, localId in zip(self.fields[index:], localIds):
            if not isinstance(localId, basestring):
                raise exceptions.BadIdentifierNotStringException(localId)
            encodedLocalId = self.encode(localId)
            setattr(self, field, encodedLocalId)
        if len(localIds) != len(self.fields) - index:
            raise ValueError(
                "Incorrect number of fields provided to instantiate ID")
        for idFieldName, prefix in self.containerIds:
            values = [getattr(self, f) for f in self.fields[:prefix + 1]]
            containerId = self.join(values)
            obfuscated = self.obfuscate(containerId)
            setattr(self, idFieldName, obfuscated)

    def __str__(self):
        values = [getattr(self, f) for f in self.fields]
        compoundIdStr = self.join(values)
        return self.obfuscate(compoundIdStr)

    @classmethod
    def join(cls, splits):
        """
        Join an array of ids into a compound id string
        """
        segments = []
        for split in splits:
            segments.append('"{}",'.format(split))
        if len(segments) > 0:
            segments[-1] = segments[-1][:-1]
        jsonString = '[{}]'.format(''.join(segments))
        return jsonString

    @classmethod
    def split(cls, jsonString):
        """
        Split a compound id string into an array of ids
        """
        splits = json.loads(jsonString)
        return splits

    @classmethod
    def encode(cls, idString):
        """
        Encode a string by escaping problematic characters
        """
        return idString.replace('"', '\\"')

    @classmethod
    def decode(cls, encodedString):
        """
        Decode an encoded string
        """
        return encodedString.replace('\\"', '"')

    @classmethod
    def parse(cls, compoundIdStr):
        """
        Parses the specified compoundId string and returns an instance
        of this CompoundId class.

        :raises: An ObjectWithIdNotFoundException if parsing fails. This is
        because this method is a client-facing method, and if a malformed
        identifier (under our internal rules) is provided, the response should
        be that the identifier does not exist.
        """
        if not isinstance(compoundIdStr, basestring):
            raise exceptions.BadIdentifierException(compoundIdStr)
        try:
            deobfuscated = cls.deobfuscate(compoundIdStr)
        except TypeError:
            # When a string that cannot be converted to base64 is passed
            # as an argument, b64decode raises a TypeError. We must treat
            # this as an ID not found error.
            raise exceptions.ObjectWithIdNotFoundException(compoundIdStr)
        try:
            encodedSplits = cls.split(deobfuscated)
            splits = [cls.decode(split) for split in encodedSplits]
        except (UnicodeDecodeError, ValueError):
            # Sometimes base64 decoding succeeds but we're left with
            # unicode gibberish. This is also and IdNotFound.
            raise exceptions.ObjectWithIdNotFoundException(compoundIdStr)
        # pull the differentiator out of the splits before instantiating
        # the class, if the differentiator exists
        fieldsLength = len(cls.fields)
        if cls.differentiator is not None:
            differentiatorIndex = cls.fields.index(
                cls.differentiatorFieldName)
            if differentiatorIndex < len(splits):
                del splits[differentiatorIndex]
            else:
                raise exceptions.ObjectWithIdNotFoundException(
                    compoundIdStr)
            fieldsLength -= 1
        if len(splits) != fieldsLength:
            raise exceptions.ObjectWithIdNotFoundException(compoundIdStr)
        return cls(None, *splits)

    @classmethod
    def obfuscate(cls, idStr):
        """
        Mildly obfuscates the specified ID string in an easily reversible
        fashion. This is not intended for security purposes, but rather to
        dissuade users from depending on our internal ID structures.
        """
        return unicode(base64.urlsafe_b64encode(
            idStr.encode('utf-8')).replace(b'=', b''))

    @classmethod
    def deobfuscate(cls, data):
        """
        Reverses the obfuscation done by the :meth:`obfuscate` method.
        If an identifier arrives without correct base64 padding this
        function will append it to the end.
        """
        # the str() call is necessary to convert the unicode string
        # to an ascii string since the urlsafe_b64decode method
        # sometimes chokes on unicode strings
        return base64.urlsafe_b64decode(str((
            data + b'A=='[(len(data) - 1) % 4:])))

    @classmethod
    def getInvalidIdString(cls):
        """
        Return an id string that is well-formed but probably does not
        correspond to any existing object; used mostly in testing
        """
        return cls.join(['notValid'] * len(cls.fields))


# TODO these are specific to the htslib implementation. We need an abstract
# class that just encodes the variant set id, and sub classes can add in
# what they want after that.
class VariantCompoundId(CompoundId):
    """
    The compound id for a variant
    """
    fields = ['variant_set_id', 'reference_name', 'start', 'md5']


class ReadAlignmentCompoundId(CompoundId):
    """
    The compound id for a read.
    """
    fields = ['read_group_id', 'query_name']
