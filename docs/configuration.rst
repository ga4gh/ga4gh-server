.. _configuration:

*************
Configuration
*************

The GA4GH reference server has two basic elements to its configuration:
the `Data repository`_ and the `Configuration file`_.

---------------
Data repository
---------------

Each GA4GH server represents a repository of information. This repository
consists of the reference sets, datasets, readgroup sets, variant sets etc. in
the server's data model and may contain data from many different unrelated
projects. The server administrator defines and manages this repository using
the ``ga4gh_repo`` command line interface, which provides commands to manage
all of the objects represented by a GA4GH server.

The information about the objects in the GA4GH data model is stored in an SQL
database, called the "registry DB" or "registry". The registry DB does not
contain the raw bulk data but rather "registers" the information about where
the server can find this information and some metadata about the object in
question. For example, if we have a variant set that is backed by a single VCF
file, the registry DB will contain the path to this file as well as the name of
the variant set, the reference set it is defined by, and other information
necessary to implement the GA4GH protocol. This registry architecture allows us
a lot of flexibility in the sources of data that we can use.

.. todo:: Architecture diagram showing the registry DB, the server, the
    data container objects (like HtslibVariantSet, etc), etc.

.. todo:: Explain how the datamodel is structured and how names are
    used as local identifiers.

+++++++
init
+++++++

The ``init`` command initialises a new registry DB at a given
file path. This is the first command that must be issued
when creating a new GA4GH repository.

.. argparse::
    :module: ga4gh.cli.repomanager
    :func: getRepoManagerParser
    :prog: ga4gh_repo
    :path: init
    :nodefault:


**Examples:**

.. code-block:: bash

    $ ga4gh_repo init registry.db

+++++++
list
+++++++

The ``list`` command is used to print the contents of a repository
to the screen. It is an essential tool for administrators to
understand the structure of the repository that they are managing.

.. note:: The ``list`` command is under development and will
   be much more sophisticated in the future. In particular, the output
   of this command should improve considerably in the near future.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: list
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo list registry.db


+++++++
verify
+++++++

The ``verify`` command is used to check that the integrity of the
data in a repository. The command checks each container object in turn
and ensures that it can read data from it. Read errors can occur for
any number of reasons (for example, a VCF file may have been moved
to another location since it was added to the registry), and the
``verify`` command allows an administrator to check that all is
well in their repository.

.. note:: The ``verify`` command is under development and will
   be much more sophisticated in the future. In particular, the output
   of this command should improve considerably in the near future.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: verify
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo verify registry.db


+++++++++++
add-dataset
+++++++++++

Creates a new dataset in a repository. A dataset is an arbitrary collection
of ReadGroupSets, VariantSets, VariantAnnotationSets and FeatureSets. Each
dataset has a name, which is used to identify it in the repository manager.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-dataset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-dataset registry.db 1kg -d 'Example dataset using 1000 genomes data'

Adds the dataset with the name ``1kg`` and description
``'Example dataset using 1000 genomes data'`` to the
registry database ``registry.db``.

++++++++++++++++
add-referenceset
++++++++++++++++

Adds a reference set derived from a FASTA file to a repository. Each
record in the FASTA file will correspond to a Reference in the new
ReferenceSet. The input FASTA file must be compressed with ``bgzip``
and indexed using ``samtools faidx``. Each ReferenceSet contains a
number of metadata values (.e.g. ``ncbiTaxonId``) which can be set
using command line options.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-referenceset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-referenceset registry.db hs37d5.fa.gz \
        --description "NCBI37 assembly of the human genome" \
        --ncbiTaxonId 9606 --name NCBI37 \
        --sourceUri ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

Adds a reference set used in the 1000 Genomes project using the name
``NCBI37``, also setting the ``ncbiTaxonId`` to 9606 (human).

++++++++++++++++
add-biosample
++++++++++++++++

Adds a new biosample to the repository. The biosample argument is
a JSON document according to the GA4GH JSON schema.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-biosample
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-biosample registry.db dataset1 HG00096 '{"individualId": "abc"}'

Adds the biosample named HG00096 to the repository with the individual ID
"abc".

++++++++++++++++
add-individual
++++++++++++++++

Adds a new individual to the repository. The individual argument is
a JSON document following the GA4GH JSON schema.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-individual
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-individual registry.db dataset1 HG00096 '{"description": "A description"}'


++++++++++++++++
add-ontology
++++++++++++++++

Adds a new ontology to the repository. The ontology supplied must be a text
file in `OBO format
<http://owlcollab.github.io/oboformat/doc/obo-syntax.html>`_. If you wish to
serve sequence or variant annotations from a repository, a sequence ontology
(SO) instance is required to translate ontology term names held in annotations
to ontology IDs. Sequence ontology definitions can be downloaded from
the `Sequence Ontology site <https://github.com/The-Sequence-Ontology/SO-Ontologies>`_.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-ontology
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-ontology registry.db path/to/so-xp.obo

Adds the sequence ontology ``so-xp.obo`` to the repository using the
default naming rules.

+++++++++++++++
add-variantset
+++++++++++++++

Adds a variant set to a named dataset in a repository. Variant sets are
currently derived from one or more non-overlapping VCF/BCF files which
may be either stored locally or come from a remote URL. Multiple VCF
files can be specified either directly on the command line or by
providing a single directory argument that contains indexed VCF files.
If remote URLs are used then index files in the local file system must be
provided using the ``-I`` option.

.. todo:: Document adding VariantAnnotationSets using the -a option.

.. argparse::
    :module: ga4gh.cli.repomanager
    :func: getRepoManagerParser
    :prog: ga4gh_repo
    :path: add-variantset
    :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-variantset registry.db 1kg 1kgPhase1/ -R NCBI37

Adds a new variant set to the dataset named ``1kg`` in the repository defined
by the registry database ``registry.db`` using the VCF files contained in the
directory ``1kgPhase1``. Note that this directory must also contain the
corresponding indexes for these files. We associate the reference set named
``NCBI37`` with this new variant set. Because we do not provide a ``--name``
argument, a name is automatically generated using the default name generation
rules.

.. todo:: Add a paragraph to the section where we discuss names above
   where we can define the default name generation rules, and link to
   that from here.

.. code-block:: bash

    $ ga4gh_repo add-variantset registry.db 1kg \
        1kgPhase1/chr1.vcf.gz 1kg/chr2.vcf.gz -n phase1-subset -R NCBI37

Like the last example, we add a new variant set to the dataset ``1kg``,
but here we only use the VCFs for chromosomes 1 and 2. We also specify the
name for this new variant set to be ``phase1-subset``.

.. code-block:: bash

    $ ga4gh_repo add-variantset registry.db 1kg \
        --name phase1-subset-remote -R NCBI37 \
        --indexFiles ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.tbi ALL.chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.tbi \
        ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz \
        ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20110521/ALL.chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

This example performs the same task of creating a subset of the phase1
VCFs, but this time we use the remote URL directly and do not keep a
local copy of the VCF file. Because we are using remote URLs to define
the variant set, we have to download a local copy of the corresponding
index files and provide them on the command line using the ``--indexFiles``
option.

+++++++++++++++++
add-readgroupset
+++++++++++++++++

Adds a readgroup set to a named dataset in a repository.  Readgroup sets are
currently derived from a single indexed BAM file, which can be either
stored locally or based on a remote URL. If the readgroup set is based on
a remote URL, then the index file must be stored locally and specified using
the ``--indexFile`` option.

Each readgroup set must be associated with the reference set that it is aligned
to. The ``add-readgroupset`` command first examines the headers of the BAM file
to see if it contains information about references, and then looks for a
reference set with name equal to the genome assembly identifer defined in the
header. (Specifically, we read the ``@SQ`` header line and use the value of the
``AS`` tag as the default reference set name.) If this reference set exists,
then the readgroup set will be associated with it automatically. If it does not
(or we cannot find the appropriate information in the header), then the
``add-readgroupset`` command will fail. In this case, the user must provide the
name of the reference set using the ``--referenceSetName`` option.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-readgroupset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-readgroupset registry.db 1kg \
        path/to/HG00114.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam

Adds a new readgroup set for an indexed 1000 Genomes BAM file stored on the
local file system. The index file follows the usual convention and is stored in
the same directory as the BAM file and has an extra ``.bai`` extension. The
name of the readgroup set is automatically derived from the file name, and the
reference set automatically set from the BAM header.

.. code-block:: bash

    $ ga4gh_repo add-readgroupset registry.db 1kg ga4gh-example-data/HG00096.bam \
        -R GRCh37-subset -n HG0096-subset

Adds a new readgroup set based on a subset of the 1000 genomes reads for the
HG00096 sample from the example data used in the reference server. In this case
we specify that the reference set name ``GRCh37-subset`` be associated with the
readgroup set. We also override the default name generation rules and specify
the name ``HG00096-subset`` for the new readgroup set.

.. code-block:: bash

    $ ga4gh_repo add-readgroupset registry.db 1kg \
        -n HG00114-remote
        -I /path/to/HG00114.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bai
        ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/HG00114/alignment/HG00114.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam

Adds a new readgroups set based on a 1000 genomes BAM directly from the NCBI
FTP server. Because this readgroup set uses a remote FTP URL, we must specify
the location of the ``.bai`` index file on the local file system.

+++++++++++++++++++++++++
init-rnaquantificationset
+++++++++++++++++++++++++

Initializes a rnaquantification set.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: init-rnaquantificationset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo init-rnaquantificationset repo.db rnaseq.db

Initializes the RNA Quantification Set with the filename rnaseq.db.

+++++++++++++++++++++++++
init-rnaquantificationset
+++++++++++++++++++++++++

Initializes a rnaquantification set.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: init-rnaquantificationset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo init-rnaquantificationset repo.db rnaseq.db

Initializes the RNA Quantification Set with the filename rnaseq.db.

+++++++++++++++++++++
add-rnaquantification
+++++++++++++++++++++

Adds a rnaquantification to a RNA quantification set.

RNA quantification formats supported are currently kallisto and RSEM.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-rnaquantification
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-rnaquantification rnaseq.db data.tsv \
             kallisto ga4gh-example-data/registry.db brca1 \
            --bioSampleName HG00096 --featureSetNames gencodev19
            --readGroupSetName HG00096rna --transcript

Adds the data.tsv in kallisto format to the `rnaseq.db` quantification set with
optional fields for associating a quantification with a Feature Set, Read Group
Set, and BioSample.

++++++++++++++++++++++++
add-rnaquantificationset
++++++++++++++++++++++++

When the desired RNA quantification have been added to the set, use this command
to add them to the registry.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-rnaquantificationset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-rnaquantificationset registry.db brca1 rnaseq.db \
        -R hg37 -n rnaseq

Adds the RNA quantification set `rnaseq.db` to the registry under the `brca1`
dataset. The flags set the reference genome to be hg37 and the name of the
set to `rnaseq`.

+++++++++++++++++++++++++++
add-phenotypeassociationset
+++++++++++++++++++++++++++

Adds an rdf object store.  The cancer genome database
Clinical Genomics Knowledge Base http://nif-crawler.neuinfo.org/monarch/ttl/cgd.ttl,
published by the Monarch project, is the supported format for Evidence.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-phenotypeassociationset
   :nodefault:


Examples:

.. code-block:: bash

    $ ga4gh_repo add-phenotypeassociationset registry.db dataset1 /monarch/ttl/cgd.ttl -n cgd


+++++++++++++++
remove-dataset
+++++++++++++++

Removes a dataset from the repository and recursively removes all
objects (ReadGroupSets, VariantSets, etc) within this dataset.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-dataset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-dataset registry.db dataset1

Deletes the dataset with name ``dataset1`` from the repository
represented by ``registry.db``

++++++++++++++++++++
remove-referenceset
++++++++++++++++++++

Removes a reference set from the repository. Attempting
to remove a reference set that is referenced by other objects in the
repository will result in an error.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-referenceset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-referenceset registry.db NCBI37

Deletes the reference set with name ``NCBI37`` from the repository
represented by ``registry.db``

++++++++++++++++
remove-biosample
++++++++++++++++

Removes a biosample from the repository.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-biosample
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-biosample registry.db dataset1 HG00096

Deletes the biosample with name ``HG00096`` in the dataset
``dataset1`` from the repository represented by ``registry.db``

+++++++++++++++++
remove-individual
+++++++++++++++++

Removes an individual from the repository.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-individual
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-individual registry.db dataset1 HG00096

Deletes the individual with name ``HG00096`` in the dataset
``dataset1`` from the repository represented by ``registry.db``

++++++++++++++++++++
remove-ontology
++++++++++++++++++++

Removes an ontology from the repository. Attempting
to remove an ontology that is referenced by other objects in the
repository will result in an error.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-ontology
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-ontology registry.db so-xp

Deletes the ontology with name ``so-xp`` from the repository
represented by ``registry.db``

+++++++++++++++++
remove-variantset
+++++++++++++++++

Removes a variant set from the repository. This also deletes all
associated call sets and variant annotation sets from the repository.

.. argparse::
    :module: ga4gh.cli.repomanager
    :func: getRepoManagerParser
    :prog: ga4gh_repo
    :path: remove-variantset
    :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-variantset registry.db dataset1 phase3-release

Deletes the variant set named ``phase3-release`` from the dataset
named ``dataset1`` from the repository represented by ``registry.db``.

++++++++++++++++++++
remove-readgroupset
++++++++++++++++++++

Removes a read group set from the repository.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-readgroupset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-readgroupset registry.db dataset1 HG00114

Deletes the readgroup set named ``HG00114`` from the dataset named
``dataset1`` from the repository represented by ``registry.db``.

+++++++++++++++++++++++++++
remove-rnaquantificationset
+++++++++++++++++++++++++++

Removes a rna quantification set from the repository.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-rnaquantificationset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-rnaquantificationset registry.db dataset1 ENCFF305LZB

Deletes the rnaquantification set named ``ENCFF305LZB`` from the dataset named
``dataset1`` from the repository represented by ``registry.db``.

++++++++++++++++++++++++++++++
remove-phenotypeassociationset
++++++++++++++++++++++++++++++

Removes an rdf object store.

.. argparse::
   :module: ga4gh.cli.repomanager
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-phenotypeassociationset
   :nodefault:

Examples:

.. code-block:: bash

    $ ga4gh_repo remove-phenotypeassociationset registry.db dataset1  cgd

------------------
Configuration file
------------------

The GA4GH reference server is a `Flask application <http://flask.pocoo.org/>`_
and uses the standard `Flask configuration file mechanisms
<http://flask.pocoo.org/docs/0.10/config/>`_.
Many configuration files will be very simple, and will consist of just
one directive instructing the server where to find the data repository;
example, we might have

.. code-block:: python

    DATA_SOURCE = "/path/to/registry.db"

For production deployments, we shouldn't need to add any more configuration
than this, as the other keys have sensible defaults. However,
all of Flask's `builtin configuration values <http://flask.pocoo.org/docs/0.10/config/>`_
are supported, as well as the extra custom configuration values documented
here.

When debugging deployment issues, it can be very useful to turn on extra debugging
information as follows:

.. code-block:: python

    DEBUG = True

.. warning::

    Debugging should only be used temporarily and not left on by default.
    Running the server with Flask debugging enable is insecure and should
    never be used in a production environment.

++++++++++++++++++++
Configuration Values
++++++++++++++++++++

DEFAULT_PAGE_SIZE
    The default maximum number of values to fill into a page when responding
    to search queries. If a client does not specify a page size in a query,
    this value is used.

MAX_RESPONSE_LENGTH
    The approximate maximum size of the server buffer used when creating
    responses. This is somewhat smaller than the size of the JSON response
    returned to the client. When a client makes a search request with a given
    page size, the server will process this query and incrementally build
    a response until (a) the number of values in the page list is equal
    to the page size; (b) the size of the internal buffer in bytes
    is >= MAX_RESPONSE_LENGTH; or (c) there are no more results left in the
    query.

REQUEST_VALIDATION
    Set this to True to strictly validate all incoming requests to ensure that
    they conform to the protocol. This may result in clients with poor standards
    compliance receiving errors rather than the expected results.

LANDING_MESSAGE_HTML
    The server provides a simple landing page at its root. By setting this
    value to point at a file containing an HTML block element it is possible to
    customize the landing page. This can be helpful to provide support links
    or details about the hosted datasets.

OIDC_PROVIDER
    If this value is provided, then OIDC is configured and SSL is used. It is
    the URI of the OpenID Connect provider, which should return an OIDC
    provider configuration document.

OIDC_REDIRECT_URI
    The URL of the redirect URI for OIDC. This will be something like
    ``https://SERVER_NAME:PORT/oauth2callback``. During testing
    (and particularly in automated tests), if TESTING is True, we can have
    this automatically configured, but this is discouraged in production,
    and fails if TESTING is not True.

OIDC_CLIENT_ID, OIDC_CLIENT_SECRET
    These are the client id and secret arranged with the OIDC provider,
    if client registration is manual (google, for instance). If the provider
    supports automated registration they are not required or used.

OIDC_AUTHZ_ENDPOINT, OIDC_TOKEN_ENDPOINT, OIDC_TOKEN_REV_ENDPOINT
    If the authorization provider has no discovery document available, you can
    set the authorization and token endpoints here.

------------------------
OpenID Connect Providers
------------------------

The server can be configured to use OpenID Connect (OIDC) for authentication.
As an example, here is how one configures it to use Google as the provider.

Go to https://console.developers.google.com/project and in create a project.

.. image:: images/Create_project.png

Navigate to the project -> APIs & auth -> Consent Screen and enter a product
name

.. image:: images/Consent_screen_-_ga4gh.png

Navigate to project -> APIs & auth -> Credentials, and create a new client ID.

.. image:: images/Credentials_-_ga4gh.png

Create the client as follows:

.. image:: images/Credentials_-_ga4gh_2.png

Which will give you the necessary client id and secret. Use these in the OIDC
configuration for the GA4GH server, using the `OIDC_CLIENT_ID` and
`OIDC_CLIENT_SECRET` configuration variables. The Redirect URI should match
the `OIDC_REDIRECT_URI` configuration variable, with the exception that the
redirect URI shown at google does not require a port (but the configuration
variable does)

.. image:: images/Credentials_-_ga4gh_3.png
