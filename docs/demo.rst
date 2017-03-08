.. _demo:

**************
GA4GH API Demo
**************

In this demo, we'll install a copy of the GA4GH reference
implementation and run a local version of the server using some example
data. We then run some example queries on this server using various
different methods to illustrate the basics of the protocol.
The server can, of course, be run on any machine on the network,
but for simplicity we assume that the client and the server are running
on your local machine during this demo.

The instructions for installation here
are not intended to be used in a production deployment. See
the :ref:`installation` section for a detailed guide on production installation.
To run the demo, you will need a working installation of
`Python 2.7 <https://www.python.org/download/releases/2.7/>`_
and also have `virtualenv <https://virtualenv.pypa.io/en/latest/>`_
installed. We also need to have `zlib <http://www.zlib.net/>`_ and
a few other common libraries installed so that we can build some of the
packages that the reference server depends on.

On Debian/Ubuntu, for example, we can install these
packages using:

.. code-block:: bash

    $ sudo apt-get install python-dev python-virtualenv zlib1g-dev libxslt1-dev libffi-dev libssl-dev

On Fedora 22+ (current), the equivalent would be:

.. code-block:: bash

    $ sudo dnf install python-devel python-virtualenv zlib-devel libxslt-devel openssl-devel

First, we create a virtualenv sandbox to isolate the demo from the
rest of the system, and then activate it:

.. code-block:: bash

    $ virtualenv ga4gh-env
    $ source ga4gh-env/bin/activate

Now, install the `ga4gh package <https://pypi.python.org/pypi/ga4gh>`_
from the `Python package index <https://pypi.python.org/pypi>`_. This
will take some time, as some upstream packages will need to be built and
installed.

.. code-block:: bash

    (ga4gh-env) $ pip install ga4gh-server --pre

(Older versions of `pip <https://pip.pypa.io/en/latest/>`_ might not recognise
the ``--pre`` argument; if not, it is safe to remove it.)

Now we can download some example data, which we'll use for our demo:

.. code-block:: bash

    (ga4gh-env) $ wget https://github.com/ga4gh/ga4gh-server/releases/download/data/ga4gh-example-data_4.6.tar
    (ga4gh-env) $ tar -xvf ga4gh-example-data_4.6.tar

After extracting the data, we can then run the ``ga4gh_server`` application:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_server
    * Running on http://127.0.0.1:8000/ (Press CTRL+C to quit)
    * Restarting with stat

(The server is using a default configuration which assumes the
existence of the ``ga4gh-example-data`` directory for simplicity here; see
the :ref:`configuration` section for detailed information on how we configure the
server.) We now have a server running in the foreground. When it receives requests,
it will print out log entries to the terminal. A summary of the server's
configuration and data is available in HTML format at
``http://locahost:8000``, which can be viewed in a web browser.
Leave the server running and open another terminal to complete the
rest of the demo.

To try out the server, we must send some requests to it using the `GA4GH
protocol <http://ga4gh.org/#/api>`_. One way in which we can do this is to
manually create the `JSON <http://json.org/>`_ requests, and send these to the
server using `cURL <http://curl.haxx.se/>`_:

.. code-block:: bash

    $ curl --data '{}' --header 'Content-Type: application/json' \
    http://localhost:8000/datasets/search | jq .


In this example, we used the `search_datasets
<http://ga4gh.org/documentation/api/v0.5.1/ga4gh_api.html#/schema/org.ga4gh.search_datasets>`_
method to ask the server for all the Datasets on the server. It responded
by sending back some JSON, which we piped into the `jq <https://stedolan.github.io/jq/>`_
JSON processor to make it easier to read. We get the following result:

.. code-block:: json

    {
      "nextPageToken": null,
      "datasets": [
        {
          "description": null,
          "name": "1kg-p3-subset",
          "id": "MWtnLXAzLXN1YnNldA=="
        }
      ]
    }

In this example we sent a SearchDatasetsRequest object to the server
and received a SearchDatasetsResponse object in return. This response object
contained one Dataset object, which is contained in the ``datasets`` array.
This approach to interacting with the server is tedious and error prone, as
we have to hand-craft the request objects. It is also quite inconvenient, as
we may have to request many pages of objects to get all the objects
that satisfy our search criteria.

To simplify interacting with the server and to abstract away the low-level
network-level details of the server, we provide a client application.
To try this out, we start another instance of our virtualenv, and then send
the equivalent command using:

.. code-block:: bash

    $ source ga4gh-env/bin/activate
    (ga4gh-env) $ ga4gh_client datasets-search http://localhost:8000

::

    MWtnLXAzLXN1YnNldA==    1kg-p3-subset

The output of this command is a summary of the Datasets on that are present on the
server. We can also get the output in JSON form such that each
object is written on one line:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client datasets-search -O json http://localhost:8000

::

    {"description": null, "name": "1kg-p3-subset", "id": "MWtnLXAzLXN1YnNldA=="}

This format is quite useful for larger queries, and can be piped into jq
to extract fields of interest, pretty printing and so on.

We can perform similar queries for variant data using the
`search_variants
<http://ga4gh.org/documentation/api/v0.5.1/ga4gh_api.html#/schema/org.ga4gh.search_variants>`_
API call. First, we find the IDs of the VariantSets on the server using the
`search_variant_sets
<http://ga4gh.org/documentation/api/v0.5.1/ga4gh_api.html#/schema/org.ga4gh.search_variant_sets>`_
method:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client variantsets-search http://localhost:8000

::

    MWtnLXAzLXN1YnNldDptdm5jYWxs    mvncall

This tells us that we have one VariantSet on the server, with ID ``MWtnLXAzLXN1YnNldDptdm5jYWxs``
and name ``mvncall``. We can then search for variants overlapping a given interval in a VariantSet
as follows:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client variants-search http://localhost:8000 \
    --referenceName=1 --start=45000 --end=50000

The output of the client program is a summary of the data received in a
free text form. This is not intended to be used as the input to other
programs, and is simply a data exploration tool for users.
To really *use* our data, we should use a GA4GH client library.

Part of the GA4GH reference implementation is a client library. This makes sending requests to the server and using the
responses very easy. For example, to run the same query as we
performed above, we can use the following code:

.. code-block:: python

    from __future__ import print_function

    from ga4gh.client import client

    httpClient = client.HttpClient("http://localhost:8000")
    # Get the datasets on the server.
    datasets = list(httpClient.search_datasets())
    # Get the variantSets in the first dataset.
    variantSets = list(httpClient.search_variant_sets(
        dataset_id=datasets[0].id))
    # Now get the variants in the interval [45000, 50000) on chromosome 1
    # in the first variantSet.
    iterator = httpClient.search_variants(
        variant_set_id=variantSets[0].id,
        reference_name="1", start=45000, end=50000)
    for variant in iterator:
        print(
            variant.reference_name, variant.start, variant.end,
            variant.reference_bases, variant.alternate_bases, sep="\t")


If we save this script as ``ga4gh-demo.py`` we can then run it
using:

.. code-block:: bash

    (ga4gh-env) $ python ga4gh-demo.py


Host the 1000 Genomes VCF
=============================

The GA4GH reference server uses a registry of files and URLs to
populate its data repository. In this tutorial we will use the
command-line client to create a registry similar to that used by
1kgenomes.ga4gh.org. Your system should have samtools installed, and at
least 30GB to host the VCF and reference sets.

Repo administrator CLI
----------------------

The CLI has methods for adding and removing Feature Sets, Read Group
Sets, Variant Sets, etc. Before we can begin adding files we must first
initialize an empty registry database. The directory that this database
is in should be readable and writable by the current user, as well as the
user running the server.

.. code-block:: bash

    $ ga4gh_repo init registry.db

This command will create a file ``registry.db`` in the current working
directory. This file should stay relatively small (a few MB for
thousands of files).

Now we will add a dataset to the registry, which is a logical container
for the genomics data we will later add. You can optionally provide a
description using the ``--description`` flag.

.. code-block:: bash

    $ ga4gh_repo add-dataset registry.db 1kgenomes \
        --description "Variants from the 1000 Genomes project and GENCODE genes annotations"

Add a Reference Set
-------------------

It is possible for a server to host multiple reference assemblies. Here
we will go through all the steps of downloading and adding the FASTA
used for the 1000 Genomes VCF.

.. code-block:: bash

    $ wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

This file is provided in ``.gz`` format, which we will decompress, and
then with samtools installed on the system, recompress it using
``bgzip``.

.. code-block:: bash

    $ gunzip hs37d5.fa.gz
    $ bgzip hs37d5.fa

This may take a few minutes depending on your system as this file is
around 3GB. Next, we will add the reference set.

.. code-block:: bash

    $ ga4gh_repo add-referenceset registry.db /full/path/to/hs37d5.fa.gz \
      -d “NCBI37 assembly of the human genome” \
      --species '{"id": "9606", "term": "Homo sapiens", "source_name": "NCBI", "source_version: "1.0"}' \
      --name NCBI37 \
      --sourceUri "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"

A number of optional command line flags have been added. We will be
referring to the name of this reference set ``NCBI37`` when we later add
the variant set.

Add an ontology
---------------

Ontologies provide a source for parsing variant annotations, as well as
organizing feature types into ontology terms. A `sequence ontology
<http://www.sequenceontology.org/>`_ instance must be added to the repository
to translate ontology term names in sequence and variant annotations to IDs.
Sequence ontology definitions can be downloaded from the `Sequence Ontology
site <https://github.com/The-Sequence-Ontology/SO-Ontologies>`_.

.. code-block:: bash

    $ wget https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/so-xp-dec.obo
    $ ga4gh_repo add-ontology registry.db /full/path/to/so-xp.obo -n so-xp

Add sequence annotations
------------------------

The GENCODE Genes dataset provides annotations for features on the
reference assembly. The server uses a custom storage format for sequence
annotations, you can download a prepared set
`here <https://ga4ghstore.blob.core.windows.net/testing/gencode_v24lift37.db>`__.
It can be added to the registry using the following command. Notice
we have told the registry to associate the reference set added above
with these annotations.

.. code-block:: bash

    $ wget https://ga4ghstore.blob.core.windows.net/testing/gencode_v24lift37.db
    $ ga4gh_repo add-featureset registry.db 1kgenomes /full/path/to/gencode.v24lift37.annotation.db \
        --referenceSetName NCBI37 --ontologyName so-xp


Add the 1000 Genomes VCFs
--------------------------

The 1000 Genomes are publicly available on the EBI server. This
command uses ``wget`` to download the "release" VCFs to a directory named
release.

.. code-block:: bash

    $ wget -m ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ -nd -P release -l 1
    $ rm release/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz

These files are already compressed and indexed. For the server to make use
of the files in this directory we must move the `wgs` file, since it covers
chromosomes that are represented elsewhere and overlapping VCF are not
currently supported. This file could be added as a separate variant set.

We can now add the directory to the registry using the following command.
Again, notice we have referred to the reference set by name.

.. code-block:: bash

    $ ga4gh_repo add-variantset registry.db 1kgenomes /full/path/to/release/ \
        --name phase3-release --referenceSetName NCBI37

Add a BAM as a Read Group Set
-----------------------------

Read Group Sets are the logical equivalent to BAM files within the
server. We will add a BAM hosting by the 1000 Genomes S3 bucket.
We will first download the index and then add it to the registry.

.. code-block:: bash

    $ wget http://s3.amazonaws.com/1000genomes/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bai
    $ ga4gh_repo add-readgroupset registry.db 1kgenomes \
        -I HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bai \
        --referenceSetName NCBI37 \
        http://s3.amazonaws.com/1000genomes/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam \

This might take a moment as some metadata about the file will be
retrieved from S3.

Start the server
----------------

Assuming you have set up your server to run using the registry file just
created, you can now start or restart the server to see the newly added
data. If the server is running via apache issue
``sudo service apache2 restart``. You can then visit the landing page of
the running server to see the newly added data.


Use the client package
=============================

If you only want to use the client and don't need the server functionality,
there is a seperate pypi package, `ga4gh-client
<https://pypi.python.org/pypi/ga4gh-client>`_, which includes only the
client.  It is also much quicker to install.  To install, simply run:

.. code-block:: bash

    (ga4gh-env) $ pip install --pre ga4gh_client

This installs the ``ga4gh_client`` command line program, which provides
identical functionality to the ``ga4gh_client`` which is installed via the
``ga4gh`` package:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client datasets-search http://1kgenomes.ga4gh.org

Installing the ``ga4gh_client`` package also gives you access to the
client's libraries for use in your own programs:

.. code-block:: python

    >>> from ga4gh.client import client
    >>> client.HttpClient
    <class 'ga4gh_client.client.HttpClient'>

For more examples of using the GA4GH client visit 
`this iPython notebook <https://github.com/BD2KGenomics/bioapi-examples/blob/master/python_notebooks/1kg.ipynb>`_.



OIDC Demonstration
==================

If we want authentication, we must have an OIDC authentication provider.
One can be found in ``oidc-provider``, and run with the ``run.sh`` script.
We can then use this with the ``LocalOidConfig`` server configuration. So:

.. code-block:: bash

  $ cd oidc-provider && ./run.sh

In another shell on the same machine

.. code-block:: bash

  $ python server_dev.py -c LocalOidConfig

Make sure you know the hostname the server is running on. It can be found with

.. code-block:: bash

  $ python -c 'import socket; print socket.gethostname()'

With a web browser, go to ``https://<server hostname>:<server port>``. You may
need to accept the security warnings as there are probably self-signed
certificates. You will be taken through an authentication flow. When asked
for a username and password, try ``upper`` and ``crust``. You will find
yourself back at the ga4gh server homepage. On the homepage will be a
'session token' This is the key to access the server with the client tool
as follows:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client --key <key from homepage> variantsets-search https://localhost:8000/current
    MWtnLXAzLXN1YnNldDptdm5jYWxs    mvncall
