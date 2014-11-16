==============================
GA4GH Reference Implementation
==============================

A reference implementation of the APIs defined in the schemas repository.

*************************
Initial skeleton overview
*************************

This is a proposed skeleton layout for the GA4GH reference client and
server applications. As such, nothing is finalised and all aspects of
the design and implementation are open for discussion and debate. The overall
goals of the project are:

Simplicity/clarity
    The main goal of this implementation is to provide an easy to understand
    and maintain implementation of the GA4GH API. Design choices
    are driven by the goal of making the code as easy to understand as
    possible, with performance being of secondary importance. With that
    being said, it should be possible to provide a functional implementation
    that is useful in many cases where the extremes of scale are not
    important.

Portability
    The code is written in Python for maximum portability, and it
    should be possible to run on any modern computer/operating system (Windows
    compatibility should be possible, although this has not been tested).  We use a
    subset of Python 3 which is backwards compatible with Python 2 following the
    current `best practices <http://python-future.org/compatible_idioms.html>`_.
    In this way, we fully support both Python 2 and 3.

Ease of use
    The code follows the `Python Packaging User Guide
    <http://python-packaging-user-guide.readthedocs.org/en/latest/>`_. This will
    make installing the ``ga4gh`` reference code very easy across a range of
    operating systems.



*************
Trying it out
*************

The project is designed to be published as a `PyPI <https://pypi.python.org/pypi>`_
package, so ultimately installing the reference client and server programs
should be as easy as::

    $ pip install ga4gh

However, the code is currently only a proposal, so it has not been uploaded to
the Python package index. The best way to try out the code right now is to
use `virtualenv <http://virtualenv.readthedocs.org/en/latest/>`_. After cloning
the git repo, and changing to the project directory, do the following::

    $ virtualenv testenv
    $ source testenv/bin/activate
    $ python setup.py install

This should install the ``ga4gh_server`` and ``ga4gh_client`` scripts into the
virtualenv and update your ``PATH`` so that they are available. When you have
finished trying out the programs you can leave the virtualenv using::

    $ deactivate

The virtualenv can be restarted at any time, and can also be deleted
when you no longer need it.

********************************
Serving variants from a VCF file
********************************

Two implementations of the variants API is available that can serve data based
on existing VCF files. This backends are based on tabix and `wormtable
<http://www.biomedcentral.com/1471-2105/14/356>`_, which is a Python library
to handle large scale tabular data. See `Wormtable backend`_ for instructions
on serving VCF data from the GA4GH API.

*****************
Wormtable backend
*****************

The wormtable backend allows us to serve variants from an arbitrary VCF file.
The VCF file must first be converted to wormtable format using the ``vcf2wt``
utility (the `wormtable tutorial
<http://pythonhosted.org/wormtable/tutorial.html>`_ discusses this process).
A subset (1000 rows for each chromosome) of the 1000 Genomes VCF data (20110521
and 20130502 releases) has been prepared and converted to wormtable format
and made available `here <http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar.gz>`_.
See `Converting 1000G data`_ for more information on converting 1000 genomes
data into wormtable format.

To run the server on this example dataset, create a virtualenv and install
wormtable::

    $ virtualenv testenv
    $ source testenv/bin/activate
    $ pip install wormtable

See the `wormtable PyPI page <https://pypi.python.org/pypi/wormtable>`_ for
detailed instructions on installing wormtable and its dependencies.

Now, download and unpack the example data, ::

    $ wget http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar.gz
    $ tar -zxvf ga4gh-example-data.tar.gz

and install the client and server scripts into the virtualenv (assuming
you are in the project root directory)::

    $ python setup.py install

We can now run the server, telling it to serve variants from the sets in
the downloaded datafile::

    $ ga4gh_server wormtable ga4gh-example-data

To run queries against this server, we can use the ``ga4gh_client`` program;
for example, here we run the ``variants/search`` method over the
``1000g_2013`` variant set, where the reference name is ``1``, the end coordinate
is 60000 and we only want calls returned for call set ID HG03279::

    $ ga4gh_client variants-search http://localhost:8000 1000g_2013 -r1 -e 60000 -c HG03279 | less -S

We can also query against the *variant name*; here we return the variant that
has variant name ``rs75454623``::

    $ ga4gh_client variants-search http://localhost:8000 1000g_2013 -r1 -e 60000 -n rs75454623  | less -S

+++++++++++++++++++++
Converting 1000G data
+++++++++++++++++++++

To duplicate the data for the above example, we must first create VCF files
that contain the entire variant set of interest. The VCF files for the set
mentioned above have been made `available
<http://www.well.ox.ac.uk/~jk/ga4gh-example-source.tar.gz>`_. After downloading
and extracting these files, we can build the wormtable using ``vcf2wt``::

    $ vcf2wt 1000g_2013-subset.vcf -s schema-1000g_2013.xml -t 1000g_2013

Schemas for the 2011 and 2013 1000G files have been provided as these do a
more compact job of storing the data than the default auto-generated schemas.
We must also truncate and remove some columns because of a current limitation
in the length of strings that wormtable can handle.
After building the table, we must create indexes on the ``POS`` and ``ID`` columns::

    $ wtadmin add 1000g_2013 CHROM+POS
    $ wtadmin add 1000g_2013 CHROM+ID

The ``wtadmin`` program supports several
commands to administer and examine the dataset; see ``wtadmin help`` for details.
These commands and schemas also work for the full 1000G data; however, it is
important to specify a sufficiently large `cache size
<http://pythonhosted.org/wormtable/performance.html#cache-tuning>`_ when
building and indexing such large tables.

*****************
Tabix backend
*****************

The tabix backend allows us to serve variants from an arbitrary VCF file.
The VCF file must first be indexed with `tabix <http://samtools.sourceforge.net/tabix.shtml>`_.
Many projects, including the `1000 genomes project
<http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>`_, release files with tabix
indicies already precomputed.  This backend can serve such datasets without any
preprocessing via the command:

    $ python ga4gh/scripts/server.py tabix DATADIR

where DATADIR is a directory that contains folders of tabix-indexed VCF file(s).  There cannot
be more than one VCF file in any subdirectory that has data for the same reference contig.

******
Layout
******

The code for the project is held in the ``ga4gh`` package, which corresponds
to the ``ga4gh`` directory in the project root. Within this package,
the functionality is split between the ``client``, ``server`` and
``protocol`` modules. There is also a subpackage called ``scripts``
which holds the code defining the command line interfaces for the
``ga4gh_client`` and ``ga4gh_server`` programs.

For development purposes, it is useful to be able to run the command
line programs directly without installing them. To do this, make hard links
to the files in the scripts directory to the project root and run them
from there; e.g::

    $ ln ga4gh/scripts/server.py .
    $ python server.py
    usage: server.py [-h] [--port PORT] [--verbose] {help,simulate} ...
    server.py: error: too few arguments

++++++++++++
Coding style
++++++++++++

The code follows the guidelines of `PEP 8
<http://legacy.python.org/dev/peps/pep-0008>`_ in most cases. The only notable
difference is the use of camel case over underscore delimited identifiers; this
is done for consistency with the GA4GH API. The code was checked for compliance
using the `pep8 <https://pypi.python.org/pypi/pep8>`_ tool.

