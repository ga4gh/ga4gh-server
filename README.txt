
.. image:: http://genomicsandhealth.org/files/logo_ga.png

==============================
GA4GH Reference Implementation
==============================

This is a prototype for the GA4GH reference client and
server applications. It is under heavy development, and many aspects of
the layout and APIs will change as requirements are better understood.
If you would like to help, please check out our list of
`issues <https://github.com/ga4gh/server/issues>`_!

Our aims for this implementation are:

**Simplicity/clarity**
    The main goal of this implementation is to provide an easy to understand
    and maintain implementation of the GA4GH API. Design choices
    are driven by the goal of making the code as easy to understand as
    possible, with performance being of secondary importance. With that
    being said, it should be possible to provide a functional implementation
    that is useful in many cases where the extremes of scale are not
    important.

**Portability**
    The code is written in Python for maximum portability, and it
    should be possible to run on any modern computer/operating system (Windows
    compatibility should be possible, although this has not been tested). Our coding
    guidelines specify using a subset of Python 3 which is backwards compatible with Python 2
    following the current `best practices <http://python-future.org/compatible_idioms.html>`_.
    The project currently does not yet support Python 3, as support for it is lacking in several
    packages that we depend on. However, our eventual goal is to support both Python 2
    and 3.

**Ease of use**
    The code follows the `Python Packaging User Guide
    <http://python-packaging-user-guide.readthedocs.org/en/latest/>`_.
    Specifically, pip is used to handle python package dependencies (see below
    for details). This allows for easy installation of the ``ga4gh`` reference code
    across a range of operating systems.


*******************************
Installing the reference server
*******************************

Prerequisites:

* Python 2.7,
* Berkeley DB together with include and lib files (version 4.8 or higher),
* Virtualenv (or another python sandboxing tool) is highly recommended.

General installation procedure:

* Install Berkeley DB (version 4.8 or higher) using your system's preferred
  package manager, see the `wormtable help page
  <https://pypi.python.org/pypi/wormtable>`_ for platform-specific details.

* (On MacOS X, make sure the LDFLAGS and CFLAGS environment variables are set to
  include the lib and include directories for the Berkeley DB install of your choice.
  The wormtable help page cited above provides more detailed instructions, or
  see the `Appendix`_ to this README for an example install on that platform.)

* Create a python sandbox directory using virtualenv, preferably
  *not* inside the ga4gh server directory. For an good introduction
  to using virtualenv, see the `Python Guide page
  <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_.
  On some systems, you may need to specify the --no-site-packages
  option to ensure a clean dependency install. For example, to
  create a sandbox named ``testenv``::

  $ virtualenv --no-site-packages testenv

* Make the virtualenv sandbox created above active::

  $ source testenv/bin/activate

* cd to the ga4gh server directory, and load the dependencies via pip::

  $ cd [your ga4gh server directory]
  $ pip install -r requirements.txt

* Finally, run the install script, and run nosetests to confirm the install::

  $ python setup.py install
  $ nosetests

A successfull install should result in a clean run of all the tests,
resulting in a line of dots followed by ``OK``. If this still isn't working,
you may want to check the `Appendix`_ at the end of this README for
system-specific installation examples.

********************************
Serving variants from a VCF file
********************************

Two implementations of the variants API are available that can serve data based
on existing VCF files. These backends are based on tabix and `wormtable
<http://www.biomedcentral.com/1471-2105/14/356>`_, which is a Python library to
handle large scale tabular data. See `Wormtable backend`_ for instructions on
serving VCF data from the GA4GH API.

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

To run the server on this example dataset, follow the steps on
installing the server, then download and unpack the example data ::

    $ wget http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar.gz
    $ tar -zxvf ga4gh-example-data.tar.gz

An easier way to download and upack the data is to run the following
script, which will do these steps for you::

    $ python scripts/update_data.py

You can now run the server, telling it to serve variants from the sets in
the downloaded datafile::

    $ ga4gh_server ga4gh-example-data

To run queries against this server, we can use the ``ga4gh_client`` program;
for example, here we run the ``variants/search`` method over the
``1000g_2013.wt`` variant set, where the reference name is ``1``
and we only want calls returned for call set ID HG03279::

    $ ga4gh_client variants-search http://localhost:8000/v0.5.1 -V 1000g_2013.wt -r 1 -c HG03279 | less -S

We can also query against the *variant name*; here we return the variant that
has variant name ``rs75454623``::

    $ ga4gh_client variants-search http://localhost:8000/v0.5.1 -V 1000g_2013.wt -r 1 -n rs75454623  | less -S

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

*************
Tabix backend
*************

The tabix backend allows us to serve variants from an arbitrary VCF file.  The
VCF file must first be indexed with `tabix
<http://samtools.sourceforge.net/tabix.shtml>`_.  Many projects, including the
`1000 genomes project
<http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>`_, release files
with tabix indices already precomputed.  This backend can serve such datasets
without any preprocessing via the command::

    $ ga4gh_server DATADIR

where DATADIR is a directory that contains subdirectories of tabix-indexed VCF
file(s).  There cannot be more than one VCF file in any subdirectory that has
data for the same reference contig.

******
Layout
******

The code for the project is held in the ``ga4gh`` package, which corresponds to
the ``ga4gh`` directory in the project root. Within this package, the
functionality is split between the ``client``, ``server``, ``protocol`` and
``cli`` modules.  The ``cli`` module contains the definitions for the
``ga4gh_client`` and ``ga4gh_server`` programs.

For development purposes, it is useful to be able to run the command line
programs directly without installing them. To do this, use the
``server_dev.py`` and ``client_dev.py`` scripts. (These are just shims to
facilitate development, and are not intended to be distributed.  The
distributed versions of the programs are packaged using the setuptools
``entry_point`` key word; see ``setup.py`` for details). For example, the run
the server command simply run::

    $ python server_dev.py
    usage: server_dev.py [-h] [--port PORT] [--verbose] {help,wormtable,tabix} ...
    server_dev.py: error: too few arguments

++++++++++++
Coding style
++++++++++++

The code follows the guidelines of `PEP 8
<http://legacy.python.org/dev/peps/pep-0008>`_ in most cases. The only notable
difference is the use of camel case over underscore delimited identifiers; this
is done for consistency with the GA4GH API. Code should be checked for compliance
using the `pep8 <https://pypi.python.org/pypi/pep8>`_ tool.


**********
Deployment
**********

*TODO* Give simple instructions for deploying the server on common platforms
like Apache and Nginx.

Configuration parameters are specified in the file ga4gh/server/config.py;
they can be overridden by setting the absolute path of a file containing
new values in the environment variable GA4GH_CONFIGURATION.

********
Appendix
********

**system specific install examples**

MacOS X (with MacPorts)::

  $ sudo port install db48
  $ export CFLAGS=-I/opt/local/include/db48/  LDFLAGS=-L/opt/local/lib/db48/
  $ cd [some working directory outside the ga4gh server directory tree]
  $ virtualenv --no-site-packages testenv
  $ source testenv/bin/activate
  $ cd [your ga4gh server directory]
  $ pip install -r requirements.txt
  $ python setup.py install
  $ nosetests

*TODO* Append examples of installs (using package managers if possible, no dependency
installs from source) on the target platform of your choice.


