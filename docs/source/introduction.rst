
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

*************************************
Configuration file and data hierarchy
*************************************

The GA4GH reference server is a `Flask application <http://flask.pocoo.org/>`_
and uses the standard `Flask configuration file mechanisms
<http://flask.pocoo.org/docs/0.10/config/>`_. An example configuration file
might look like::

    DATA_SOURCE = "/path/to/data/root"
    # TODO other example config

Data is input to the GA4GH server as a directory hierarchy, in which
the structure of data to be served is represented by the filesystem. For now,
we support only one dataset, but this will be generalised to multiple
datasets in later releases. An example data layout might be::

    ga4gh-data/
        /variants/
            variantSet1/
                chr1.vcf.gz
                chr1.vcf.gz.tbi
                chr2.vcf.gz
                chr2.vcf.gz.tbi
                # More VCFs
            variantSet2/
                chr1.bcf
                chr1.bcf.csi
                chr2.bcf
                chr2.bcf.csi
                # More VCFs
        /reads/
            readGroupSet1
                # TODO fill in details for read data.

