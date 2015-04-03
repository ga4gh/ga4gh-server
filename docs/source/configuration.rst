.. _configuration:

*************
Configuration
*************

The GA4GH reference server has two basic elements to its configuration:
the `Data hierarchy`_ and the `Configuration file`_.

--------------
Data hierarchy
--------------

Data is input to the GA4GH server as a directory hierarchy, in which
the structure of data to be served is represented by the file system. For now,
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
                # More BCFs
        /reads/
            readGroupSet1
                sample1.bam
                sample1.bam.bai
                sample2.bam
                sample2.bam.bai
                # More BAMS

------------------
Configuration file
------------------

The GA4GH reference server is a `Flask application <http://flask.pocoo.org/>`_
and uses the standard `Flask configuration file mechanisms
<http://flask.pocoo.org/docs/0.10/config/>`_.
Many configuration files will be very simple, and will consist of just
one directive instructing the server where to look for data; for
example, we might have

.. code-block:: python

    DATA_SOURCE = "/path/to/data/root"

For production deployments, we shouldn't need to add any more configuration
than this, as the all other keys have sensible defaults. However,
all of Flask's `builtin configuration values <http://flask.pocoo.org/docs/0.10/config/>`_
are supported, as well as the extra custom configuration values documented
here.

When debugging deployment issues, it can be very useful to turn on extra debugging
information as follows:

.. code-block:: python

    DEBUG = True

.. warning::

    Debugging should only be used temporarily and not left on by default.

++++++++++++++++++++
Configuration Values
++++++++++++++++++++

DEFAULT_PAGE_SIZE
    The default maximum number of values to fill into a page when responding
    to search queries. If a client does not specify a page size in a query,
    this value is used.

MAX_RESPONSE_LENGTH
    The approximate maximum size of a response sent to a client in bytes. This
    is used to control the amount of memory that the server uses when
    creating responses. When a client makes a search request with a given
    page size, the server will process this query and incrementally build
    a response until (a) the number of values in the page list is equal
    to the page size; (b) the size of the serialised response in bytes
    is >= MAX_RESPONSE_LENGTH; or (c) there are no more results left in the
    query.

REQUEST_VALIDATION
    Set this to True to strictly validate all incoming requests to ensure that
    they conform to the protocol. This may result in clients with poor standards
    compliance receiving errors rather than the expected results.

RESPONSE_VALIDATION
    Set this to True to strictly validate all outgoing responses to ensure
    that they conform to the protocol. This should only be used for development
    purposes.


