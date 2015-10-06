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
