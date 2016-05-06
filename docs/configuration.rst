.. _configuration:

*************
Configuration
*************

The GA4GH reference server has two basic elements to its configuration:
the `Data repository`_ and the `Configuration file`_.

---------------
Data repository
---------------

The repository in the GA4GH reference server defines how your data is organised. The
repository itself is a SQLite database, which contains information about your
datasets, reference sets and so on. Bulk data (such as variants and reads)
is not stored in database, but instead accessed directly from the primary
data files at run time. The locations of these data files is entirely up
to the administrator.

The repository manager provides an administration interface to the the data
repository. It can be accessed via ``ga4gh_repo`` (or ``python repo_dev.py`` if
developing). The ``ga4gh_repo`` commands provides a suite of commands to
help administer a GA4GH data repository.

+++++++
init
+++++++

.. todo:: Description of init.

.. argparse::
    :module: ga4gh.cli
    :func: getRepoManagerParser
    :prog: ga4gh_repo
    :path: init
    :nodefault:


**Examples:**

.. code-block:: bash

    $ ga4gh_repo init path/to/repo.db

+++++++
verify
+++++++

Performs some consistency checks on the given data repository to ensure it is
well-formed.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: verify
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo verify path/to/repo.db

+++++++
list
+++++++

Lists the contents of the given data repository.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: list
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo list path/to/repo.db

+++++++++++
add-dataset
+++++++++++

Creates a dataset in the given repository with a given name.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-dataset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-dataset path/to/repo.db aDataset

+++++++++++++++
remove-dataset
+++++++++++++++

Destroys a dataset in the given repository with a given name.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-dataset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-dataset path/to/repo.db aDataset

++++++++++++++++
add-referenceset
++++++++++++++++

Adds a given reference set file to a given data repository.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-referenceset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-referenceset path/to/repo.db path/to/aReferenceSet.fa.gz

++++++++++++++++++++
remove-referenceset
++++++++++++++++++++

Removes a given reference set from a given data repository.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-referenceset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-referenceset path/to/repo.db aReferenceSet

++++++++++++++++
add-ontology
++++++++++++++++

Adds an Ontology Map, which maps identifiers to ontology terms, to
the repository. Ontology maps are tab delimited files with an
identifier/term pair per row.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-ontology
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-ontology path/to/repo.db path/to/aOntoMap.txt

++++++++++++++++++++
remove-ontology
++++++++++++++++++++

Removes a given Ontology Map from a given data repository.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-ontology
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-ontology path/to/repo.db aOntoMap

+++++++++++++++++
add-readgroupset
+++++++++++++++++

Adds a given read group set file to a given data repository and dataset.  The
file must have the extension ``.bam``.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: add-readgroupset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-readgroupset path/to/repo.db aDataset path/to/aReadGroupSet.bam

++++++++++++++++++++
remove-readgroupset
++++++++++++++++++++

Removes a read group set from a given data repository and dataset.

.. argparse::
   :module: ga4gh.cli
   :func: getRepoManagerParser
   :prog: ga4gh_repo
   :path: remove-readgroupset
   :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-readgroupset path/to/repo.db aDataset aReadGroupSet

+++++++++++++++
add-variantset
+++++++++++++++

Adds a variant set directory to a given data repository and dataset.  The
directory should contain file(s) with extension ``.vcf.gz``. If a variant set
is annotated it will be added as both a variant set and a variant annotation set.

.. argparse::
    :module: ga4gh.cli
    :func: getRepoManagerParser
    :prog: ga4gh_repo
    :path: add-variantset
    :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo add-variantset path/to/repo.db aDataset path/to/aVariantSet

+++++++++++++++++
remove-variantset
+++++++++++++++++

Removes a variant set from a given data repository and dataset.

.. argparse::
    :module: ga4gh.cli
    :func: getRepoManagerParser
    :prog: ga4gh_repo
    :path: remove-variantset
    :nodefault:

**Examples:**

.. code-block:: bash

    $ ga4gh_repo remove-variantset path/to/repo.db aDataset aVariantSet

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

    DATA_SOURCE = "/path/to/repo.db"

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
