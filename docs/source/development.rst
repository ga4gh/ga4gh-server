.. _development:

-----------
Development
-----------

A development installation of the GA4GH reference implementation is a local
copy of the ``server`` repo, along with all of the tools required for development.
Please ensure that all the system requirements (as listed above) are installed, and
clone a local copy of the repo. Install all of the required Python libraries
into your Python user installation::

    $ pip install -r requirements.txt --user

All of the command line interface utilities have local scripts
that simplify development: for example, we can run the local version of the
``ga2sam`` program by using::

    $ python ga2sam_dev.py

To run the server locally in development mode, we can use the ``server_dev.py``
script, e.g.::

    $ python server_dev.py

will run a server using the default configuration. This default configuration
expects a data hierarchy to exist in the ``ga4gh-example-data`` directory.
This default configuration can be changed by providing a (fully qualified)
path to a configuration file (see the :ref:`configuration`
section for details).

++++++
Layout
++++++

The code for the project is held in the ``ga4gh`` package, which corresponds to
the ``ga4gh`` directory in the project root. Within this package, the
functionality is split between the ``client``, ``server``, ``protocol`` and
``cli`` modules.  The ``cli`` module contains the definitions for the
``ga4gh_client`` and ``ga4gh_server`` programs.
