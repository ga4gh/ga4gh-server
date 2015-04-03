.. _development:

-----------
Development
-----------

Thanks for your interest in helping us develop the GA4GH reference
implementation! There are lots of ways to contribute, and it's easy
to get up and running. This page should provide the basic information
required to get started; if you encounter any difficulties
`please let us know <https://github.com/ga4gh/server/issues>`_

.. warning::

    This guide is a work in progress, and is incomplete.

***********************
Development environment
***********************

We need a development Python 2.7 installation, Git, and some basic
libraries. On Debian or Ubuntu, we can install these using

.. code-block:: bash

    $ sudo apt-get install python-dev zlib1g-dev git

.. note::
    TODO: Document this basic step for other platforms? We definitely want
    to tell people how to do this with Brew or ports on a Mac.

If you don't have admin access to your machine, please contact your system
administrator, and ask them to install the development version of Python 2.7
and the development headers for `zlib <http://www.zlib.net/>`_.

Once these basic prerequisites are in place, we can then bootstrap our
local Python installation so that we have all of the packages we require
and we can keep them up to date. Because we use the functionality
of the recent versions of ``pip`` and other tools, it is important to
use our own version of it and not any older versions that may be
already on the system.

.. code-block:: bash

    $ wget https://bootstrap.pypa.io/get-pip.py
    $ python get-pip.py --user

This creates a `user specific <https://www.python.org/dev/peps/pep-0370/>`_
site-packages installation for Python, which is based in your ``~/.local``
directory. This means that you can now install any Python packages you like
without needing to either bother your sysadmin or worry about breaking your
system Python installation. To use this, you need to add the newly installed
version of ``pip`` to your ``PATH``. This can be done by adding something
like

.. code-block:: bash

    export PATH=$HOME/.local/bin:$PATH

to your ``~/.bashrc`` file. (This will be slightly different if you use
another shell like ``csh`` or ``zsh``.)

We then need to activate this configuration by logging out, and logging back in.
Then, test this by running:

.. code-block:: bash

    $ pip --version
    pip 6.0.8 from /home/username/.local/lib/python2.7/site-packages (python 2.7)

We are now ready to start developing!

***************
GitHub workflow
***************

First, go to https://github.com/ga4gh/server and click on the 'Fork'
button in the top right-hand corner. This will allow you to create
your own private fork of the server project where you can work.
See the `GitHub documentation <https://help.github.com/articles/fork-a-repo/>`_
for help on forking repositories.
Once you have created your own fork on GitHub, you'll need to clone a
local copy of this repo. This might look something like:

.. code-block:: bash

    $ git clone git@github.com:username/server.git

We can then install all of the packages that we need for developing the
GA4GH reference server:

.. code-block:: bash

    $ cd server
    $ pip install -r requirements.txt --user

This will take a little time as the libraries that we require are
fetched from PyPI and built.

It is also important to set up an
`upstream remote <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`_
for your repo so that you can sync up with the changes that other people
are making:

.. code-block:: bash

    $ git remote add upstream https://github.com/ga4gh/server.git


**TODO**

- Walk a new developer through the process of fetching ``upstream/master``
  and creating a new topic branch.
- Walk them through pushing to their fork and making a pull request.
- Explain why we do a lot of rebasing, and walk them through the
  process of squashing commits.


************
Contributing
************

See the files ``CONTRIBUTING.md`` and ``STYLE.md`` for an overview of
the processes for contributing code and the style guidelines that we
use.


*********************
Development utilities
*********************

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

************
Organisation
************

The code for the project is held in the ``ga4gh`` package, which corresponds to
the ``ga4gh`` directory in the project root. Within this package, the
functionality is split between the ``client``, ``server``, ``protocol`` and
``cli`` modules.  The ``cli`` module contains the definitions for the
``ga4gh_client`` and ``ga4gh_server`` programs.
