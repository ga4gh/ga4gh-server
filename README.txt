
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

************
Installation
************

There are four different types of installation that we deal with here:
`Deployment`_, `Deployment with Docker`_, `Client tools`_ and `Development`_ installations. A
deployment installation is a production server, usually using Apache
or another web server on a dedicated machine. A client tools installation
creates a sandbox in which a user can easily try out the GA4GH client
utilities, and run queries against arbitrary servers (specifically,
any server running the correct version of the GA4GH API; not necessarily
this implementation). Finally, a development installation is a local
installation used for either development of GA4GH reference server itself,
or client applications depending on the reference Python client libraries.

----------
Deployment
----------

To deploy on Apache on Debian/Ubuntu platforms, do the following.

- Install some basic pre-requisite packages::

  $ sudo apt-get install python-dev zlib1g-dev

- Install Apache and mod_wsgi, and enable mod_wsgi::

  $ sudo apt-get install apache2 libapache2-mod-wsgi
  $ sudo a2enmod wsgi

- Create the python egg cache directory, and make it writable by
  www-data::

  $ sudo mkdir /var/cache/apache2/python-egg-cache
  $ sudo chown www-data:www-data /var/cache/apache2/python-egg-cache/

- Create a directory to hold the GA4GH server code, configuration
  and data. For convenience, we make this owned by the current user
  (but make sure all the files are world-readable).::

  $ sudo mkdir /srv/ga4gh
  $ sudo chown $USER /srv/ga4gh
  $ cd /srv/ga4gh

- Make a virtualenv, and install the ga4gh package::

  $ virtualenv ga4gh-server-env
  $ source ga4gh-server-env/bin/activate
  $ pip install --pre ga4gh  # We need the --pre because ga4gh is pre-release
  $ deactivate

- Download and unpack the example data::

  $ wget http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar.gz
  $ tar -zxf ga4gh-example-data.tar.gz

- Create the WSGI file at ``/srv/ga4gh/application.wsgi`` and write the following
  contents::

    from ga4gh.frontend import app as application
    import ga4gh.frontend as frontend
    frontend.configure("/srv/ga4gh/config.py")

- Create the configuration file at ``/srv/ga4gh/config.py``, and write the
  following contents::

    DATA_SOURCE = "/srv/ga4gh/ga4gh-example-data"

- Configure Apache. Edit the file ``/etc/apache2/sites-enabled/000-default.conf``
  and insert the following contents towards the end of the file
  (*within* the ``<VirtualHost:80>...</VirtualHost>`` block)::

    WSGIDaemonProcess ga4gh python-path=/srv/ga4gh/ga4gh-server-env/lib/python2.7/site-packages python-eggs=/var/cache/apache2/python-egg-cache
    WSGIScriptAlias /ga4gh /srv/ga4gh/application.wsgi

    <Directory /srv/ga4gh>
        WSGIProcessGroup ga4gh
        WSGIApplicationGroup %{GLOBAL}
        Require all granted
    </Directory>

- Restart Apache::

  $ sudo service apache2 restart

- Test the installation by pointing a web-browser at the root URL; for example,
  to test on the installation server use::

    $ links http://localhost/ga4gh

  To test the API, make a development installation and try running some
  `Example client queries`_ against the server

These instructions are just one way in which we can achieve the same thing.
There are any number of different ways in which we can set up a WSGI
application under Apache, which may be preferable in different installations.
(In particular, the Apache configuration here may be specific to
Ubuntu 14.04, where this was tested.)
See the `mod_wsgi documentation <https://code.google.com/p/modwsgi/>` for
more details. These instructions are also specific to Debian/Ubuntu and
different commands and directory structures will be required on
different platforms.

The server can be deployed on any WSGI compliant web server. See the
instructions in the `Flask documentation
<http://flask.pocoo.org/docs/0.10/deploying/>` for more details on
how to deploy on various other servers.

+++++++++++++++
Troubleshooting
+++++++++++++++

If you are encountering difficulties getting the above to work, it is helpful
to turn on debugging output. Do this by adding the following line to your
config file::

    DEBUG = True

When an error occurs, the details of this will then be printed to the web server's
error log (in Apache on Debian/Ubuntu, for example, this is ``/var/log/apache2/error.log``).


----------------------
Deployment with Docker
----------------------

It is also possible to deploy the server using Docker.

- First build the Docker Image.  From the server/scripts folder run::

  $ docker build .

- Find the id of the recently built image by running::

  $ docker images

- Use that id to launch the server::

  $ docker run -d -p 8888:80 <id>

This will run the docker container in the background, and translate calls from your host environment
port 8888 to the docker container port 80. At that point you should be able to access it like a normal website.


------------
Client tools
------------

Prerequisites:

* Python 2.7,
* Virtualenv (or another python sandboxing tool) is highly recommended.

General installation procedure:

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
you may want to check the `System specific install examples`_ section.


To run the server on this example dataset, follow the steps on
installing the server, then download and unpack the example data ::

    $ wget http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar.gz
    $ tar -zxvf ga4gh-example-data.tar.gz

An easier way to download and upack the data is to run the following
script, which will do these steps for you::

    $ python scripts/update_data.py

You can now run the server, which will by default serve variants from the sets in
the downloaded datafile::

    $ ga4gh_server

To change the data that is served,  a configuration file can be specified using
the ``-f <config_file>`` command line argument. Run::

    $ ga4gh_server --help

for details on the options for this program.


++++++++++++++++++++++
Example client queries
++++++++++++++++++++++

To run queries against this server, we can use the ``ga4gh_client`` program;
for example, here we run the ``variants/search`` method over the
``1000g_2013`` variant set, where the reference name is ``1``
and we only want calls returned for call set ID `1000g_2013.HG03279`::

    $ ga4gh_client variants-search http://localhost:8000/v0.5.1 -V 1000g_2013 -r 1 -c 1000g_2011.HG03279 | less -S

++++++++++++++++++++++++++++++++
System specific install examples
++++++++++++++++++++++++++++++++

MacOS X (with MacPorts)::

  $ cd [some working directory outside the ga4gh server directory tree]
  $ virtualenv --no-site-packages testenv
  $ source testenv/bin/activate
  $ cd [your ga4gh server directory]
  $ pip install -r requirements.txt
  $ python setup.py install
  $ nosetests

*TODO* Append examples of installs (using package managers if possible, no dependency
installs from source) on the target platform of your choice.

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
path to a configuration file (see the `Configuration file and data hierarchy`_
section for details).

++++++
Layout
++++++

The code for the project is held in the ``ga4gh`` package, which corresponds to
the ``ga4gh`` directory in the project root. Within this package, the
functionality is split between the ``client``, ``server``, ``protocol`` and
``cli`` modules.  The ``cli`` module contains the definitions for the
``ga4gh_client`` and ``ga4gh_server`` programs.
