.. _installation:

************
Installation
************

This section documents the process of deploying the GA4GH reference
server in a production setting. The intended audience is therefore
server administrators. If you are looking for a quick demo of the
GA4GH API using a local installation of the reference server
please check out the :ref:`demo`. If you are looking for
instructions to get a development system up and running, then
please go to the :ref:`development` section.

--------------------
Deployment on Apache
--------------------

To deploy on Apache on Debian/Ubuntu platforms, do the following.

- Install some basic pre-requisite packages:


.. code-block:: bash

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


These instructions are just one way in which we can achieve the same thing.
There are any number of different ways in which we can set up a WSGI
application under Apache, which may be preferable in different installations.
(In particular, the Apache configuration here may be specific to
Ubuntu 14.04, where this was tested.)
See the `mod_wsgi documentation <https://code.google.com/p/modwsgi/>`_ for
more details. These instructions are also specific to Debian/Ubuntu and
different commands and directory structures will be required on
different platforms.

The server can be deployed on any WSGI compliant web server. See the
instructions in the `Flask documentation
<http://flask.pocoo.org/docs/0.10/deploying/>`_ for more details on
how to deploy on various other servers.

**TODO**

1. Add more detail on how we can test out the API by making some client
   queries.
2. Add links to the Configuration section to give details on how we
   configure the server.
3. Change the example blocks above to use explicit code blocks with
   the correct syntax highlighting.

+++++++++++++++
Troubleshooting
+++++++++++++++

If you are encountering difficulties getting the above to work, it is helpful
to turn on debugging output. Do this by adding the following line to your
config file::

    DEBUG = True

When an error occurs, the details of this will then be printed to the web server's
error log (in Apache on Debian/Ubuntu, for example, this is ``/var/log/apache2/error.log``).

--------------------
Deployment on Docker
--------------------

**TODO**
