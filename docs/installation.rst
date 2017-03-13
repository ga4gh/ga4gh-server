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

First, we install some basic pre-requisite packages:

.. code-block:: bash

  sudo apt-get install python-dev python-virtualenv zlib1g-dev libxslt1-dev libffi-dev libssl-dev libcurl4-openssl-dev

Install Apache and mod_wsgi, and enable mod_wsgi:

.. code-block:: bash

  sudo apt-get install apache2 libapache2-mod-wsgi
  sudo a2enmod wsgi

Create the Python egg cache directory, and make it writable by
the ``www-data`` user:

.. code-block:: bash

  sudo mkdir /var/cache/apache2/python-egg-cache
  sudo chown www-data:www-data /var/cache/apache2/python-egg-cache/

Create a directory to hold the GA4GH server code, configuration
and data. For convenience, we make this owned by the current user
(but make sure all the files are world-readable).:

.. code-block:: bash

  sudo mkdir /srv/ga4gh
  sudo chown $USER /srv/ga4gh
  cd /srv/ga4gh

Make a virtualenv, and install the ga4gh package:

.. code-block:: bash

  virtualenv ga4gh-server-env
  source ga4gh-server-env/bin/activate
  pip install ga4gh-server
  deactivate

Download and unpack the example data:

.. code-block:: bash

  wget https://github.com/ga4gh/ga4gh-server/releases/download/data/ga4gh-example-data_4.6.tar
  tar -xf ga4gh-example-data_4.6.tar

Create the WSGI file at ``/srv/ga4gh/application.wsgi`` and write the following
contents:

.. code-block:: python

    from ga4gh.server.frontend import app as application
    import ga4gh.server.frontend as frontend
    frontend.configure("/srv/ga4gh/config.py")

Create the configuration file at ``/srv/ga4gh/config.py``, and write the
following contents:

.. code-block:: python

    DATA_SOURCE = "/srv/ga4gh/ga4gh-example-data/registry.db"

Note that it is expected that the user running the server, `www-data`, 
have write and read access to the directories containing data files.

(Many more configuration options are available --- see the :ref:`configuration`
section for a detailed discussion on the server configuration and input data.)

Configure Apache. Note that these instructions are for Apache 2.4 or greater.
Edit the file ``/etc/apache2/sites-available/000-default.conf``
and insert the following contents towards the end of the file
(*within* the ``<VirtualHost:80>...</VirtualHost>`` block):

.. code-block:: apacheconf

    WSGIDaemonProcess ga4gh \
        processes=10 threads=1 \
        python-path=/srv/ga4gh/ga4gh-server-env/lib/python2.7/site-packages \
        python-eggs=/var/cache/apache2/python-egg-cache
    WSGIScriptAlias /ga4gh /srv/ga4gh/application.wsgi

    <Directory /srv/ga4gh>
        WSGIProcessGroup ga4gh
        WSGIApplicationGroup %{GLOBAL}
        Require all granted
    </Directory>

.. warning::

    Be sure to keep the number of threads limited to 1 in the WSGIDaemonProcess
    setting. Performance tuning should be done using the processes setting.

The instructions for configuring Apache 2.2 (on Ubuntu 14.04) are the same as
above with thee following exceptions:

You need to edit
``/etc/apache2/sites-enabled/000-default``

instead of
``/etc/apache2/sites-enabled/000-default.conf``

And while in that file, you need to set permissions for the directory to

.. code-block:: apacheconf

    Allow from all

instead of

.. code-block:: apacheconf

    Require all granted



Now restart Apache:

.. code-block:: bash

  sudo service apache2 restart

We will now test to see the server started properly by requesting the
landing page.

.. code-block:: bash

    curl http://localhost/ga4gh/ --silent | grep GA4GH
    #         <title>GA4GH reference server 0.2.3.dev4+nge0b07f3</title>
    #    <h2>GA4GH reference server 0.2.3.dev4+nge0b07f3</h2>
    # Welcome to the GA4GH reference server landing page! This page describes

We can also test the server by running some API commands. Please refer to
the instructions in the :ref:`demo` for how to access data made available
by this server.

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

+++++++++++++++
Troubleshooting
+++++++++++++++

Server errors will be output to the web server's error log by default (in Apache on
Debian/Ubuntu, for example, this is ``/var/log/apache2/error.log``). Each client
request will be logged to the web server's access log (in Apache on Debian/Ubuntu
this is ``/var/log/apache2/access.log``).

For more server configuration options see :ref:`Configuration`

--------------------
Deployment on Docker
--------------------
It is also possible to deploy the server using Docker.

First, you need an environment running the docker daemon. For non-production use, we recommend `boot2docker <http://boot2docker.io/>`_. For production use you should install docker on a stable linux distro.
Please reference the `platform specific Docker installation instructions <https://docs.docker.com/installation/>`_. OSX and Windows are instructions for boot2docker.

**Local Dataset Mounted as Volume**

If you already have a dataset on your machine, you can download and deploy the apache server in one command:

.. code-block:: bash

  docker run -e GA4GH_DATA_SOURCE=/data -v /my/ga4gh_data/:/data:ro -d -p 8000:80 --name ga4gh_server ga4gh/ga4gh-server:latest

Replace ``/my/ga4gh_data/`` with the path to your data.

This will:

* pull the automatically built image from `Dockerhub <https://registry.hub.docker.com/u/ga4gh/ga4gh-server/>`_
* start an apache server running mod_wsgi on container port 80
* mount your data read-only to the docker container
* assign a name to the container
* forward port 8000 to the container.

For more information on docker run options, see the `run reference <https://docs.docker.com/reference/run/>`_.

**Demo Dataset Inside Container**

If you do not have a dataset yet, you can deploy a container which includes the demo data:

.. code-block:: bash

  docker run -d -p 8000:80 --name ga4gh_demo ga4gh/ga4gh-server:latest

This is identical to the production container, except that a copy of the demo data is included and appropriate defaults are set.

**Developing Client Code: Run a Client Container and a Server**

In this example you run a server as a daemon in one container, and the client as an ephemeral instance in another container.
From the client, the server is accessible at ``http://server/``, and the ``/tmp/mydev`` directory is mounted at ``/app/mydev/``. Any changes you make to scripts in ``mydev`` will be reflected on the host and container and persist even after the container dies.

.. code-block:: bash

  # make a development dir and place the example client script in it
  mkdir /tmp/mydev
  curl https://raw.githubusercontent.com/ga4gh/ga4gh-server/master/scripts/demo_example.py > /tmp/mydev/demo_example.py
  chmod +x /tmp/mydev/demo_example.py

  # start the server daemon
  # assumes the demo data on host at /my/ga4gh_data
  docker run -e GA4GH_DEBUG=True -e GA4GH_DATA_SOURCE=/data -v /my/ga4gh_data/:/data:ro -d --name ga4gh_server ga4gh/ga4gh-server:latest

  # start the client and drop into a bash shell, with mydev/ mounted read/write
  # --link adds a host entry for server, and --rm destroys the container when you exit
  docker run -e GA4GH_DEBUG=True -v /tmp/mydev/:/app/mydev:rw -it --name ga4gh_client --link ga4gh_server:server --entrypoint=/bin/bash --rm ga4gh/ga4gh-server:latest

  # call the client code script
  root@md5:/app# ./mydev/demo_example.py

  # call the command line client
  root@md5:/app# ga4gh_client variantsets-search http://server/current

  #exit and destroy the client container
  root@md5:/app# exit

**Ports**

The ``-p 8000:80`` argument to ``docker run`` will run the docker container in the background, and translate calls from your host environment
port 8000 to the docker container port 80. At that point you should be able to access it like a normal website, albeit on port 8000.
Running in `boot2docker <http://boot2docker.io/>`_, you will need to forward the port from the boot2docker VM to the host.
From a terminal on the host to forward traffic from localhost:8000 to the VM 8000 on OSX:

.. code-block:: bash

  VBoxManage controlvm boot2docker-vm natpf1 "ga4gh,tcp,127.0.0.1,8000,,8000"

For more info on port forwarding see `the VirtualBox manual <https://www.virtualbox.org/manual/ch06.html#natforward>`_ and this `wiki article <https://github.com/CenturyLinkLabs/panamax-ui/wiki/How-To%3A-Port-Forwarding-on-VirtualBox>`_.

++++++++
Advanced
++++++++

If you want to build the images yourself, that is possible. The `ga4gh/ga4gh-server repo <https://registry.hub.docker.com/u/ga4gh/ga4gh-server/>`_
builds automatically on new commits, so this is only needed if you want to modify the Dockerfiles, or build from a different source.

The prod and demo builds are based off of `mod_wsgi-docker <https://github.com/GrahamDumpleton/mod_wsgi-docker>`_, a project from the author of mod_wsgi.
Please reference the Dockerfiles and documentation for that project during development on these builds.

**Examples**

Build the code at server/ and run for production, serving a dataset on local host located at ``/my/dataset``

.. code-block:: bash

 cd server/
 docker build -t my-repo/my-image .
 docker run -e GA4GH_DATA_SOURCE=/dataset -v /my/dataset:/dataset:ro -itd -p 8000:80 --name ga4gh_server my-repo/my-image

Build and run the production build from above, with the demo dataset in the container
(you will need to modify the FROM line in ``/deploy/variants/demo/Dockerfile`` if you want to use your image from above as the base):


++++++++++++++++++++++
Troubleshooting Docker
++++++++++++++++++++++

**DNS**

The docker daemon's DNS may be corrupted if you switch networks, especially if run in a VM.
For boot2docker, running udhcpc on the VM usually fixes it.
From a terminal on the host:

.. code-block:: bash

  eval "$(boot2docker shellinit)"
  boot2docker ssh
  >	sudo udhcpc
  (password is tcuser)

**DEBUG**

To enable DEBUG on your docker server, call docker run with ``-e GA4GH_DEBUG=True``

.. code-block:: bash

  docker run -itd -p 8000:80 --name ga4gh_demo -e GA4GH_DEBUG=True ga4gh/ga4gh-server:latest

This will set the environment variable which is read by config.py

You can then get logs from the docker container by running ``docker logs (container)`` e.g. ``docker logs ga4gh_demo``

.. _macosinstall:

----------------------------------------------
Installing the development version on Mac OS X
----------------------------------------------

**Prerequisites**

First install libraries and header code for
`Python 2.7 <https://www.python.org/download/releases/2.7/>`_.
It will be a lot easier if you have `Homebrew <http://brew.sh/index.html>`_,
the "missing package manager" for OS X, installed first.
To install Homebrew, paste the following at a Terminal prompt ($):

.. code-block:: bash

  /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Now use ``brew install`` to install Python if you don't have Python 2.7
installed and then ``pip install``, which comes with Python, can be used to
install virtual environment:

.. code-block:: bash

  brew install python
  pip install virtualenv

**Install**

Download source code from GitHub to the project target folder, here assumed to 
be ``~/ga4gh``: (If you haven't already done so, 
`set up github <https://help.github.com/articles/set-up-git/>`_ 
to work from your command line.)

.. code-block:: bash

  git clone https://github.com/ga4gh/ga4gh-server.git

Before installing Python library dependencies, create a virtualenv sandbox to 
isolate it from the rest of the system, and then activate it:

.. code-block:: bash

  cd server
  virtualenv ga4gh-env
  source ga4gh-env/bin/activate

Install Python dependencies:

.. code-block:: bash

  pip install -r dev-requirements.txt -c constraints.txt

**Test and run**

Run tests to verify the install:

.. code-block:: bash

  ga4gh_run_tests

Please refer to the instructions in the :ref:`demo` for how to access
data made available by this server.

