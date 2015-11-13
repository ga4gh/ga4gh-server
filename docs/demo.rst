.. _demo:

**************
GA4GH API Demo
**************

In this demo, we'll install a copy of the GA4GH reference
implementation and run a local version of the server using some example
data. We then run some example queries on this server using various
different methods to illustrate the basics of the protocol.
The server can, of course, be run on any machine on the network,
but for simplicity we assume that the client and the server are running
on your local machine during this demo.

The instructions for installation here
are not intended to be used in a production deployment. See
the :ref:`installation` section for a detailed guide on production installation.
To run the demo, you will need a working installation of
`Python 2.7 <https://www.python.org/download/releases/2.7/>`_
and also have `virtualenv <https://virtualenv.pypa.io/en/latest/>`_
installed. We also need to have `zlib <http://www.zlib.net/>`_ and
a few other common libraries installed so that we can build some of the
packages that the reference server depends on.

On Debian/Ubuntu, for example, we can install these
packages using:

.. code-block:: bash

    $ sudo apt-get install python-dev python-virtualenv zlib1g-dev libxslt1-dev libffi-dev libssl-dev

On Fedora 22+ (current), the equivalent would be:

.. code-block:: bash

    $ sudo dnf install python-devel python-virtualenv zlib-devel libxslt-devel openssl-devel

First, we create a virtualenv sandbox to isolate the demo from the
rest of the system, and then activate it:

.. code-block:: bash

    $ virtualenv ga4gh-env
    $ source ga4gh-env/bin/activate

Now, install the `ga4gh package <https://pypi.python.org/pypi/ga4gh>`_
from the `Python package index <https://pypi.python.org/pypi>`_. This
will take some time, as some upstream packages will need to be built and
installed.

.. code-block:: bash

    (ga4gh-env) $ pip install ga4gh --pre

(Older versions of `pip <https://pip.pypa.io/en/latest/>`_ might not recognise
the ``--pre`` argument; if not, it is safe to remove it.)

Now we can download some example data, which we'll use for our demo:

.. code-block:: bash

    (ga4gh-env) $ wget http://www.well.ox.ac.uk/~jk/ga4gh-example-data-v3.0.tar
    (ga4gh-env) $ tar -xvf ga4gh-example-data-v3.0.tar

After extracting the data, we can then run the ``ga4gh_server`` application:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_server
    * Running on http://127.0.0.1:8000/ (Press CTRL+C to quit)
    * Restarting with stat

(The server is using a default configuration which assumes the
existence of the ``ga4gh-example-data`` directory for simplicity here; see
the :ref:`configuration` section for detailed information on how we configure the
server.) We now have a server running in the foreground. When it receives requests,
it will print out log entries to the terminal. A summary of the server's
configuration and data is available in HTML format at
``http://locahost:8000``, which can be viewed in a web browser.
Leave the server running and open another terminal to complete the
rest of the demo.

To try out the server, we must send some requests to it using the `GA4GH
protocol <http://ga4gh.org/#/api>`_. One way in which we can do this is to
manually create the `JSON <http://json.org/>`_ requests, and send these to the
server using `cURL <http://curl.haxx.se/>`_:

.. code-block:: bash

    $ curl --data '{}' --header 'Content-Type: application/json' \
    http://localhost:8000/datasets/search | jq .


In this example, we used the `searchDatasets
<http://ga4gh.org/documentation/api/v0.5.1/ga4gh_api.html#/schema/org.ga4gh.searchDatasets>`_
method to ask the server for all the Datasets on the server. It responded
by sending back some JSON, which we piped into the `jq <https://stedolan.github.io/jq/>`_
JSON processor to make it easier to read. We get the following result:

.. code-block:: json

    {
      "nextPageToken": null,
      "datasets": [
        {
          "description": null,
          "name": "1kg-p3-subset",
          "id": "MWtnLXAzLXN1YnNldA=="
        }
      ]
    }

In this example we sent a SearchDatasetsRequest object to the server
and received a SearchDatasetsResponse object in return. This response object
contained one Dataset object, which is contained in the ``datasets`` array.
This approach to interacting with the server is tedious and error prone, as
we have to hand-craft the request objects. It is also quite inconvenient, as
we may have to request many pages of objects to get all the objects
that satisfy our search criteria.

To simplify interacting with the server and to abstract away the low-level
network-level details of the server, we provide a client application.
To try this out, we start another instance of our virtualenv, and then send
the equivalent command using:

.. code-block:: bash

    $ source ga4gh-env/bin/activate
    (ga4gh-env) $ ga4gh_client datasets-search http://localhost:8000

::

    MWtnLXAzLXN1YnNldA==    1kg-p3-subset

The output of this command is a summary of the Datasets on that are present on the
server. We can also get the output in JSON form such that each
object is written on one line:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client datasets-search -O json http://localhost:8000

::

    {"description": null, "name": "1kg-p3-subset", "id": "MWtnLXAzLXN1YnNldA=="}

This format is quite useful for larger queries, and can be piped into jq
to extract fields of interest, pretty printing and so on.

We can perform similar queries for variant data using the
`searchVariants
<http://ga4gh.org/documentation/api/v0.5.1/ga4gh_api.html#/schema/org.ga4gh.searchVariants>`_
API call. First, we find the IDs of the VariantSets on the server using the
`searchVariantSets
<http://ga4gh.org/documentation/api/v0.5.1/ga4gh_api.html#/schema/org.ga4gh.searchVariantSets>`_
method:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client variantsets-search http://localhost:8000

::

    MWtnLXAzLXN1YnNldDptdm5jYWxs    mvncall

This tells us that we have one VariantSet on the server, with ID ``MWtnLXAzLXN1YnNldDptdm5jYWxs``
and name ``mvncall``. We can then search for variants overlapping a given interval in a VariantSet
as follows:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client variants-search http://localhost:8000 \
    --referenceName=1 --start=45000 --end=50000

The output of the client program is a summary of the data received in a
free text form. This is not intended to be used as the input to other
programs, and is simply a data exploration tool for users.
To really *use* our data, we should use a GA4GH client library.

Part of the GA4GH reference implementation is a :ref:`client-library`. This makes sending requests to the server and using the
responses very easy. For example, to run the same query as we
performed above, we can use the following code:

.. code-block:: python

    from __future__ import print_function

    import ga4gh.client as client

    httpClient = client.HttpClient("http://localhost:8000")
    # Get the datasets on the server.
    datasets = list(httpClient.searchDatasets())
    # Get the variantSets in the first dataset.
    variantSets = list(httpClient.searchVariantSets(
        datasetId=datasets[0].id))
    # Now get the variants in the interval [45000, 50000) on chromosome 1
    # in the first variantSet.
    iterator = httpClient.searchVariants(
        variantSetId=variantSets[0].id,
        referenceName="1", start=45000, end=50000)
    for variant in iterator:
        print(
            variant.referenceName, variant.start, variant.end,
            variant.referenceBases, variant.alternateBases, sep="\t")


If we save this script as ``ga4gh-demo.py`` we can then run it
using:

.. code-block:: bash

    (ga4gh-env) $ python ga4gh-demo.py

.. todo:: Add more examples of using the reads API and give
   examples of using the references API. We should aim to have
   a single complete example, where we start with a given
   variant, and drill down into the reads in question programatically.
   values as parameters which have sensible defaults.

---------
With OIDC
---------

.. todo:: Should we move the OIDC documentation into its own section?
    It is quite a lot of complication to add here to a beginners HOWTO.



If we want authentication, we must have an OIDC authentication provider.
One can be found in ``oidc-provider``, and run with the ``run.sh`` script.
We can then use this with the ``LocalOidConfig`` server configuration. So:

.. code-block:: bash

  $ cd oidc-provider && ./run.sh

In another shell on the same machine

.. code-block:: bash

  $ python server_dev.py -c LocalOidConfig

Make sure you know the hostname the server is running on. It can be found with

.. code-block:: bash

  $ python -c 'import socket; print socket.gethostname()'

With a web browser, go to ``https://<server hostname>:<server port>``. You may
need to accept the security warnings as there are probably self-signed
certificates. You will be taken through an authentication flow. When asked
for a username and password, try ``upper`` and ``crust``. You will find
yourself back at the ga4gh server homepage. On the homepage will be a
'session token' This is the key to access the server with the client tool
as follows:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client --key <key from homepage> variantsets-search https://localhost:8000/current
    MWtnLXAzLXN1YnNldDptdm5jYWxs    mvncall


