.. _demo:

**************
GA4GH API Demo
**************

In this demo, we'll install a copy of the GA4GH reference
implementation and run a local copy of the server using some example
data. We then run some example queries on this server using various
different methods to illustrate the basics of the protocol.
The server can, of course, be run on any machine on network
but for simplicity, we assume that the client and the server are running
on your local machine during this demo.

The instructions for installation here
are not intended to be used in a production deployment. See [link to
deployment] for a detailed guide on production installation.
To run the demo, you will need a working installation of
`Python 2.7 <https://www.python.org/download/releases/2.7/>`_
and also have `virtualenv <https://virtualenv.pypa.io/en/latest/>`_
installed. We also need to have `zlib <http://www.zlib.net/>`
installed so that we can build some of the packages that the
reference server depends on.

On Debian/Ubuntu, for example, we can install these
packages using:

.. code-block:: bash

    $ sudo apt-get install python-dev python-virtualenv zlib1g-dev

First, we create a virtualenv sandbox to isolate the demo from the
rest of the system, and then activate it:

.. code-block:: bash

    $ virtualenv ga4gh-env
    $ source ga4gh-env/bin/activate

Now, install the `ga4gh package <https://pypi.python.org/pypi/ga4gh>`_
from the `Python package index <https://pypi.python.org/pypi>`_. This
will take some time, as some upstream packages will need to be built and
installed:

.. code-block:: bash

    (ga4gh-env) $ pip install ga4gh --pre

(Older versions of `pip <https://pip.pypa.io/en/latest/>`_ might not recognise
the ``--pre`` argument; if not, it is safe to remove it.)

Now we can download some example data, which we'll use for our demo:

.. code-block:: bash

    (ga4gh-env) $ wget http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar
    (ga4gh-env) $ tar -xvf ga4gh-example-data.tar

After extracting the data, we can then run the ``ga4gh_server`` application:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_server
    * Running on http://0.0.0.0:8000/ (Press CTRL+C to quit)
    * Restarting with stat

(The server is using a default configuration which assumes the
existence of the ``ga4gh-example-data`` for simplicity here; see
[configuration] for detailed information on how we configure the
server.) We now have a server running in the foreground. When it receives requests,
it will print out log entries to the terminal.
Now, leave the server running and open another terminal to complete the
rest of the demo.

To try out the server, we must send some requests to it using the `GA4GH
protocol <http://ga4gh.org/#/api>`_. One way in which we can do this is to
manually create the `JSON <http://json.org/>`_ requests, and send these to the
server using `cURL <http://curl.haxx.se/>`_:

.. code-block:: bash

    $ curl --data '{"datasetIds":[], "name":null}' --header 'Content-Type: application/json' \
    http://localhost:8000/v0.5.1/readgroupsets/search

In this example, we used the `searchReadGroupSets
<http://ga4gh.org/documentation/api/v0.5.1/ga4gh_api.html#/schema/org.ga4gh.searchReadGroupSets>`_
method to ask the server for all the ReadGroupSets on the server. It responded
by sending back some JSON, which cURL then printed to the terminal.

Creating these JSON requests by hand is tedious and error prone, and
so there is a client application to do this for us. To try this out, we
start another instance of our virtualenv, and then send the
equivalent command using:

.. code-block:: bash

    $ source ga4gh-env/bin/activate
    (ga4gh-env) $ ga4gh_client readgroupsets-search http://localhost:8000/v0.5.1

The output of this command is a simple summary of the ReadGroupSets that
are present on the server. We can also see the JSON messages passing
between the client and the server if we increase the verbosity level:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client -vv readgroupsets-search http://localhost:8000/v0.5.1

We can perform similar queries for variant data using the
`searchVariants
<http://ga4gh.org/documentation/api/v0.5.1/ga4gh_api.html#/schema/org.ga4gh.searchVariants>`_
API call. First, we find the IDs of the VariantSets on the server using the
`searchVariantSets
<http://ga4gh.org/documentation/api/v0.5.1/ga4gh_api.html#/schema/org.ga4gh.searchVariantSets>`_
method:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client variantsets-search http://localhost:8000/v0.5.1
    1kg-phase1
    1kg-phase3

This tells us that we have two VariantSets on the server, with IDs ``1kg-phase1``
and ``1kg-phase3``. In our example data, these correspond to a subset of the
data from `1000 Genomes <http://www.1000genomes.org/>`_ phases 1 and 3.

We can then search for variants overlapping a given interval in a VariantSet
as follows:

.. code-block:: bash

    (ga4gh-env) $ ga4gh_client variants-search http://localhost:8000/v0.5.1 \
    --variantSetIds=1kg-phase1 --referenceName=2 --start=33100 --end=34000

The output of the client program is a summary of the data received in a
free text form. This is not intended to be used as the input to other
programs, and is simply a data exploration tool for users.
To really *use* our data, we should use a GA4GH client library.

Part of the GA4GH reference implementation is a Python client-side
library. This makes sending requests to the server and using the
responses very easy. For example, to run the same query as we
performed above, we can use the following code:

.. code-block:: python

    from __future__ import print_function

    import ga4gh.client as client
    import ga4gh.protocol as protocol

    httpClient = client.HttpClient("http://localhost:8000/v0.5.1")
    request = protocol.GASearchVariantsRequest()
    request.variantSetIds = ["1kg-phase1"]
    request.referenceName = "2"
    request.start = 33100
    request.end = 34000
    for variant in httpClient.searchVariants(request):
        print(
            variant.referenceName, variant.start, variant.end,
            variant.referenceBases, variant.alternateBases, sep="\t")


If we save this script as ``ga4gh-demo.py`` we can then run it
using:

.. code-block:: bash

    (ga4gh-env) $ python ga4gh-demo.py


**TODO**

1. Add more examples of using the reads API and give
   examples of using the references API. We should aim to have
   a single complete example, where we start with a given
   variant, and drill down into the reads in question programatically.
2. Update the client API to be more user-friendly. We shouldn't need
   to create an instance of ``GASearchVariantsRequest`` to call
   ``searchVariants``. Rather, ``searchVariants`` should have the corresponding
   values as parameters which have sensible defaults.
