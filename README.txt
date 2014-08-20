==============================
GA4GH Reference Implementation
==============================

A reference implementation of the APIs defined in the schemas repository.

*************************
Initial skeleton overview
*************************

This is a proposed skeleton layout for the GA4GH reference client and
server applications. As such, nothing is finalised and all aspects of
the design and implementation are open for discussion and debate. The overall
goals of the project are:

Simplicity/clarity
    The main goal of this implementation is to provide an easy to understand
    and maintain implementation of the GA4GH API. Design choices
    are driven by the goal of making the code as easy to understand as
    possible, with performance being of secondary importance. With that
    being said, it should be possible to provide a functional implementation
    that is useful in many cases where the extremes of scale are not
    important.

Portability
    The code is written in Python for maximum portability, and it
    should be possible to run on any modern computer/operating system (Windows
    compatibility should be possible, although this has not been tested).  We use a
    subset of Python 3 which is backwards compatible with Python 2 following the
    current `best practices <http://python-future.org/compatible_idioms.html>`_.
    In this way, we fully support both Python 2 and 3.

Ease of use
    The code follows the `Python Packaging User Guide
    <http://python-packaging-user-guide.readthedocs.org/en/latest/>`_. This will
    make installing the ``ga4gh`` reference code very easy across a range of
    operating systems.

Currently, the project includes a minimal working version of the ``variants/search``
method using simulated data, and simple client and server applications to
exercise this API.

*******************
Simulating Variants
*******************

To provide a working outline of the project as quickly as possible, we took the
approach of serving simulated data from the ``variants/search`` method.  We
currently use a very naive and simplistic model for generating variants, in
which SNPs are distributed over all query regions uniformly at a fixed density.
The model is simplistic, but it is  *consistent*: regardless of
the context, any query that intersects with a given position will always return
precisely the same variant. Thus, from a client's perspective, it would appear
that  the server is returning (not very plausible) variants drawn from an
infinitely large file.

Only SNPs are returned presently, but more complex variants could also be
added. For example, large indels could be included by generating a list of
variants with the required properties when the server starts.  We then build a
lookup table of the extremities of these regions. When we are generating a
variant for a particular site, we first search to see if it intersects with any
of the generated indels; if it does, we return the indel. If not, we generate a
SNP as before. This can be implemented quite efficiently if the number of
large variants is reasonably small and they don't overlap.

The implementation will certainly need to serve from existing data
files, but there are distinct advantages in also being able to simulate
data:

- We can simulate arbitrarily large datasets without worrying about
  memory usage and needing to deal with difficult problems of scale.
  This is useful for the initial development of the protocol,
  allowing us to (for example) quickly benchmark protocol throughput.
  This  would also be a useful tool for developers writing client code in the
  longer term.
- Devising a reasonable model for simulating data forces us to be explicit
  about the assumptions we make about the underlying data. If we just
  work on translating data from existing formats, there is the danger that
  we inherit their assumptions.
- Working exclusively from existing data files also biases the development
  towards the types of data that are currently easily available. To consider
  the full implications of choices made in the protocol design, we should
  be able to serve data under a wide variety of assumptions (e.g.
  different ploidy levels, etc.).
- The availability of a good source of random data is essential for large
  scale testing.

Generating variants randomly in a consistent manner is relatively
straightforward; however, implementing the reads API in a similar
fashion would be much more difficult. Thus, while the arguments above
in favour of implementing a simulated backend for reads still hold,
they may not be sufficiently strong to justify the effort required.


*************
Trying it out
*************

The project is designed to be published as a `PyPI <https://pypi.python.org/pypi>`_
package, so ultimately installing the reference client and server programs
should be as easy as::

    $ pip install ga4gh

However, the code is currently only a proposal, so it has not been uploaded to
the Python package index. The best way to try out the code right now is to
use `virtualenv <http://virtualenv.readthedocs.org/en/latest/>`_. After cloning
the git repo, and changing to the project directory, do the following::

    $ virtualenv testenv
    $ source testenv/bin/activate
    $ python setup.py install

This should install the ``ga4gh_server`` and ``ga4gh_client`` scripts into the
virtualenv and update your ``PATH`` so that they are available. When you have
finished trying out the programs you can leave the virtualenv using::

    $ deactivate

The virtualenv can be restarted at any time, and can also be deleted
when you no longer need it.

The easiest way to run the client and server programs is to keep both
on your local machine, and use one shell for each. For example, we
can start the server running in one shell using::

    $ ga4gh_server simulate

We can then send messages to this server using, e.g.::

    $ ga4gh_client variants-search --start=400 --end=405

This tells the client to sent a ``variants/search`` message in which we
set the ``start`` attribute to ``400`` and the ``end`` attribute to ``405``.
The results are printed out in a crude VCF-like manner. We can see more
detail about the protocol messages being exchanged by adding ``-v`` options
to turn up the verbosity.

To get help on the various options for the client and server programs
use the ``help`` command or the ``-h`` option to get help for a specific
command. For example::

    $ ga4gh_client variants-search -h
    usage: ga4gh_client variants-search [-h] [--referenceName REFERENCENAME]
                                        [--variantName VARIANTNAME]
                                        [--start START] [--end END]
                                        [--maxResults MAXRESULTS]

    Search for variants

    optional arguments:
      -h, --help            show this help message and exit
      --referenceName REFERENCENAME, -r REFERENCENAME
                            Only return variants on this reference.
      --variantName VARIANTNAME, -n VARIANTNAME
                            Only return variants which have exactly this name.
                            TODO
      --start START, -s START
                            The start of the search range (inclusive).
      --end END, -e END     The end of the search range (exclusive).
      --maxResults MAXRESULTS, -m MAXRESULTS
                            The maximum number of variants returned in one
                            response.

******
Layout
******

The code for the project is held in the ``ga4gh`` package, which corresponds
to the ``ga4gh`` directory in the project root. Within this package,
the functionality is split between the ``client``, ``server`` and
``protocol`` modules. There is also a subpackage called ``scripts``
which holds the code defining the command line interfaces for the
``ga4gh_client`` and ``ga4gh_server`` programs.

For development purposes, it is useful to be able to run the command
line programs directly without installing them. To do this, make hard links
to the files in the scripts directory to the project root and run them
from there; e.g::

    $ ln ga4gh/scripts/server.py .
    $ python server.py
    usage: server.py [-h] [--port PORT] [--verbose] {help,simulate} ...
    server.py: error: too few arguments

++++++++++++
Coding style
++++++++++++

The code follows the guidelines of `PEP 8
<http://legacy.python.org/dev/peps/pep-0008>`_ in most cases. The only notable
difference is the use of camel case over underscore delimited identifiers; this
is done for consistency with the GA4GH API. The code was checked for compliance
using the `pep8 <https://pypi.python.org/pypi/pep8>`_ tool.

****
TODO
****

Here is a partial list of the outstanding issues with the current
implementation and discussion of what might be done to address these
issues.

++++++++++++++++++++++++++++++++++++
Completing the variant/search method
++++++++++++++++++++++++++++++++++++

The variants/search method currently only handles the very simplest case
in which we search for variants between ``start`` and ``end``.

++++++++++++++
Error handling
++++++++++++++

No attempt is currently made to handle errors within the client or
server. This can be extended fairly easily, but there are some issues
that need to be addressed by the protocol. For example, are there standard
errors to be raised when a non-conforming protocol message is received by
the server? Does the protocol specify some classes of error or leave this
entirely up to the implementation? When do we return HTTP error status codes
and when do we return a HTTP success along with a ``GAException`` object?

We also need to add functionality to detect protocol errors in the
``ga4gh.protocol`` module. These fall into several classes:

1. Missing or additional attributes in the JSON definition of a
   ``ProtocolElement``;
2. Missing mandatory values, or mutually contradictory values;
3. Type errors;
4. Bounds errors (e.g. maxResults=-1);
5. More?

+++++++++++++++++++++++++++
Adding all Protocol classes
+++++++++++++++++++++++++++

Presently, only a small subset of the GA4GH API has been implemented in the
``ga4gh.protocol`` module. These classes are simple copies of the classes
defined in Avro, and inherit from a superclass to provide the JSON
serialisation/deserialisation functionality. The classes contain an instance
variable for all of the attributes defined in JSON, plus a simple annotation
to help define the acceptable types.

It is a simple (if tedious) job to convert the Avro definitions into the
corresponding Python classes. However, this is not a very satisfactory
approach since the definitions and documentation of the classes will
inevitably become out of sync with the Avro definitions. It would be
much better if we could devise a way to automatically generate the Python
code from the Avro definitions. This could then be regenerated and checked
into git each time the definitions go through a version change.

++++++++++
Test cases
++++++++++

There are currently no test cases. A comprehensive test suite must be included
with a reference implementation such as this. There would be considerable
overlap between the process of rigorously testing this implementation
and in defining an API compliance test suite. Therefore, one approach might
be to add a ``ga4gh.compliance`` module which would provide the definitions
of what a compliant API must be. This could be utilised by the internal
test suite, as well as forming the basis for a compliance checking
script.

++++++++++++++++++++++++++++++++
Adding support for the reads API
++++++++++++++++++++++++++++++++

Only the variants API is supported at the moment, so the reads API must be
added.

+++++++++++++++++++++++++++++++++++++++++++++
Adding support for serving from VCF/SAM files
+++++++++++++++++++++++++++++++++++++++++++++

Using a simulator to generate variants and reads is useful for development
purposes, but it is essential that we support real data files also. There
are some technical challenges in this, if we wish to support realistic
file  sizes.


