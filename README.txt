==============================
GA4GH Reference Implementation
==============================

A reference implementation of the APIs defined in the schemas repository.

*************************
Initial skeleton overview
*************************

Goals:

Simplicity/clarity
    The main goal of this implementation is to provide an easy to understand
    and maintain implementation of the GA4GH API. As such, design choices 
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
    In this way, we fully support both Python 2 and 3, giving excellent 
    portability.

Ease of use
    The code is packaged as a `Python package index <https://pypi.python.org/pypi>`_ 
    package following the `Python Packaging User Guide
    <http://python-packaging-user-guide.readthedocs.org/en/latest/>`_. This will
    make installing the ``ga4gh`` reference code very easy across a range of 
    operating systems.

*************
Trying it out
*************

The project is designed to be published as a `PyPI <https://pypi.python.org/pypi>`_
package, so ultimately installing the reference client and server programs 
(and API) should be as easy as::
    
    pip install ga4gh

However, the code is currently only a proposal, so it has not been uploaded to
the Python package index. The best way to try out the code right now is to 
use `virtualenv <http://virtualenv.readthedocs.org/en/latest/>`_. After cloning
the git repo, and changing to the project directory, do the following::
    
    virtualenv testenv
    source testenv/bin/activate
    python setup.py install

This should install the ``ga4gh_server`` and ``ga4gh_client`` scripts into 
your ``PATH``. Once you are done trying out the programs you can leave
the virtualenv using::

    deactivate

The virtualenv can be restarted at any time, and can also be deleted
when you no longer need it.

The easiest way to run the client and server programs is to keep both 
on your local machine, and use one shell for each. For example, we 
can start the server running in one shell using::

    ga4gh_server simulate

We can then send messages to this server using, e.g.::

    ga4gh_client variants-search --start=400 --end=405

This tells the client to sent a ``variants/search`` message in which we 
set the ``start`` attribute to ``400`` and the ``end`` attribute to ``405``.
The results are printed out in a crude VCF-like manner. We can see more 
detail about the protocol messages being exchanged by adding ``-v`` options
to turn up the verbosity.

******
Layout
******

The 

****
TODO
****

Here is a partial list of the outstanding issues with the current 
implementation and discussion of what might be done to address these
issues.

++++++++++++++++++++++++++++++++++++
Completing the variant/search method
++++++++++++++++++++++++++++++++++++

The variants/search method currently only handles the very simplest case of
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


