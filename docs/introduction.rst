.. _introduction:

************
Introduction
************

The `Data Working Group <http://ga4gh.org/#/>`_ of the
`Global Alliance for Genomics and Health <http://genomicsandhealth.org/>`_
has defined an
`API <https://ga4gh-schemas.readthedocs.org/en/latest/>`_
to facilitate interoperable exchange of genomic data.
This is the the documentation for the reference implementation of the API.

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
    for details). This provides easy installation of the ``ga4gh`` reference code
    across a range of operating systems.
