.. _status:

------
Status
------

The GA4GH server is currently under active development, and several
features are not yet fully functional.  These mostly involve the
reads API. Some missing features are:

- Unmapped reads. We do not support searching for unmapped reads.

- Searching over multiple ReadGroups.

For more detail on individual development issues, please see the project's
`issue page <https://github.com/ga4gh/server/issues>`_.

+++++++++++++
Release Notes
+++++++++++++

*****
0.2.1
*****

Bug release to fix broken upstream packages.

*****
0.2.0
*****

No changes from 0.2.0a1


*******
0.2.0a1
*******

Alpha pre-release supporting major schema update. This release is backwards
incompatible with previous releases, and requires a revised data directory
layout.

- Schema version changed from v0.5 to v0.6.0a1

- Various backwards incompatible changes to the data directory layout

- Almost complete support for the API.

- Numerous code layout changes.

*****
0.1.2
*****

This bugfix release addresses the following issues:

- #455: bugs in reads/search call (pysam calls not sanitized, wrong
  number of arguments to getReadAlignments)

- #433: bugs in paging code for reads and variants

*****
0.1.1
*****

- Fixes dense variants not being correctly handled by the server (#391)

- Removes unused paths (thus they won't confusingly show up in the HTML
  display at / )

*****
0.1.0
*****

Just bumping the version number to 0.1.0.

*******
0.1.0b1
*******

This is a beta pre-release of the GA4GH reference implementation. It includes

- A fully functional client API;

- A set of client side tools for interacting with any conformant server;

- A partial implementation of the server API, providing supports for variants and
  reads from native  file formats.


*******
0.1.0a2
*******

This is an early alpha release to allow us to test the PyPI package and
the README. This is not intended for general use.
