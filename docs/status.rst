.. _status:

------
Status
------

The GA4GH server is currently under active development, and several
features are not yet fully functional.  These mostly involve the
reads API. Some missing features are:

- Unmapped reads. We do not support searching for unmapped reads.

- Searching over multiple ReadGroups in different ReadGroupSets.

For more detail on individual development issues, please see the project's
`issue page <https://github.com/ga4gh/server/issues>`_.

+++++++++++++
Release Notes
+++++++++++++

*****
0.3.5
*****

Alpha pre-release supporting minor update. We have done some restructuring,
some bug fixes and some minor protocol updates in this release.

- Restructuring: We have reorganized our codebase to separate out the schemas,
   client, and server components as separately installable Python components.
   Running 'pip install ga4gh --pre' will still automatically install all
   required packages. But now developers who wish to leverage the Python
   libraries around the schemas can just use 'pip install ga4gh-schemas',
   and developers of client code can gain access to the client API library
   my just running 'pip install ga4gh-client'.
- Bug fixes and improvements
    - Fix to be able to handle VCFs with genotype == './.' (server #1389)
    - Fixed bug with OIDC Authentication (server #1452, #1508)
    - Repo manager list out of range error (server #1472)
    - G2P features search was returning a 404, this has been fixed 
      (server #1379)
    - Sped up readGroupSet generator (server #1316)
    - Enhanced ENCODE RNA downloader
    - Added support for Auth0 (server #1452)
    - Refactored usage of the term 'biosample' to be consistent 
      (server #1394)
    - Moved to Protobuf version 3.1 (server #1471)
- Documentation updates
    - Github usage
    - Added sections for new Python packages: client and schemas
    - Split out repository manager docs
    - Docker file documentation updated
    - Updated the Configuration section to document the Auth0 settings

*****
0.3.4
*****

Alpha pre-release supporting major feature update.

- G2P functionality added to support the following API endpoints:
   - POST `/phenotypeassociationsets/search`
   - POST `/phenotypes/search`
   - POST `/featurephenotypeassociations/search`
- Biometadata tags for RNA quantifications.
- Improvements to the RNA quantification ingestion pipeline.
- Migrated CLI related code to `cli` module.
- Add demonstration RNA quantification data.
- Minor doc fixes

Known Issues

- When searching using a wildcard, `*`, an Internal Server Error 
  occurs. #1379
- When listing many Read Group Sets, responses can be quite slow
  causing timeouts. #1316


*****
0.3.3
*****

Alpha pre-release supporting major feature update.

- RNA functionality added to support the following API endpoints:
   - POST /rnaquantificationsets/search
   - GET /rnaquantificationsets/{id}
   - POST /rnaquantifications/search
   - GET /rnaquantifications/{id}
   - POST /expressionlevels/search
   - GET /expressionlevels/{id}

- Fixed bug where transcript effects would be repeated within a 
  search result.


*****
0.3.2
*****

Alpha pre-release supporting major feature update.

- Metadata functionality added to support biosample and individual metadata
  capabilities.

- Now support searching features by 'name' and 'gene_symbol'. These fields
  have been promoted to facilitate the future integration of the RNA and
  G2P modules.


*****
0.3.1
*****

Alpha pre-release supporting major feature update. This release is not
backwards compatible with previous releases due to several changes to 
the schemas that were required to move to protocol buffers.

- This release includes the code changes necessary for the migration 
  to protocol buffers from Avro.

- Client applications will need to be rebuilt to the new schemas and 
  use the protobuf json serialization libraries to be compatible 
  with this version of the server. 


*****
0.3.0
*****

Alpha pre-release supporting major feature update. This release is not
backwards compatible with previous releases, and requires the data files
be re-imported.

- File locations are now tracked in a registry.db registry such that the
  files can be located anywhere. The information from the json sidecar
  files are also included in the database.

- Ontology terms are now imported via an OBO file instead of the old
  pre-packaged sequence_ontology.txt file. A sample OBO file has been
  added to the sample data set for the reference server.

- Added a custom landing page option for the main page of the server.

- Performance improvement for variant search when calls are set to an empty
  string.

- Improved server configuration including Apache configuration and
  robots.txt file.

*****
0.2.2
*****

Alpha pre-release supporting major feature update. This release is backwards
incompatible with previous releases, and requires a revised data directory
layout.

- Added sequence and variant annotations (which introduces a sqlite
  database component)

- Added repo manager, a command line tool to manage data files and
  import them into the server's data repository

- Supported searching over multiple ReadGroups, so long as they are
  all in the same ReadGroupSet and all of the ReadGroups in the
  ReadGroupSet are specified

*****
0.2.1
*****

Bugfix release that fixes a problem introduced by upstream package changes

*****
0.2.0
*****

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
