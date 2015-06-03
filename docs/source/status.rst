.. _status:

---------------------
Implementation Status
---------------------

The GA4GH server is currently under active development, and several
features are not yet fully functional.  Here are some things to keep in
mind about the current server's operation:

* There is good supports for the variants API.  This should be production
  ready. 
* We have initial support for the reads API.  This works, but needs
  expert review (especially for verifying the conversion between BAM
  data fields and GA4GH objects).
* Datasets are currently not implemented.  As such all `datasetId` fields
  in requests are ignored.
* References are currently not supported.

For more detail on individual development issues, please see the project's
`issue page <https://github.com/ga4gh/server/issues>`_.
