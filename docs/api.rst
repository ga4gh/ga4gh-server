.. _client-library:

*********************
Python Client Library
*********************

This is the GA4GH client API library. This is a convenient wrapper for
the low-level HTTP GA4GH API, and abstracts away network centric
details such as paging.  The methods and types used by the client
library are defined by the
`GA4GH schema <https://ga4gh-schemas.readthedocs.org/en/latest/>`_.


.. warning:: This client API should be considered early alpha quality,
          and may change in arbitrary ways. In particular, the current
          camelCase convention for identifiers is scheduled to change
          in the near future.

.. todo:: A full description of this API and links to a tutorial on how
       to use it, as well as a quickstart showing the basic usage.

-----
Types
-----

.. todo:: Add links to the upstream documentation for the GA4GH types.

----------
Client API
----------

.. todo:: Add overview documentation for the client API.

.. autoclass:: ga4gh.client.HttpClient
    :members: get_reference_set, get_reference,
        get_dataset, get_variant_set, get_variant,
        get_read_group_set, get_read_group,
        get_bio_sample, get_individual,
        search_datasets, search_reference_sets, search_references,
        search_variant_sets, search_variants, search_read_group_sets,
        search_reads, search_bio_samples, search_individuals

