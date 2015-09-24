.. _client-library:

*********************
Python Client Library
*********************

This is the GA4GH client API library. This is a convenient wrapper for
the low-level HTTP GA4GH API, and abstracts away network centric
details such as paging.


.. warning:: This client API should be considered early alpha quality,
          and may change in arbitrary ways. In particular, the current
          camelCase convention for identifiers may change to snake_case
          in future.

.. todo:: A full description of this API and links to a tutorial on how
       to use it, as well as a quickstart showing the basic usage.

-----
Types
-----

.. todo:: A short overview of the types and links to the high-level
        docs.

++++++++++
References
++++++++++

.. autoclass:: ga4gh.protocol.ReferenceSet
    :members:

.. autoclass:: ga4gh.protocol.Reference
    :members:


++++++++
Datasets
++++++++

.. autoclass:: ga4gh.protocol.Dataset
    :members:

++++++++
Variants
++++++++

.. autoclass:: ga4gh.protocol.VariantSet
    :members:

.. autoclass:: ga4gh.protocol.CallSet
    :members:

.. autoclass:: ga4gh.protocol.Variant
    :members:

+++++
Reads
+++++

.. autoclass:: ga4gh.protocol.ReadGroupSet
    :members:

.. autoclass:: ga4gh.protocol.ReadGroup
    :members:

.. autoclass:: ga4gh.protocol.ReadAlignment
    :members:

.. autoclass:: ga4gh.protocol.Position
    :members:

----------
Client API
----------

.. todo:: Add overview documentation for the client API.

.. autoclass:: ga4gh.client.HttpClient
    :members: getReferenceSet, getReference,
        getDataset, getVariantSet, getVariant,
        getReadGroupSet, getReadGroup,
        searchDatasets, searchReferenceSets, searchReferences,
        searchVariantSets, searchVariants, searchReadGroupSets,
        searchReads

