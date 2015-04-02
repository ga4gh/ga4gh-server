.. _configuration:

*************
Configuration
*************

The GA4GH reference server has two basic elements to its configuration:
the `Configuration file`_ and the `Data hierarchy`_.

------------------
Configuration file
------------------

The GA4GH reference server is a `Flask application <http://flask.pocoo.org/>`_
and uses the standard `Flask configuration file mechanisms
<http://flask.pocoo.org/docs/0.10/config/>`_. An example configuration file
might look like::

    DATA_SOURCE = "/path/to/data/root"

--------------
Data hierarchy
--------------

Data is input to the GA4GH server as a directory hierarchy, in which
the structure of data to be served is represented by the filesystem. For now,
we support only one dataset, but this will be generalised to multiple
datasets in later releases. An example data layout might be::

    ga4gh-data/
        /variants/
            variantSet1/
                chr1.vcf.gz
                chr1.vcf.gz.tbi
                chr2.vcf.gz
                chr2.vcf.gz.tbi
                # More VCFs
            variantSet2/
                chr1.bcf
                chr1.bcf.csi
                chr2.bcf
                chr2.bcf.csi
                # More BCFs
        /reads/
            readGroupSet1
                sample1.bam
                sample1.bam.bai
                sample2.bam
                sample2.bam.bai
                # More BAMS

