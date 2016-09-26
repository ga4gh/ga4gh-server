"""
Shim for running the ga2vcf tool during development
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.cli.ga2vcf as cli_ga2vcf

if __name__ == "__main__":
    cli_ga2vcf.ga2vcf_main()
