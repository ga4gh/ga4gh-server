"""
Shim for running the ga2sam tool during development
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.cli_ga2sam as cli_ga2sam

if __name__ == "__main__":
    cli_ga2sam.ga2sam_main()
