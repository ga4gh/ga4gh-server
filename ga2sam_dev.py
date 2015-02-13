"""
Shim for running the ga2sam tool during development
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.cli

if __name__ == "__main__":
    ga4gh.cli.ga2sam_main()
