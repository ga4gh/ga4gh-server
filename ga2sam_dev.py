"""
Shim for running the ga2sam tool during development
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.cli.ga2sam as ga2sam

if __name__ == "__main__":
    ga2sam.ga2sam_main()
