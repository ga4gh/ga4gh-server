"""
Simple shim for running the repo program during development.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import dev_glue  # NOQA
import ga4gh.server.cli.repomanager as cli_repomanager

if __name__ == "__main__":
    cli_repomanager.repo_main()
