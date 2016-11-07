"""
Simple shim for running the configuration testing program during development.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import dev_glue  # NOQA
import ga4gh.server.cli.configtest as cli_configtest

if __name__ == "__main__":
    cli_configtest.configtest_main()
