"""
Simple shim for running the server program during development.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import dev_glue  # NOQA
import ga4gh.server.cli.server as cli_server

if __name__ == "__main__":
    cli_server.server_main()
