"""
Simple shim for running the server program during development.
"""
import ga4gh.cli.server as cli_server

if __name__ == "__main__":
    cli_server.server_main()
