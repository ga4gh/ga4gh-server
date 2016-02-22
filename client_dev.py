"""
Simple shim for running the client program during development.
"""
import ga4gh.cli_client as cli_client

if __name__ == "__main__":
    cli_client.client_main()
