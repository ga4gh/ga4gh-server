"""
Simple shim for running the configuration testing program during development.
"""
import ga4gh.cli.configtest as cli_configtest

if __name__ == "__main__":
    cli_configtest.configtest_main()
