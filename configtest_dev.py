"""
Simple shim for running the configuration testing program during development.
"""
import ga4gh.cli

if __name__ == "__main__":
    ga4gh.cli.configtest_main()
