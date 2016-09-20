"""
Simple shim for running the repo program during development.
"""
import ga4gh.cli.repomanager as cli_repomanager

if __name__ == "__main__":
    cli_repomanager.repo_main()
