# This file needs to be copied to /app/ga4gh/config.py by default
import os

#TODO this logic could move to frontend.configure() or BaseConfig
# For docker, if/when a default is set in serverconfig.py, run -v /localdata:/default-path
# If the env variable GA4GH_DATA_SOURCE is set, use that path. Otherwise, use the default path
DATA_SOURCE = os.getenv('GA4GH_DATA_SOURCE', "/ga4gh-example-data")

# If the env variable GA4GH_DEBUG is set, use that. Otherwise, use the empty string (False)
# Enable with True
DEBUG = os.getenv('GA4GH_DEBUG', "")
