"""
Sets of configuration values for a GA4GH server.
Any of these values may be overridden by a file on the path
given by the environment variable GA4GH_CONFIGURATION.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


class DefaultConfig(object):
    """
    Simplest default server configuration.
    """
    MAX_CONTENT_LENGTH = 2 * 1024 * 1024  # 2MB
    REQUEST_VALIDATION = False
    RESPONSE_VALIDATION = False
    DATA_SOURCE = "__EMPTY__"

    # Options for the simulated backend.
    SIMULATED_BACKEND_RANDOM_SEED = 0
    SIMULATED_BACKEND_NUM_CALLS = 1
    SIMULATED_BACKEND_VARIANT_DENSITY = 0.5
    SIMULATED_BACKEND_NUM_VARIANT_SETS = 1


class DevelopmentConfig(DefaultConfig):
    """
    Configuration used for development.
    """
    DATA_SOURCE = "ga4gh-example-data"


class TestConfig(DefaultConfig):
    """
    Configuration used in frontend unit tests.
    """
    TESTING = True
    REQUEST_VALIDATION = True
    RESPONSE_VALIDATION = True


class TestWithoutValidationConfig(TestConfig):
    REQUEST_VALIDATION = False
    RESPONSE_VALIDATION = False
