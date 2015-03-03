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
    REQUEST_VALIDATION = True
    RESPONSE_VALIDATION = True


class TestConfig(DefaultConfig):
    """
    Configuration used in frontend unit tests.
    """
    TESTING = True


class TestWithoutValidationConfig(TestConfig):
    REQUEST_VALIDATION = False
    RESPONSE_VALIDATION = False
