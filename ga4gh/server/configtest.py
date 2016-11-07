"""
Tests for the server configuration. Run via configtest_dev.py

Tests are standard python unittest tests.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import flask
import os


class TestConfig(unittest.TestCase):
    """
    This class implements configuration checks for the GA4GH server. It is
    structured as a set of unit tests, and relies on the unittest module.
    However, as some of the tests require state injection (through command
    line arguments) we need to break the usual unittest pattern and expose
    some configuration as class attributes that can be set before tests are
    run.
    """

    configStr = None
    configEnv = None
    configFile = None

    def setUp(self):
        """
        Create a flask app to be used in the tests. Use the supplied class
        attributes for configuration
        """
        self.app = flask.Flask(__name__)
        try:
            self.app.config.from_object(self.configStr)
        except Exception:
            self.fail("Cannot create app from {0}".format(self.configStr))
        if os.environ.get(self.configEnv) is not None:
            self.app.config.from_envvar(self.configEnv)
        if self.configFile is not None:
            self.assertTrue(os.path.exists(self.configFile),
                            "config file does not exist")
            self.app.config.from_pyfile(self.configFile)

    def test_config_parameters(self):
        """
        A simple test that configuration parameters are of the correct types
        and ranges
        """
        self.assertIsInstance(self.app.config["REQUEST_VALIDATION"],
                              bool,
                              'REQUEST_VALIDATION is not a boolean')
        self.assertIsNotNone(self.app.config["DATA_SOURCE"],
                             'DATA_SOURCE must be set')
