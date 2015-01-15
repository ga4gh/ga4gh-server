from flask.ext.api import FlaskAPI
from flask.ext.cors import CORS
import os

app = FlaskAPI(__name__)

app.config.from_object('ga4gh.server.config:DefaultConfig')
if os.environ.get('GA4GH_CONFIGURATION') is not None:
    app.config.from_envvar('GA4GH_CONFIGURATION')

CORS(app, allow_headers='Content-Type')

import views

__all__ = [views]
