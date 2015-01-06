from flask import Flask
from flask.ext.api import FlaskAPI
from flask.ext.cors import CORS

app = FlaskAPI(__name__)
CORS(app, allow_headers='Content-Type')

# TODO: Replace hard-coded value with a flag or conifg setting.
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024 # 2MB 

import views
