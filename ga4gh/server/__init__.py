from flask import Flask
from flask.ext.cors import CORS

app = Flask(__name__)
CORS(app, allow_headers='Content-Type')

app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024 # 2MB 

import views
