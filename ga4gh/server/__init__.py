"""
The GA4GH HTTP frontend.

TODO: document this properly.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from flask.ext.api import FlaskAPI
from flask.ext.cors import CORS

app = FlaskAPI(__name__)
CORS(app, allow_headers='Content-Type')

# TODO: Replace hard-coded value with a flag or conifg setting.
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024  # 2MB

import views

__all__ = [views]
