import flask
import requests
import functools
import json
import base64

import jwt

from werkzeug.contrib.cache import SimpleCache
cache = SimpleCache()

# In order to make it so an endpoint is protected we need to do a couple of things.
# The first is a little helper function that will kick back authorization failures
# with hopefully helpful error codes.
def authenticate(error):
  resp = flask.jsonify(error)
  resp.status_code = 401
  return resp

def authorize_email(email='davidcs@ucsc.edu'):
    cache.set(email, {'authorized': True})

# This decorator wraps any function that needs to have an authorized header to be
# served from. It just checks to see if the key that is being used is still good.
def auth_decorator(app=None):
    def requires_auth(f):
        @functools.wraps(f)
        def decorated(*args, **kwargs):
            if app.config.get('AUTH0_ENABLED', False):
                client_id = app.config.get("AUTH0_CLIENT_ID")
                client_secret = app.config.get("AUTH0_CLIENT_SECRET")
                auth = flask.request.headers.get('Authorization', None)
                if not auth:
                  return authenticate({'code': 'authorization_header_missing', 'description': 'Authorization header is expected'})
                parts = auth.split()

                if parts[0].lower() != 'bearer':
                  return authenticate({'code': 'invalid_header', 'description': 'Authorization header must start with Bearer'})
                elif len(parts) == 1:
                  return authenticate({'code': 'invalid_header', 'description': 'Token not found'})
                elif len(parts) > 2:
                  return authenticate({'code': 'invalid_header', 'description': 'Authorization header must be Bearer + \s + token'})

                token = parts[1]
                try:
                    payload = jwt.decode(
                        token,
                        base64.b64decode(client_secret.replace("_","/").replace("-","+")),
                        audience=client_id)
                except jwt.ExpiredSignature:
                    return authenticate({'code': 'token_expired', 'description': 'token is expired'})
                except jwt.InvalidAudienceError:
                    return authenticate({'code': 'invalid_audience', 'description': 'incorrect audience, expected: {}'.format(
                        client_id)})
                except jwt.DecodeError:
                    return authenticate({'code': 'token_invalid_signature', 'description': 'token signature is invalid'})
                # We store the token in the session so that later stages can use it to connect identity and
                # authorization.
                flask.session['key'] = token
                # Now we need to make sure that on top of having a good token
                # They are authorized, and if not kick out useful error messages
                flask._request_ctx_stack.top.current_user = user = payload
                if not cache.get(user['email']):
                    return authenticate({'code': 'Not authorized', 'description': '{} is not authorized to access this resource'.format(user['email'])})
                if not cache.get(token):
                    return authenticate({'code': 'Not authorized', 'description': 'The token is good, but you are not logged in.'})
            return f(*args, **kwargs)
    return requires_auth

def logout():
    cache.set(flask.session['key'], None)
    flask.session.clear()
    # TODO guard against key replay attacks

def callback_maker(domain='david4096.auth0.com',
                   client_id='7PyT5eQRRdaSxgM4hEYz04wzncNiXpqH',
                   client_secret='2OkzpU7LMYMllXf7cd027WFcbBpb4b_p3iNRaUfsYOQgxiYM_puwRxCddCRy_RtV',
                   redirect_uri='http://localhost:5000/callback'):
    """
    This function will generate the callback handler that can be used to handle the return from Auth0.
    It sets a value in the cache that sets the current user to being logged in.

    :param domain:
    :param client_id:
    :param client_secret:
    :param redirect_uri:
    :return:
    """
    def callback_handling():
      code = flask.request.args.get('code')
      json_header = {'content-type': 'application/json'}
      token_url = "https://{domain}/oauth/token".format(domain=domain)
      token_payload = {
        'client_id':     client_id,
        'client_secret': client_secret,
        'redirect_uri':  redirect_uri,
        'code':          code,
        'grant_type':    'authorization_code'
      }
      token_info = requests.post(token_url, data=json.dumps(token_payload), headers = json_header).json()
      # TODO failures against auth0 should be properly handled here
      print(token_info)
      user_url = "https://{domain}/userinfo?access_token={access_token}" \
          .format(domain=domain, access_token=token_info['access_token'])
      user_info = requests.get(user_url).json()

      # there are failure modes for the callback in which case we won't provide a key
      # or redirect
      user = cache.get(user_info['email'])
      print("attempting to login as {} {}".format(user_info['email'], user))
      print(user_info)
      if user and user['authorized']:
          cache.set(token_info['id_token'], user_info)

      # Redirect to a page that presents the code.
      return flask.redirect('/?code={}'.format(token_info['id_token']))

    return callback_handling

def render_login(app=None, scopes='', redirect_uri='', domain='', client_id=''):
    return app.jinja_env.from_string(LOGIN_HTML).render(
        scopes=scopes,
        redirect_uri=redirect_uri,
        domain=domain,
        client_id=client_id)

LOGIN_HTML = """<html>
<head>

  <title>Log in</title></head><body><div>
    <script src="https://cdn.auth0.com/js/lock/10.0/lock.min.js"></script>
    <script type="text/javascript">
      var lock = new Auth0Lock('{{ client_id }}', '{{ domain }}', {
        auth: {
          redirectUrl: '{{ redirect_uri }}',
          responseType: 'code',
          params: {
            scope: '{{ scopes }}' // Learn about scopes: https://auth0.com/docs/scopes
          }
        }
      });
    lock.show();
    </script>
    </div>"""