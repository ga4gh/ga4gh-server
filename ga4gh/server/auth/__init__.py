import flask
import requests
import functools
import json
import base64

import jwt

from werkzeug.contrib.cache import SimpleCache
cache = SimpleCache()

def authenticate(error):
    """
    Helper function for generating unauthorized responses.
    :param error: Dictionary with details of error
    :return:
    """
    resp = flask.jsonify(error)
    resp.status_code = 401
    return resp

def authorize_email(email='davidcs@ucsc.edu'):
    """
    Adds an email address to the list of authorized users stored in an
    ephemeral cache.
    :param email:
    """
    cache.set(email, {'authorized': True})

def auth_decorator(app=None):
    """
    This decorator wraps a view function so that it is protected when Auth0
    is enabled. This means that any request will be expected to have a signed
    token in the authorization header if the `AUTH0_ENABLED` configuration
    setting is True.

    The authorization header will have the form:

    "authorization: Bearer eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9....."

    If a request is not properly signed, an attempt is made to provide the
    client with useful error messages. This means that if a request is not
    authorized the underlying view function will not be executed.

    When `AUTH0_ENABLED` is false, this decorator will simply execute the
    decorated view without observing the authorization header.
    :param app:
    :return: Flask view decorator
    """
    def requires_auth(f):
        @functools.wraps(f)
        def decorated(*args, **kwargs):
            if app.config.get('AUTH0_ENABLED', False):
                client_id = app.config.get("AUTH0_CLIENT_ID")
                client_secret = app.config.get("AUTH0_CLIENT_SECRET")
                auth = flask.request.headers.get('Authorization', None)
                if not auth:
                  return authenticate(
                      {'code': 'authorization_header_missing',
                       'description': 'Authorization header is expected'})
                parts = auth.split()

                if parts[0].lower() != 'bearer':
                  return authenticate(
                      {'code': 'invalid_header',
                       'description': 'Authorization header must'
                                      ' start with Bearer'})
                elif len(parts) == 1:
                  return authenticate(
                      {'code': 'invalid_header',
                       'description': 'Token not found'})
                elif len(parts) > 2:
                  return authenticate(
                      {'code': 'invalid_header',
                       'description': 'Authorization header must'
                                      ' be Bearer + \s + token'})

                token = parts[1]
                try:
                    b64secret = client_secret.replace("_","/").replace("-","+")
                    payload = jwt.decode(
                        token,
                        base64.b64decode(b64secret),
                        audience=client_id)
                except jwt.ExpiredSignature:
                    return authenticate(
                        {'code': 'token_expired',
                         'description': 'token is expired'})
                except jwt.InvalidAudienceError:
                    message = 'incorrect audience, expected: {}'.format(
                        client_id)
                    return authenticate(
                        {'code': 'invalid_audience',
                         'description': message})
                except jwt.DecodeError:
                    return authenticate(
                        {'code': 'token_invalid_signature',
                         'description': 'token signature is invalid'})
                # We store the token in the session so that later
                # stages can use it to connect identity and
                # authorization.
                flask._request_ctx_stack.top.current_user = user = payload
                flask.session['key'] = token
                # Now we need to make sure that on top of having a good token
                # They are authorized, and if not kick out useful error messages
                if not cache.get(user['email']):
                    message = '{} is not authorized to ' \
                              'access this resource'.format(user['email'])
                    return authenticate(
                        {'code': 'Not authorized','description': message})
                if not cache.get(token):
                    return authenticate(
                        {'code': 'Not authorized',
                         'description': 'The token is good, '
                                        'but you are not logged in.'})
            return f(*args, **kwargs)
        return decorated
    return requires_auth

def logout():
    """
    Logs out the current session by removing it from the cache.
    """
    cache.set(flask.session['key'], None)
    flask.session.clear()
    return True

def callback_maker(domain='david4096.auth0.com',
                   client_id='7PyT5eQRRdaSxgM4hEYz04wzncNiXpqH',
                   client_secret='2OkzpU7LMYMllXf7cd027WFcbBpb4b_'
                                 'p3iNRaUfsYOQgxiYM_puwRxCddCRy_RtV',
                   redirect_uri='http://localhost:5000/callback'):
    """
    This function will generate a view function that can be used to handle
    the return from Auth0. The "callback" is a redirected session from auth0
    that includes the token we can use to authenticate that session.

    If the session is properly authenticated Auth0 will provide a code so our
    application can identify the session. Once this has been done we ask
    for more information about the identified session from Auth0. We then use
    the email of the user logged in to Auth0 to authorize their token to make
    further requests by adding it to the application's cache.

    It sets a value in the cache that sets the current session as logged in. We
    can then refer to this id_token to later authenticate a session.

    :param domain:
    :param client_id:
    :param client_secret:
    :param redirect_uri:
    :return View function:
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
      try:
          token_info = requests.post(
              token_url,
              data=json.dumps(token_payload),
              headers = json_header).json()
          id_token = token_info['id_token']
          access_token = token_info['access_token']
      except Exception as e:
          authenticate({'code': 'bad_callback',
                         'description': 'The callback from Auth0 did not'
                                        'include the expected tokens: \n'
                                        '{}'.format(e.message)})

      try:
          user_url = \
              "https://{domain}/userinfo?access_token={access_token}".format(
                  domain=domain, access_token=access_token)
          user_info = requests.get(user_url).json()
          email = user_info['email']
      except Exception as e:
            authenticate({'code': 'no_user_profile',
                 'description': 'The user profile from Auth0 did not contain'
                                'the expected data: \n'
                                '{}'.format(e.message)})

      user = cache.get(email)
      if user and user['authorized']:
          # TODO set timeout in cache here
          cache.set(id_token, user_info)

      # Redirect to a page that presents the code.
      # TODO make the homepage show the code
      return flask.redirect('/?code={}'.format(id_token))
    return callback_handling

def render_login(app=None, scopes='', redirect_uri='', domain='', client_id=''):
    """
    This function will generate a view function that can be used to handle
    the return from Auth0. The "callback" is a redirected session from auth0
    that includes the token we can use to authenticate that session.

    If the session is properly authenticated Auth0 will provide a code so our
    application can identify the session. Once this has been done we ask
    for more information about the identified session from Auth0. We then use
    the email of the user logged in to Auth0 to authorize their token to make
    further requests by adding it to the application's cache.

    It sets a value in the cache that sets the current session as logged in. We
    can then refer to this id_token to later authenticate a session.

    :param app:
    :param scopes:
    :param redirect_uri:
    :param domain:
    :param client_id:
    :return : Rendered login template
    """
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