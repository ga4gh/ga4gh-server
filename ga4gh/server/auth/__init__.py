"""
Helps to implement authentication and authorization using Auth0.

Offers functions for generating the view functions needed to implement Auth0,
a login screen, callback maker, and a function decorator for protecting
endpoints.
"""
import flask
import requests
import functools
import json

import jwt

import ga4gh.server.exceptions as exceptions


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
            # This decorator will only apply with AUTH0_ENABLED set to True.
            if app.config.get('AUTH0_ENABLED', False):
                client_id = app.config.get("AUTH0_CLIENT_ID")
                client_secret = app.config.get("AUTH0_CLIENT_SECRET")
                auth_header = flask.request.headers.get('Authorization', None)
                # Each of these functions will throw a 401 is there is a
                # problem decoding the token with some helpful error message.
                if auth_header:
                    token, profile = decode_header(
                        auth_header, client_id, client_secret)
                else:
                    raise exceptions.NotAuthorizedException()
                # We store the token in the session so that later
                # stages can use it to connect identity and authorization.
                flask.session['auth0_key'] = token
                # Now we need to make sure that on top of having a good token
                # They are authorized, and if not provide an error message
                is_authorized(app.cache, profile['email'])
                is_active(app.cache, token)
            return f(*args, **kwargs)
        return decorated
    return requires_auth


def decode_header(auth_header, client_id, client_secret):
    """
    A function that threads the header through decoding and returns a tuple
    of the token and payload if successful. This does not fully authenticate
    a request.
    :param auth_header:
    :param client_id:
    :param client_secret:
    :return: (token, profile)
    """
    return _decode_header(
        _well_formed(
            _has_token(_has_bearer(_has_header(auth_header)))),
        client_id, client_secret)


def logout(cache):
    """
    Logs out the current session by removing it from the cache. This is
    expected to only occur when a session has
    """
    cache.set(flask.session['auth0_key'], None)
    flask.session.clear()
    return True


def callback_maker(
        cache=None, domain='', client_id='',
        client_secret='', redirect_uri=''):
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
    :return : View function
    """
    def callback_handling():
        code = flask.request.args.get('code')
        if code is None:
            raise exceptions.NotAuthorizedException(
                'The callback expects a well '
                'formatted code, {} was provided'.format(code))
        json_header = {'content-type': 'application/json'}
        # Get auth token
        token_url = "https://{domain}/oauth/token".format(domain=domain)
        token_payload = {
            'client_id':     client_id,
            'client_secret': client_secret,
            'redirect_uri':  redirect_uri,
            'code':          code,
            'grant_type':    'authorization_code'}
        try:
            token_info = requests.post(
                token_url,
                data=json.dumps(token_payload),
                headers=json_header).json()
            id_token = token_info['id_token']
            access_token = token_info['access_token']
        except Exception as e:
            raise exceptions.NotAuthorizedException(
                'The callback from Auth0 did not'
                'include the expected tokens: \n'
                '{}'.format(e.message))
        # Get profile information
        try:
            user_url = \
              "https://{domain}/userinfo?access_token={access_token}".format(
                  domain=domain, access_token=access_token)
            user_info = requests.get(user_url).json()
            email = user_info['email']
        except Exception as e:
            raise exceptions.NotAuthorizedException(
                'The user profile from Auth0 did '
                'not contain the expected data: \n {}'.format(e.message))
        # Log token in
        user = cache.get(email)
        if user and user['authorized']:
            cache.set(id_token, user_info)
            return flask.redirect('/login?code={}'.format(id_token))
        else:
            return flask.redirect('/login')
    return callback_handling


def render_login(
        app=None, scopes='', redirect_uri='', domain='', client_id=''):
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


def render_key(app, key=""):
    """
    Renders a view from the app and a key that lets the current session grab
    its token.
    :param app:
    :param key:
    :return: Rendered view
    """
    return app.jinja_env.from_string(KEY_HTML).render(
        key=key)


def authorize_email(email='davidcs@ucsc.edu', cache=None):
    """
    Adds an email address to the list of authorized emails stored in an
    ephemeral cache.
    :param email:
    """
    # TODO safely access cache
    cache.set(email, {'authorized': True})


def _has_header(auth_header):
    if not auth_header:
        raise exceptions.NotAuthorizedException(
            'Authorization header is expected.')
    return auth_header


def _has_bearer(auth_header):
    parts = auth_header.split()
    if parts[0].lower() != 'bearer':
        raise exceptions.NotAuthorizedException(
            'Authorization header must start with "Bearer".')
    return auth_header


def _has_token(auth_header):
    parts = auth_header.split()
    if len(parts) == 1:
        raise exceptions.NotAuthorizedException(
           'Token not found in header.')
    return auth_header


def _well_formed(auth_header):
    parts = auth_header.split()
    if len(parts) > 2:
        raise exceptions.NotAuthorizedException(
            'Authorization header must be Bearer + \s + token.')
    return auth_header


def _decode_header(auth_header, client_id, client_secret):
    """
    Takes the header and tries to return an active token and decoded
    payload.
    :param auth_header:
    :param client_id:
    :param client_secret:
    :return: (token, profile)
    """
    try:
        token = auth_header.split()[1]
        payload = jwt.decode(
            token,
            client_secret,
            audience=client_id)
    except jwt.ExpiredSignature:
        raise exceptions.NotAuthorizedException(
            'Token has expired, please log in again.')
    # is valid client
    except jwt.InvalidAudienceError:
        message = 'Incorrect audience, expected: {}'.format(
            client_id)
        raise exceptions.NotAuthorizedException(message)
    # is valid token
    except jwt.DecodeError:
        raise exceptions.NotAuthorizedException(
            'Token signature could not be validated.')
    except Exception as e:
        raise exceptions.NotAuthorizedException(
            'Token signature was malformed. {}'.format(e.message))
    return token, payload


def is_authorized(cache, email):
    if not cache.get(email):
        message = '{} is not authorized to ' \
                  'access this resource'.format(email)
        raise exceptions.NotAuthenticatedException(message)
    return email


def is_active(cache, token):
    """
    Accepts the cache and ID token and checks to see if the profile is
    currently logged in. If so, return the token, otherwise throw a
    NotAuthenticatedException.
    :param cache:
    :param token:
    :return:
    """
    profile = cache.get(token)
    if not profile:
        raise exceptions.NotAuthenticatedException(
            'The token is good, but you are not logged in. Please '
            'try logging in again.')
    return profile


# This HTML string is used to render the login page. It is a jinja template.
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
            scope: '{{ scopes }}' // https://auth0.com/docs/scopes
          }
        }
      });
    lock.show();
    </script>
    </div>"""

KEY_HTML = """<html>
<head>
<title>GA4GH Server API Token</title></head><body><div>
<h1>Your API Token</h1>
<p>Your token is now active, add it as your "Authorization: bearer $TOKEN" header
when making requests to protected endpoints</p>
<textarea cols=120 rows=5 onClick='this.select()' readonly>{{ key }}</textarea>
<h3><a href="/?key={{ key }}">Visit landing page</a></h3>
</div>
""" # noqa
