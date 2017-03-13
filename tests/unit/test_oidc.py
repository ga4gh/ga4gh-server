"""
Unit tests for the oidc code in the frontend.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import socket
import mock
import urlparse
import logging
import mimetypes
import shutil

import oic
import oic.oic.message as message

import ga4gh.server.frontend as frontend
import ga4gh.common.utils as utils


RANDSTR = 'abc'
OICCODE = (
    "5GzAsW4Gfz8Km8PTuayFi6j2tA2VnPrHiWybL%2BkryEpr05D0OUKTGTi49Z54yl%2"
    "FIB5KsUzGwSg%2BFpx%2FpBt0aVc3QsvPlsJu6WGHIchJ55UU%3D")

OICRESPONSE = {
    'scope': 'openid profile',
    'state': RANDSTR,
    'code': '5GzAsW4Gfz8Km8PTuayFi6j2tA2VnPrHiWybL+'
            'kryEpr05D0OUKTGTi49Z54yl/IB5KsUzGwSg+'
            'Fpx/pBt0aVc3QsvPlsJu6WGHIchJ55UU='}

# see https://urllib3.readthedocs.org/en/latest/security.html
# #insecurerequestwarning
# for explanation
logging.captureWarnings(True)


def mockAuthorizationRequest(self, request_args, state):
    class mockResult:
        def __init__(self):
            self.url = None

    result = mockResult()
    urlString = 'http://auth.com/auth?nonce={0}&state={1}&redirect_uri={2}'
    result.url = urlString.format(
        request_args['nonce'], state, request_args['redirect_uri'])
    return result


def mockRndstr(size):
    return RANDSTR


ARESP = {'scope': u'openid profile', 'state': RANDSTR,
         'code': OICRESPONSE['code']}


def mockParseResponse(self, response_class, info, sformat):
    response = message.AuthorizationResponse()
    response.from_dict(ARESP)
    return response


ATR = {
    'access_token': '5GzAsW4Gfz8Km8PTuayFi6j2tA2VnPrHiWybL+kryEpr05D0OUKT'
                    'GTi49Z54yl/ImAwqKJ28BWaX6Bvzwj+7uD4udUz'
                    'gFh7bRlFXaQsjz+M=',
    'id_token': {'nonce': RANDSTR,
                 'sub': '-2c68a9e99c6eb6da',
                 'iss': 'https://localhost:8443/',
                 'acr': 'password',
                 'exp': 1435138548,
                 'auth_time': 1435055251,
                 'iat': 1435055748,
                 'aud': ['aKRswYgOAiye']},
    'expires_in': 3600,
    'token_type': 'Bearer',
    'state': RANDSTR,
    'scope': 'openid profile',
    'refresh_token': '5GzAsW4Gfz8Km8PTuayFi6j2tA2VnPrHiWybL+kryEpr05D0OUK'
                     'TGTi49Z54yl/IRkog8/WippmBGCn5UR1'
                     '/W3fxa3wuWQxcjL+4B9BoqHY='}


def mockAccessTokenRequest(self, scope, state, request_args):
    response = message.AccessTokenResponse()
    response.from_dict(ATR)
    return response


def mockAccessTokenRequestError(self, scope, state, request_args):
    return 'random string that is not a message.AccessTokenResponse'


class TestFrontendOidc(unittest.TestCase):
    """
    Tests the OIDC support for the Flask app.
    """

    @classmethod
    def setUpClass(cls):
        config = {
            "SIMULATED_BACKEND_RANDOM_SEED": 1111,
            "SIMULATED_BACKEND_NUM_CALLS": 1,
            "SIMULATED_BACKEND_VARIANT_DENSITY": 1.0,
            "SIMULATED_BACKEND_NUM_VARIANT_SETS": 1,
            "DATA_SOURCE": "simulated://",
            "OIDC_CLIENT_ID": "123",
            "OIDC_CLIENT_SECRET": RANDSTR,
            "OIDC_PROVIDER": "https://accounts.example.com",
            "SECRET_KEY": "secret"
            # "DEBUG" : True
        }
        shutil.rmtree('/tmp/ga4gh', True)
        frontend.reset()
        frontend.configure(
            baseConfig="TestOidcConfig", extraConfig=config, port=8001)
        cls.app = frontend.app.test_client()

    @classmethod
    def tearDownClass(cls):
        cls.app = None
        shutil.rmtree('/tmp/ga4gh', True)

    def sendPostRequest(self, path, request):
        """
        Sends the specified GA request object and returns the response.
        """
        versionedPath = utils.applyVersion(path)
        headers = {
            'Content-type': 'application/json',
            'Origin': self.exampleUrl,
        }
        return self.app.post(
            versionedPath, headers=headers,
            data=request.toJsonString())

    @mock.patch('oic.oauth2.rndstr', mockRndstr)
    @mock.patch('oic.oic.Client.do_authorization_request',
                mockAuthorizationRequest)
    def testAuthorizationRequest(self):
        """
        Test that when the index is referenced, we get a redirect to the
        authorization provider, of the correct form.
        """
        response = self.app.get('/')
        scheme, netloc, path, params, query, fragment = urlparse.urlparse(
            response.location)
        self.assertEqual(netloc, 'auth.com')
        args = urlparse.parse_qs(query, strict_parsing=True)
        self.assertEqual(path, '/auth')
        self.assertEqual(args['redirect_uri'][0],
                         'https://{}:8001/oauth2callback'.format(
                             socket.gethostname()))

    @mock.patch('oic.oauth2.rndstr', lambda x: RANDSTR)
    @mock.patch('oic.oic.Client.parse_response',
                mockParseResponse)
    @mock.patch('oic.oic.Client.do_access_token_request',
                mockAccessTokenRequest)
    def testOidcCallbackSucceeds(self):
        """
        Test that when the authorization provider calls back to us, we can
        handle the reply correctly. This test is the 'succeeds' case.
        """
        url = '/oauth2callback?scope=openid+profile&state={0}&code={1}'.format(
            oic.oauth2.rndstr(0), OICCODE
        )
        with self.app.session_transaction() as sess:
            sess['state'] = oic.oauth2.rndstr(0)
            sess['nonce'] = oic.oauth2.rndstr(0)
        result = self.app.get(url)
        self.assertEqual(result.status_code, 302)
        self.assertEqual(result.location, 'http://{0}:8001/'.format(
            socket.gethostname()))

    def testSessionKeyAllowsIndex(self):
        """
        Test that if we have a valid session in play, we retrieve the index
        page
        """
        with self.app as app:
            with app.session_transaction() as sess:
                sess['key'] = 'xxx'
            app.application.cache.set('xxx', RANDSTR)
            result = app.get('/')
            self.assertEqual(result.status_code, 200)
            self.assertEqual("text/html", result.mimetype)
            self.assertGreater(len(result.data), 0)

    def testKeyParamAllowsIndex(self):
        """
        Test that if we have a valid key in play, we retrieve the index
        page
        """
        with self.app as app:
            app.application.cache.set('xxx', RANDSTR)
            result = app.get('/?key=xxx')
            self.assertEqual(result.status_code, 200)
            self.assertEqual("text/html", result.mimetype)
            self.assertGreater(len(result.data), 0)

    @mock.patch('oic.oauth2.rndstr', mockRndstr)
    @mock.patch('oic.oic.Client.do_authorization_request',
                mockAuthorizationRequest)
    def testInvalidSessionKeyRedirects(self):
        """
        Test that if we have an invalid session in play, we redirected to OP
        """
        with self.app as app:
            with app.session_transaction() as sess:
                sess['key'] = 'xxx'
            result = app.get('/')
            scheme, netloc, path, params, query, fragment = urlparse.urlparse(
                result.location)
            self.assertEqual(netloc, 'auth.com')
            args = urlparse.parse_qs(query, strict_parsing=True)
            self.assertEqual(path, '/auth')
            self.assertEqual(args['redirect_uri'][0],
                             'https://{}:8001/oauth2callback'.format(
                                 socket.gethostname()))

    @mock.patch('oic.oauth2.rndstr', mockRndstr)
    @mock.patch('oic.oic.Client.do_authorization_request',
                mockAuthorizationRequest)
    def testInvalidKeyParamIsUnauthorized(self):
        """
        Test that if we have an invalid key in play, we get an exception
        """
        result = self.app.get('/?key=xxx')
        self.assertEqual(result.status_code, 403)

    @mock.patch('oic.oauth2.rndstr', lambda x: RANDSTR)
    @mock.patch('oic.oic.Client.parse_response',
                mockParseResponse)
    @mock.patch('oic.oic.Client.do_access_token_request',
                mockAccessTokenRequest)
    def testOidcCallbackBadNonce(self):
        """
        Test that when the authorization provider calls back to us, we can
        handle the reply correctly.
        In this case, the nonce returned does not match that in the session.
        """
        url = '/oauth2callback?scope=openid+profile&state={0}&code={1}'.format(
            oic.oauth2.rndstr(0), OICCODE
        )
        with self.app as app:
            with app.session_transaction() as sess:
                sess['state'] = oic.oauth2.rndstr(0)
                sess['nonce'] = 'other'
            result = app.get(url)
            self.assertEqual(result.status_code, 403)

    @mock.patch('oic.oauth2.rndstr', lambda x: RANDSTR)
    @mock.patch('oic.oic.Client.parse_response',
                mockParseResponse)
    @mock.patch('oic.oic.Client.do_access_token_request',
                mockAccessTokenRequest)
    def testOidcCallbackBadState(self):
        """
        Test that when the authorization provider calls back to us, we can
        handle the reply correctly.
        In this case, the state returned does not match that in the session.
        """
        url = '/oauth2callback?scope=openid+profile&state={0}&code={1}'.format(
            oic.oauth2.rndstr(0), OICCODE
        )
        with self.app as app:
            with app.session_transaction() as sess:
                sess['state'] = 'other'
                sess['nonce'] = oic.oauth2.rndstr(0)
            result = app.get(url)
            self.assertEqual(result.status_code, 403)

    @mock.patch('oic.oauth2.rndstr', lambda x: RANDSTR)
    @mock.patch('oic.oic.Client.parse_response',
                mockParseResponse)
    @mock.patch('oic.oic.Client.do_access_token_request',
                mockAccessTokenRequestError)
    def testOidcCallbackBadTokenResponse(self):
        """
        Test that when the authorization provider calls back to us, we can
        handle the reply correctly.
        In this case, the access token request returns an invalid response.
        """
        url = '/oauth2callback?scope=openid+profile&state={0}&code={1}'.format(
            oic.oauth2.rndstr(0), OICCODE
        )
        with self.app as app:
            with app.session_transaction() as sess:
                sess['state'] = oic.oauth2.rndstr(0)
                sess['nonce'] = oic.oauth2.rndstr(0)
            result = app.get(url)
            self.assertEqual(result.status_code, 403)

    def testMimetype(self):
        """
        The  `oidc-provider` relies on mimetypes.guess_type(full_path)
        `mimetypes` in turn relies on the host OS configuration,
        typically ["/etc/mime.types","/etc/httpd/conf/mime.types",...]
        If this test fails, see the results of `mimetypes.knownfiles`
        on your host and ensure configuration set up
        see https://github.com/ga4gh/ga4gh-server/issues/501
        """
        content_type, encoding = mimetypes.guess_type("foo.json")
        self.assertEqual(content_type, "application/json")
