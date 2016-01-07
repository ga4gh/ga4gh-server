"""
Reference implementation of the GA4GH APIs.
"""

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass
