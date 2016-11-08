"""
Reference implementation of the GA4GH APIs.
"""
# Don't include future imports here; we don't want to export them as
# part of the package

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass
