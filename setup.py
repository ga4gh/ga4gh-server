import re
import sys
from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup

# Following the recommendations of PEP 396 we parse the version number 
# out of the module.
def parse_version(module_file):
    """
    Parses the version string from the specified file.
    
    This implementation is ugly, but there doesn't seem to be a good way
    to do this in general at the moment.
    """ 
    f = open(module_file)
    s = f.read()
    f.close()
    match = re.findall("__version__ = '([^']+)'", s)
    return match[0]

ga4gh_version = parse_version("ga4gh/__init__.py") 

requirements = []
v = sys.version_info[:2]
if v < (2, 7) or v == (3, 0) or v == (3, 1):
    requirements.append("argparse")

setup(
    name = "ga4gh",
    version = ga4gh_version,
    packages = ["ga4gh"], 
    scripts = ["scripts/ga4gh_ref.py"],
    install_requires = requirements,
)
