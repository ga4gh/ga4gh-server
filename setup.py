import re
from ez_setup import use_setuptools
from setuptools import setup

MIN_SETUPTOOLS_VERSION = "0.7"
use_setuptools(version=MIN_SETUPTOOLS_VERSION)


# Following the recommendations of PEP 396 we parse the version number
# out of the module.
def parseVersion(moduleFile):
    """
    Parses the version string from the specified file.

    This implementation is ugly, but there doesn't seem to be a good way
    to do this in general at the moment.
    """
    f = open(moduleFile)
    s = f.read()
    f.close()
    match = re.findall("__version__ = '([^']+)'", s)
    return match[0]

f = open("README.txt")
ga4ghReadme = f.read()
f.close()
ga4ghVersion = parseVersion("ga4gh/__init__.py")
# Flask must come after all other requirements that have "flask" as a prefix
# due to a setuptools bug.
requirements = ["avro", "Flask-API", "flask-cors", "flask", "humanize",
                "pysam>=0.8.2", "requests"]

setup(
    name="ga4gh",
    version=ga4ghVersion,
    description="A reference implementation of the ga4gh API",
    license='Apache License 2.0',
    long_description=ga4ghReadme,
    packages=["ga4gh", "ga4gh.datamodel", "ga4gh.templates"],
    include_package_data=True,
    zip_safe=False,
    author="Global Alliance for Genomics and Health",
    author_email="theglobalalliance@genomicsandhealth.org",
    url="https://github.com/ga4gh/server",
    entry_points={
        'console_scripts': [
            'ga4gh_client=ga4gh.cli:client_main',
            'ga4gh_server=ga4gh.cli:server_main',
            'ga2vcf=ga4gh.cli:ga2vcf_main',
            'ga2sam=ga4gh.cli:ga2sam_main',
        ]
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',

        # We should add other versions that we can confirm pass the tests
        # (2.6?)
        'Programming Language :: Python :: 2.7',

        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='genomics reference',
    install_requires=requirements,
)
