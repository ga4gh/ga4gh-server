# First, we try to use setuptools. If it's not available locally,
# we fall back on ez_setup.
try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup

with open("README.pypi.rst") as readmeFile:
    long_description = readmeFile.read()

setup(
    name="ga4gh",
    description="A reference implementation of the ga4gh API",
    license='Apache License 2.0',
    long_description=long_description,
    packages=["ga4gh", "ga4gh.datamodel", "ga4gh.templates"],
    include_package_data=True,
    zip_safe=False,
    author="Global Alliance for Genomics and Health",
    author_email="theglobalalliance@genomicsandhealth.org",
    url="https://github.com/ga4gh/server",
    entry_points={
        'console_scripts': [
            'ga4gh_client=ga4gh.cli:client_main',
            'ga4gh_configtest=ga4gh.cli:configtest_main',
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
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='genomics reference',
    # Flask must come after all other requirements that have "flask" as a
    # prefix due to a setuptools bug.
    install_requires=[
        "Flask-Cors==2.0.1",
        "Flask==0.10.1",
        "avro==1.7.7",
        "humanize==0.5.1",
        "pysam==0.8.4",
        "oic==0.7.6",
        "requests==2.7.0",
    ],
    # Use setuptools_scm to set the version number automatically from Git
    setup_requires=['setuptools_scm'],
    use_scm_version={
        "write_to": "ga4gh/_version.py"
    },
)
