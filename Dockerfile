############################################################
## Dockerfile to build the ga4gh server on mod_wsgi-express
## Configurable to use a local dataset
############################################################
FROM ubuntu

# Originally created by Steve Hershman GitHub @hershman
# previously maintained by Alastair Firth, and Maciek Smuga-Otto of the
# UCSC Genomics Institute
MAINTAINER David Steinberg <david@resium.com>

# Update the sources list
RUN apt-get update  --fix-missing
RUN apt-get upgrade --yes

# Install packages
RUN apt-get install -y tar git curl libcurl4-openssl-dev wget dialog \
    net-tools build-essential python python-dev python-distribute \
    python-pip zlib1g-dev apache2 libapache2-mod-wsgi libxslt1-dev \
    libffi-dev libssl-dev

# Enable wsgi module
RUN a2enmod wsgi

# Create cache directories
RUN mkdir /var/cache/apache2/python-egg-cache && \
    chown www-data:www-data /var/cache/apache2/python-egg-cache/

# build the GA4GH server
RUN mkdir -p /srv/ga4gh/server
WORKDIR /srv/ga4gh/server

# Configure the python requirements
# Do this as a separate step prior to the build so that changes
# to the GA4GH Server codebase do not trigger a full rebuild of the
# pip requirements.
COPY requirements.txt /srv/ga4gh/server/
COPY constraints.txt /srv/ga4gh/server/
RUN pip install -r requirements.txt -c constraints.txt

# Install the code
COPY . /srv/ga4gh/server/
RUN pip install . -c constraints.txt

# Write new apache config
COPY deploy/001-ga4gh.conf /etc/apache2/sites-available/001-ga4gh.conf

# Write application.wsgi
COPY deploy/application.wsgi /srv/ga4gh/application.wsgi
COPY deploy/config.py /srv/ga4gh/config.py

# Configure apache to serve GA4GH site
WORKDIR /etc/apache2/sites-enabled
RUN a2dissite 000-default
RUN a2ensite 001-ga4gh

# Open port 80 for HTTP
EXPOSE 80

# Prepare container for deployment
# The directory that the user will land in when executing an interactive shell
WORKDIR /srv/ga4gh/server
RUN python scripts/prepare_compliance_data.py -o ../ga4gh-compliance-data

# Default action: Bring up a webserver instance to run as a daemon
CMD ["/usr/sbin/apache2ctl", "-D", "FOREGROUND"]
