# Base image
FROM ubuntu

# Originally created by Steve Hershman GitHub @hershman
MAINTAINER Steve Hershman

# Update the sources list
RUN apt-get update

# Install packages
RUN apt-get install -y tar git curl nano wget dialog net-tools build-essential python python-dev python-distribute python-pip zlib1g-dev python-virtualenv apache2 libapache2-mod-wsgi

# Enable wsgi module
RUN a2enmod wsgi

# Create cache directories
RUN mkdir /var/cache/apache2/python-egg-cache && \
    chown www-data:www-data /var/cache/apache2/python-egg-cache/

# Set up GA4GH server
RUN mkdir /srv/ga4gh
WORKDIR /srv/ga4gh
RUN virtualenv ga4gh-server-env
RUN /bin/bash -c "source ga4gh-server-env/bin/activate"
RUN pip install --pre ga4gh

# Install relevant sample data
RUN wget http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar && \
    tar -xvf ga4gh-example-data.tar

# Write application.wsgi
RUN echo "from ga4gh.frontend import app as application\nimport ga4gh.frontend as frontend\nfrontend.configure(\"/srv/ga4gh/config.py\")" > application.wsgi

# Write config.py
RUN echo "DATA_SOURCE = \"/srv/ga4gh/ga4gh-example-data\"" > config.py

# Write new apache config
WORKDIR /etc/apache2/sites-available
RUN echo "<VirtualHost *:80>\n  ServerAdmin webmaster@localhost\n  DocumentRoot /var/www/html\n  ErrorLog ${APACHE_LOG_DIR}/error.log\n  CustomLog ${APACHE_LOG_DIR}/access.log combined\n  WSGIDaemonProcess ga4gh python-path=/srv/ga4gh/ga4gh-server-env/lib/python2.7/site-packages python-eggs=/var/cache/apache2/python-egg-cache\n  WSGIScriptAlias /ga4gh /srv/ga4gh/application.wsgi\n  <Directory /srv/ga4gh>\n    WSGIProcessGroup ga4gh\n    WSGIApplicationGroup %{GLOBAL}\n    Require all granted\n  </Directory>\n</VirtualHost>" > 001-ga4gh.conf

# Configure apache to serve GA4GH site
WORKDIR /etc/apache2/sites-enabled
RUN rm -f 000-default.conf && ln -s /etc/apache2/sites-available/001-ga4gh.conf 001-ga4gh.conf

# Open port 80
EXPOSE 80

# Start server when container starts
CMD ["/usr/sbin/apache2ctl", "-D", "FOREGROUND"]
