############################################################
# Dockerfile to build the ga4gh demo 
# Based on CentOS7
############################################################
FROM centos:7

# File Author / Maintainer
MAINTAINER Alastair Firth


################## BEGIN INSTALLATION ######################
# Install ga4gh demo (non-apache)
# NOT for production use
# Ref: http://ga4gh-reference-implementation.readthedocs.org/en/stable/demo.html#demo

#TODO add yum key

# Install prereqs
RUN yum -y update && yum -y install \
  curl \
  git \
  python-backports \
  python-backports-ssl_match_hostname \
  python-devel \
  python-setuptools \
  tar \
  wget \
  zlib-devel
RUN yum -y groupinstall "Development Tools"
RUN yum clean all

# Create the default data directory
RUN mkdir -p /data/
WORKDIR /data/

#wget the example data
RUN wget http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar -nv -O - | tar -x

# Install setuptools and pip
RUN easy_install setuptools
RUN easy_install pip

# Install ga4gh with pip
RUN pip install ga4gh --pre

##################### INSTALLATION END #####################

# Expose the default port
EXPOSE 8000

ENTRYPOINT ["ga4gh_server"]

