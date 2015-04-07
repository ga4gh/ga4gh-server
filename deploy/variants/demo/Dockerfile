############################################################
## Dockerfile to build the ga4gh server on mod_wsgi-express
## Downloads and uses the demo dataset
## Based on the prod build
## Results of this build are available from Dockerhub as afirth/ga4gh_server_apache:demo
############################################################
FROM afirth/ga4gh_server_apache:prod

# File Author / Maintainer
MAINTAINER Alastair Firth

# Create the data directory
RUN mkdir -p /data/

# Fetch the example data and extract
# On some terminals the progress bar freezes during build
WORKDIR /data/
RUN curl http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar | tar -x
WORKDIR /

# Set the dataset location explicitly
ENV GA4GH_DATA_SOURCE=/data/ga4gh-example-data
