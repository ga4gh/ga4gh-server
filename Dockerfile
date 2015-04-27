############################################################
## Dockerfile to build the ga4gh server on mod_wsgi-express
## Configurable to use a local dataset
## Based on mod_wsgi-docker
## Results of this build are available from Dockerhub as afirth/ga4gh_server_apache:prod
############################################################
FROM grahamdumpleton/mod-wsgi-docker:python-2.7-onbuild

# File Author / Maintainer
MAINTAINER Alastair Firth

# Place the config, if an existing config is not in place
ADD deploy/config.py /app/ga4gh/docker_config.py
RUN cp --no-clobber /app/ga4gh/docker_config.py /app/ga4gh/config.py

# Pass '-e GA4GH_DATA_SOURCE=/container/data/path' and '-v /host/data/path:/container/data/path:ro' to docker run to mount local data
# See docs for more info

CMD [ "--working-directory", "ga4gh", \
      "--log-to-terminal", \
      "ga4gh/application.wsgi" ]
