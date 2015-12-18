#!/bin/bash

HOSTNAME=`python -c 'import socket; print socket.gethostname()'`
cd simple_op && python src/run.py --base https://${HOSTNAME}:8443 -p 8443 -d settings.yaml

# This is how you would register a client. Uses the httpie package.
#C=$(http --verify=no POST https://localhost:8443/registration redirect_uris:='["https://localhost:8080/oauth2callback"]')

