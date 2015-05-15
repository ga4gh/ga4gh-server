#!/bin/bash

virtualenv .
. bin/activate
pip install -r simple_op/requirements.txt
cd simple_op && python src/run.py -p 8443 settings.yaml

# This is how you would register a client
C=$(http --verify=no POST https://localhost:8443/registration redirect_uris:='["https://localhost:8444/cb"]')
