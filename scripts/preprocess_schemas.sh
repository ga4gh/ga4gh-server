#!/bin/bash

# remove old schemas
rm -r ../schemas

# download the avro-tools jar
wget http://www.carfab.com/apachesoftware/avro/stable/java/avro-tools-1.7.7.jar

# download schema files
wget https://raw.githubusercontent.com/ga4gh/schemas/master/src/main/resources/avro/common.avdl
wget https://raw.githubusercontent.com/ga4gh/schemas/master/src/main/resources/avro/variants.avdl

# extract schemata
java -jar avro-tools-1.7.7.jar idl2schemata variants.avdl

# cd down, and make a schemas directory
cd ..
mkdir -p schemas

# move schemas to schema directory
mv scripts/*.avdl schemas
mv scripts/*.avsc schemas

# remove avro-tools jarfile
rm scripts/avro-tools-1.7.7.jar