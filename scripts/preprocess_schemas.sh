#!/bin/bash -x

# Create tmp dir and cd into it
PRJROOT=$PWD
TMPDIR=$( mktemp -d )
cd $TMPDIR

echo $TMPDIR

# download the avro-tools jar
wget http://www.carfab.com/apachesoftware/avro/stable/java/avro-tools-1.7.7.jar

# download schema files
wget https://raw.githubusercontent.com/ga4gh/schemas/master/src/main/resources/avro/common.avdl
wget https://raw.githubusercontent.com/ga4gh/schemas/master/src/main/resources/avro/variants.avdl

# extract schemata
java -jar avro-tools-1.7.7.jar idl2schemata variants.avdl

SCHEMA_DIR=schemas
# make a schemas directory if it does not exist
cd $PRJROOT
mkdir -p $SCHEMA_DIR

# remove old schemas
rm -fr $SCHEMA_DIR/*

# move schemas to schema directory
mv $TMPDIR/*.avdl $SCHEMA_DIR 
mv $TMPDIR/*.avsc $SCHEMA_DIR

rm -fR $TMPDIR
