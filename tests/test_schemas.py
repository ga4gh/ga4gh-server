"""
Test loading the schemas.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import avro.schema
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter
import os
import subprocess
import unittest


class TestLoadSchemas(unittest.TestCase):

    def setUp(self):

        # have we already generated the schemas? if not, generate them
        if not os.path.isdir("schemas"):
            # generate schemas
            subprocess.call("./scripts/preprocess_schemas.sh")

    def testLoadSchema(self):

        # load schema
        schema = avro.schema.parse(open("schemas/GAVariant.avsc").read())

        # write schemas
        writer = DataFileWriter(
            open("variants.avro", "w"), DatumWriter(), schema)
        writer.append({"id": "myID",
                       "variantSetId": "mySetID",
                       "names": ["myName"],
                       "created": 100L,
                       "updated": None,
                       "referenceName": "chr1",
                       "start": 1L,
                       "end": 2L,
                       "referenceBases": "A",
                       "alternateBases": [],
                       "info": {"infoKey": ["myInfo"]},
                       "calls": []})
        writer.close()

        # read schemas back
        reader = DataFileReader(open("variants.avro", "r"), DatumReader())
        readCount = 0
        for variant in reader:
            assert variant["id"] == "myID"
            assert variant["variantSetId"] == "mySetID"
            assert len(variant["names"]) == 1
            assert "myName" in variant["names"]
            assert variant["created"] == 100L
            assert variant["updated"] is None
            assert variant["referenceName"] == "chr1"
            assert variant["start"] == 1L
            assert variant["end"] == 2L
            assert variant["referenceBases"] == "A"
            assert len(variant["alternateBases"]) == 0
            assert len(variant["info"]) == 1
            assert len(variant["info"]["infoKey"]) == 1
            assert "myInfo" in variant["info"]["infoKey"]
            assert len(variant["calls"]) == 0

            readCount += 1
        reader.close()

        # should have only read one element
        assert readCount == 1

        # remove avro file
        os.remove("variants.avro")
