"""
    Convert thousand genomes samples csv (converted from excel)
    into json documents
    Headers:
    Sample,Family ID,Population,Population Description,Gender,
    Relationship,Unexpected Parent/Child ,Non Paternity,
    Siblings,Grandparents,Avuncular,Half Siblings,
    Unknown Second Order,Third Order,Other Comments
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import csv
import json
import datetime
import utils


@utils.Timed()
def main():
    with open('20130606_sample_info.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            description = "{} {} {}".format(
                row['Population'],
                row['Population Description'],
                row['Gender'])
            info = {}
            for key in row:
                info[key] = [row[key]]
            biosample = {
                "name": row['Sample'],
                "description": description,
                "disease": None,  # Ontology term
                "createDateTime": datetime.datetime.now().isoformat(),
                "updateDateTime": datetime.datetime.now().isoformat(),
                "info": info
            }
            if row['Gender'] == 'male':
                sex = {
                    "id": "PATO:0020001",
                    "term": "male genotypic sex",
                    "sourceName": "PATO",
                    "sourceVersion": "2015-11-18"
                }
            elif row['Gender'] == 'female':
                sex = {
                    "id": "PATO:0020002",
                    "term": "female",
                    "sourceName": "PATO",
                    "sourceVersion": "2015-11-18"
                }
            else:
                sex = None
            individual = {
                "name": row['Sample'],
                "description": description,
                "species": {
                    "term": "Homo sapiens",
                    "id": "NCBITaxon:9606",
                    "sourceName": "http://purl.obolibrary.org/obo",
                    "sourceVersion": "2016-02-02"
                },
                "sex": sex,
                "createDateTime": datetime.datetime.now().isoformat(),
                "updateDateTime": datetime.datetime.now().isoformat(),
                "info": info
            }
            bioSampleFilepath = \
                "biodata/biosamples/" + row['Sample'] + ".json"
            individualFilePath = \
                "biodata/individuals/" + row['Sample'] + ".json"
            with open(bioSampleFilepath, 'w') as bioSampleOut:
                json.dump(biosample, bioSampleOut)

            with open(individualFilePath, 'w') as individualOut:
                json.dump(individual, individualOut)

if __name__ == "__main__":
    main()
