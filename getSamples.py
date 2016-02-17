import ga4gh.client as client
import json

c = client.HttpClient("http://localhost:8000")

callsets = []
readgroups = []

class Biosample():
    id = None
    name = None
    description = None
    def __init__(self, id, name, description):
        self.id = id
        self.name = name
        self.description = description

    def toJson(self):
        return json.dumps({
            "id": self.id,
            "name": self.name,
            "description": self.description
        })

for d in c.searchDatasets():
    for vs in c.searchVariantSets(d.id):
        for cs in c.searchCallSets(vs.id):
            callsets.append(cs)
    for rgs in c.searchReadGroupSets(d.id):
        for rg in rgs.readGroups:
            readgroups.append(rg)

biosamples = []
for cs in callsets:
    biosamples.append(Biosample(cs.name, cs.name, "vsid=" + str(cs.variantSetIds)))

for rg in readgroups:
    biosamples.append(Biosample(rg.sampleId, rg.sampleId, "rgid=" + str(rg.id)))

for bs in biosamples:
    print(bs.toJson())

# TODO use HTTP post to add to ref server

