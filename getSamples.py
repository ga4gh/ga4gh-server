import ga4gh.client as client
import ga4gh.protocol as protocol

c = client.HttpClient("http://localhost:8000")

BioSample = protocol.BioSample

callsets = []
readgroups = []

for d in c.searchDatasets():
    for vs in c.searchVariantSets(d.id):
        for cs in c.searchCallSets(vs.id):
            callsets.append(cs)
    for rgs in c.searchReadGroupSets(d.id):
        for rg in rgs.readGroups:
            readgroups.append(rg)

biosamples = []
for cs in callsets:
    biosamples.append(
        BioSample(
            name=cs.name,
            id=cs.name,
            description="vsid=" + str(cs.variantSetIds)))

for rg in readgroups:
    biosamples.append(
        BioSample(
            name=rg.bioSampleId,
            id=rg.bioSampleId,
            description="rgid=" + str(rg.id)))

for bs in biosamples:
    print(bs.toJsonDict())

# TODO use HTTP post to add to ref server
