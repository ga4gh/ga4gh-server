name.fa.gz created by hand and then compressed using bgzip:

bgzip name.fa

name.fa.gz.fai and name.fa.gz.gzi created with:

samtools fadix name.fa.gz
