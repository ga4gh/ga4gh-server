two.fa.gz created by hand and then compressed using bgzip:

bgzip two.fa

two.fa.gz.fai and two.fa.gz.gzi created with:

samtools fadix two.fa.gz
