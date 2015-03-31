Create data set in tests/faultydata/variants/inconsist_sampleid before 

1kgPhase1chr2_removed_HG00531.vcf.gz and its index file are soft links to tests/faultydata/variants/1kgPhase1chr2_removed_HG00531.vcf.gz

To generate 1kgPhase1chr3_removed_HG00533.vcf.gz and its index file, do the following:
```
$ cp ../../../data/variants/1kgPhase1/chr3.vcf.gz .
$ vcftools --remove-indv HG00533 --gzvcf chr3.vcf.gz --recode --out 1kgPhase1chr3_removed_HG00533
$ mv 1kgPhase1chr3_removed_HG00533.recode.vcf 1kgPhase1chr3_removed_HG00533.vcf
$ bgzip 1kgPhase1chr3_removed_HG00533.vcf
$ tabix 1kgPhase1chr3_removed_HG00533.vcf.gz
```
