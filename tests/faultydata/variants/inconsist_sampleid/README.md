1kgPhase1chr1.vcf.gz is a soft link to tests/data/variants/1kgPhase1/chr1.vcf.gz
Index file 1kgPhase1chr1.vcf.gz.tbi is a soft link to tests/data/variants/1kgPhase1/chr1.vcf.gz.tbi

To generate 1kgPhase1chr2_removed_HG00531.vcf.gz and its index file, do the following:
```
$ cp ../../../data/variants/1kgPhase1/chr2.vcf.gz .
$ vcftools --remove-indv HG00531 --gzvcf chr2.vcf.gz --recode --out 1kgPhase1chr2_removed_HG00531
$ mv 1kgPhase1chr2_removed_HG00531.recode.vcf 1kgPhase1chr2_removed_HG00531.vcf
$ bgzip 1kgPhase1chr2_removed_HG00531.vcf
$ tabix 1kgPhase1chr2_removed_HG00531.vcf.gz
```
