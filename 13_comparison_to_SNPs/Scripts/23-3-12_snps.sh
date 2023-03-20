#!/usr/bin/env bash 

### Modified on 23-3-14 in order to get SNPs from five chromosomes for comparison 
### I decided to take chromosomes 3, 4, 8, 10, 14 (did not hear anything bad about them)
vcf1=$1
vcf2=$2
vcf3=$3
vcf4=$4
vcf5=$5
samples=$6
# subset=~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps/23-3-12_vcf_snps_set180.vcf
filtered1=~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps/23-3-14_vcf_snps_set180_encoded1.vcf.gz
filtered2=~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps/23-3-14_vcf_snps_set180_encoded2.vcf.gz
filtered3=~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps/23-3-14_vcf_snps_set180_encoded3.vcf.gz
filtered4=~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps/23-3-14_vcf_snps_set180_encoded4.vcf.gz
filtered5=~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps/23-3-14_vcf_snps_set180_encoded5.vcf.gz
shuf=~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps/23-3-14_vcf_snps_set180_sample.vcf
# bcftools view -S $samples $vcf -o $subset
bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS -S $samples $vcf1 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 's|0/1|1|g' | sed 's|1/0|1|g' | gzip -c >> $filtered1
bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS -S $samples $vcf2 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 's|0/1|1|g' | sed 's|1/0|1|g' | gzip -c >> $filtered2
bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS -S $samples $vcf3 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 's|0/1|1|g' | sed 's|1/0|1|g' | gzip -c >> $filtered3
bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS -S $samples $vcf4 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 's|0/1|1|g' | sed 's|1/0|1|g' | gzip -c >> $filtered4
bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS -S $samples $vcf5 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 's|0/1|1|g' | sed 's|1/0|1|g' | gzip -c >> $filtered5

## sbatch -c 32 --time=08:00:00 -p skylake-himem --wrap="bash 23-3-12_snps.sh ~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/biallelic3/malawi_cichlids_v3_biallelic_chr4.vcf.gz ~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/biallelic3/malawi_cichlids_v3_biallelic_chr6.vcf.gz ~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/biallelic3/malawi_cichlids_v3_biallelic_chr10.vcf.gz ~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/biallelic3/malawi_cichlids_v3_biallelic_chr12.vcf.gz ~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/biallelic3/malawi_cichlids_v3_biallelic_chr14.vcf.gz snps/23-3-1_meta_samplesids.txt"

## bash 23-3-12_snps.sh ~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/biallelic3/malawi_cichlids_v3_biallelic_chr4.vcf.gz snps/23-3-1_meta_samplesids.txt
## bash 23-3-12_snps.sh snps/23-3-12_vcf_snps_set180.vcf snps/23-3-1_meta_samplesids.txt
## bash 23-3-12_snps.sh snps/23-3-12_vcf_snps_set180.vcf snps/23-3-1_meta_samplesids.txt

### original script
# vcf=$1
# samples=$2
# subset=~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps/23-3-12_vcf_snps_set180.vcf
# filtered=~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps/23-3-14_vcf_snps_set180_encoded1.vcf
# bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS -S $samples $vcf1 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 's|0/1|1|g' | sed 's|1/0|1|g' | gzip -c >> $filtered


