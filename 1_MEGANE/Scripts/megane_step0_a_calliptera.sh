#!/usr/bin/env bash

# Step 0 makes the k-m,er reference file from the reference genome 
## All species we are going to analyse are going to be mapped to A. calliptera
## Therefore, I am only running the script once - it would be good to understand why we are only using one script
## Not quite sure where exactly I should output this?

# first need to generate a fa.fai file 
samtools faidx \
/home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa \
-o /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa.fai

sbatch -c 16 -t 36:00:00 \
--wrap="singularity exec /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/megane_v1.0.1.beta.sif build_kmerset \
-fa /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa \
-outdir /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/" -J name

# to run 