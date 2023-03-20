#!/usr/bin/env bash

# Creating files required for step 1 in megane
# Repeat library and repeat annotation (via RepeatMasker) are already done

## files that are required:
# repeat library - from Pio
## this is /home/bw450/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph.rep
# repeat annotation - this requires repeat masker output - already present 
# list of ME repeats with poly-A tails - human-specific (Alu), not using
# list of non-ME repeats - not using (just in case some are important)
# list of main chromosomes - I will a file use with scaffolds
# BLAST DB of reference genomes - done 23/1/23
## see here if in doubt: https://github.com/shohei-kojima/MEGAnE/wiki/Things-to-check-before-analysis
## see here if in doubt: https://github.com/shohei-kojima/MEGAnE/wiki/Run-with-a-custom-repeat-library-and-annotation-(mostly-for-analysis-of-non-human-WGS) 

# note that cichlids do not have a sex chromosome so we are using -no_sex_chr option
# -p I think is the number of cores
# -skip_unmapped: "on-human WGS often contains a far larger number of unmapped reads than human. In such cases, we recommend to use the -skip_unmapped option"

### List of polyA - this can be ommitted, as it is a human-specific thing
# would be good to check if there are TEs that carry poly-A in cichlids

# Generate blastdb
singularity exec /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/megane_v1.0.1.beta.sif /usr/local/bin/ncbi-blast-2.12.0+/bin/makeblastdb \
-in /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa \
-out /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/blastdb \
-dbtype nucl \
-parse_seqids

## This was to test if the RepBase file is okay
singularity exec /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/megane_v1.0.1.beta.sif /usr/local/bin/ncbi-blast-2.12.0+/bin/makeblastdb \
-in /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph.rep \
-out /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/Lib3_blastdb \
-dbtype nucl \
-parse_seqids

# Generate chromosome names file
## Just the first column and include scaffolds 
cat /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa | awk '{print $1}' | grep ">" \
> /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mainchr.txt
sed 's/^\(>\)*//' /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mainchr.txt \
> /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mainchr_new.txt

# List of non-ME repeats
## Unknown could be something that is a TE but the classification is not great so we just do not know - maybe it is something
## Pio suggests not removing anything just in case - so I will NOT run this script but will keep it so that I know how to do it
cat /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph.rep \
| grep '>' | cut -f 2 | sort | uniq
cat << _EOT_ > /home/bw450/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_TE_remove.txt 
rRNA
tRNA-Core-RTE
_EOT_

cat /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph.rep \
| grep '>' | cut -f 2 | sort | uniq
cat << _EOT_ > /home/bw450/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_polyA.txt 
L1
L1-Tx1
_EOT_

## In case I am wondering, 2007 Nat Gen Rev says "poly(A) for LINEs is proposed but unconfirmed"
## But this paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4671821/ shows that poly(A) is required for LINE-1 retrotransposition


