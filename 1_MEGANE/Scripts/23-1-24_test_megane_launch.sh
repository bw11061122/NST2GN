
singularity exec megane_v1.0.1.beta.sif call_genotype \
  -i /home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/alignments/Malawi_cichlids_2017/cichlid6978775/fAstCal1.2/cichlid6978775.mem.crumble.cram  \
  -fa  genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa \
  -fadb  genome/blastdb  \
  -mk  genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa.mk  \
  -rep  Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph_tab.rep \
  -repout  genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa.out \
  -repremove  genome/empty.txt  \
  -pA_ME  genome/empty.txt  \
  -mainchr  genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mainchr_new.txt  \
  -no_sex_chr \
  -p 8 \
  -lowdep \
  -skip_unmapped \
  -outdir test
