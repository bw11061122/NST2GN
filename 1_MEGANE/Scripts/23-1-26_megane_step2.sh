#!/usr/bin/env bash

# Megane step 2 - for the original code, 
# see: https://github.com/shohei-kojima/MEGAnE/wiki/Run-with-a-custom-repeat-library-and-annotation-(mostly-for-analysis-of-non-human-WGS)

## 23-1-26
## I am modifying the Megane step 2 script to make automation easier
## I will include the identifier for each batch (which I used to create results_$ID folder in step 1)
## so that I can submit the same script for each batch

sif=megane_v1.0.1.beta.sif # path to Megane in my directory

# Path to output directories from Step 1) I am merging
ls -d results_$1/* > dirlist_$1.txt
## -d option is to specify directories w/in the folder
## this is useful if you need help with ls: https://www.rapidtables.com/code/linux/ls.html

# merge non-reference ME insertions
### non-reference refers to insertions which are not present in the reference genome
singularity exec ${sif} joint_calling_hs \
-merge_mei \
-f dirlist_$1.txt \
-fa genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa \
-rep Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph_tab.rep \
-no_sex_chr \
-cohort_name cichlid_$1

# merge reference ME polymorphisms
### reference here refers to "insertions" which ARE in the reference genome
### so we can understand this as deletions in the non-reference samples compared to the reference genome
### this is why it says "merge_absent_me"
singularity exec ${sif} joint_calling_hs \
-merge_absent_me \
-f dirlist_$1.txt \
-fa genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa \
-rep Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph_tab.rep \
-no_sex_chr \
-cohort_name cichlid_$1

## output is located in the directory jointcall_out

## I am using the "-no_sex_chr" option because cichlids do not have a sex chromosome
## what is cohort name?
## are our genomes sequenced on the same platform? (they do not recommend making a joint call on WGS sequences w/ different methods)

## to run: bash 23-1-26_megane_step2.sh test1_lf  ### $1 is the identifier I used to create results folder


