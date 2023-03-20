#!/usr/bin/env Rscript

### 2023-2-21
## In this script, I want to get a subsetted vcf file to be able to work with it interactively 

## load required libraries 
.libPaths( c( "~/R/x86_64-redhat-linux-gnu-library/4.2" , .libPaths() ) ) 
.libPaths()
library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library("FactoMineR") ## PCA
library("factoextra") ## PCA 
print("loaded libraries successfully")

getwd()
# "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/23-2-1_MEGANE_all_cichlids"

## read metadata
## Note that there are two metadata csv files
meta <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-16_cichlid_meta_standardized_dataset.csv",
              sep=",", header = T) ## this is the full metadata
meta <- meta[, clade:=as.factor(clade)] 
meta <- meta[, sex:=as.factor(sex)] 
meta <- meta[, location:=as.factor(location)] 
meta <- meta[, sublocation:=as.factor(sublocation)] 
meta <- meta[, species:=as.factor(species)] 
meta <- meta[, genus:=as.factor(genus)] 
meta <- meta[, names:=as.factor(names)] 
meta <- meta[, name_loc:=as.factor(name_loc)] 
meta_ids <- unlist(meta[,primary_id]) 
print("loaded metadata")

### read the whole MEI file 
vcf_all_MEI <- fread("23-2-1_MEGANE_all_cichlids/vcf_for_phasing/cichlid_all_MEI_biallelic_no_meta.vcf",
                     sep="\t", header = T)
print("loaded complete MEI vcf file - encoded, no metadata")
## I need to change colnames to be able to filer
colnames(vcf_all_MEI) <- gsub(".mem.crumble.cram", "", colnames(vcf_all_MEI))
colnames(vcf_all_MEI) <- gsub(".mem.cram", "", colnames(vcf_all_MEI))

## subset for the cichlids you are interested in
col_keep <- c("#CHROM", "POS", "ID", "REF", "INFO")
vcf_sub_MEI <- vcf_all_MEI[,c(col_keep, meta_ids), with=FALSE]
print(dim(vcf_sub_MEI)) 
print("subsetted vcf file to 180 cichlid samples of interest")

vcf_sub_MEI[,sum_0_0 := apply(vcf_sub_MEI, 1, function(x) length(which(x=="0/0")))]
vcf_sub_MEI[,sum_0_1 := apply(vcf_sub_MEI, 1, function(x) length(which(x=="0/1")))]
vcf_sub_MEI[,sum_1_1 := apply(vcf_sub_MEI, 1, function(x) length(which(x=="1/1")))]
vcf_sub_MEI[,sum_0_missing := apply(vcf_sub_MEI, 1, function(x) length(which(x=="0/.")))]
vcf_sub_MEI[,sum_1_missing := apply(vcf_sub_MEI, 1, function(x) length(which(x=="./1")))]

vcf_sub_MEI <- vcf_sub_MEI[sum_1_1 > 0 | sum_1_missing > 0| sum_1_1 > 0,]

for(col in names(vcf_sub_MEI)) set(vcf_sub_MEI, i=which(vcf_sub_MEI[[col]]=="0/0"), j=col, value="0")
for(col in names(vcf_sub_MEI)) set(vcf_sub_MEI, i=which(vcf_sub_MEI[[col]]=="0/1"), j=col, value="1")
for(col in names(vcf_sub_MEI)) set(vcf_sub_MEI, i=which(vcf_sub_MEI[[col]]=="1/1"), j=col, value="2")
for(col in names(vcf_sub_MEI)) set(vcf_sub_MEI, i=which(vcf_sub_MEI[[col]]=="0/."), j=col, value="0.5")
for(col in names(vcf_sub_MEI)) set(vcf_sub_MEI, i=which(vcf_sub_MEI[[col]]=="./1"), j=col, value="1.5")

print("MEI subset columns encoded")
write.table(vcf_sub_MEI, file="/home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-22_vcf_filtered_evidence_for_ins.txt", 
            sep = "\t", col.names=T)
print("wrote encoded file")