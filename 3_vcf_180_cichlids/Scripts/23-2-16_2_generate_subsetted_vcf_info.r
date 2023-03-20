#!/usr/bin/env Rscript

### 2023-2-15
## In this script, I am going to try to deal with the new set of cichlids that we want to look at 
## first, subset the dataset for MEI cichlids
## then further analysis - trying to quantify insertions in each sample / samples for each insertion

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
# library(umap) ## UMAP
# library(M3C)
print("loaded libraries successfully")

getwd()
# "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara"

#################################################
#################################################
## to see how I subsetted cichlids such that they match what Richard suggested on 23-2-14 
## see script 23-2-14_standardized dataset (locally and also in this folder)

#################################################
#################################################
## read metadata
## Note that there are two metadata csv files
## 23-2-14 is additionally filtered by sublocation (arbitrary chose one sublocation which had at least 12 specimens)
## we decided that the sublocations are very close so 23-2-15 is just with top seq depth
##### NOTE this got changed AGAIN by RD on 23-2-16 (meeting at 10.30) 
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

#### subset meta for specific clades
meta_astcal <- meta[clade=="AstCal",]
meta_benthic <- meta[clade=="Benthic",]
meta_deep <- meta[clade=="Deep",]
meta_diplo <- meta[clade=="Diplo",]
meta_utaka <- meta[clade=="Utaka",]
meta_mbuna <- meta[clade=="Mbuna",]
meta_rhampo <- meta[clade=="Rhampho",]
meta_all <- rbind(meta_astcal, meta_benthic, meta_deep, meta_diplo,
meta_mbuna, meta_rhampo, meta_utaka)

#### get the ids of cichlids in each clade that you will analyse 
meta_astcal_names <- unlist(meta_astcal[,primary_id])
meta_benthic_names <- unlist(meta_benthic[,primary_id])
meta_deep_names <- unlist(meta_deep[,primary_id])
meta_diplo_names <- unlist(meta_diplo[,primary_id])
meta_utaka_names <- unlist(meta_utaka[,primary_id])
meta_mbuna_names <- unlist(meta_mbuna[,primary_id])
meta_rhampo_names <- unlist(meta_rhampo[,primary_id])
meta_all_names <- unlist(meta_all[,primary_id])

### see how the full vcf was subsetted in script 23-2-15_subset_vcf_new_meta.r
## Read the subsetted vcf
vcf_sub_MEI <- fread("/home/bw450/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-16_subset_180_cichlids_encoded.txt", 
    sep = "\t", header=T)
print(dim(vcf_sub_MEI))
print("the file should have 3522039 rows and 191 columns")

#### did not project column names smh so let's add them
col_keep <- c("Lp", "#CHROM", "POS", "ID", "REF", "INFO")
c_sums <- c("sum_0_0", "sum_1_1", "sum_0_1", "sum_0_missing", "sum_1_missing")
new_colnames <- c(col_keep, meta_ids, c_sums)
vcf_sub_MEI <- setnames(vcf_sub_MEI, colnames(vcf_sub_MEI), c(new_colnames))

## split information about transposons 
vcf_sub_MEI[, c("SVTYPE", "TE", "MEPRED", "START", "END", "LENGTH", "STRAND", "AC") := tstrsplit(INFO, ";", fixed=FALSE)]
vcf_sub_MEI[, c("delete", "SVTYPE") := tstrsplit(SVTYPE, "=", fixed=FALSE)] ## I just want the number
vcf_sub_MEI[, "delete" := NULL] ## I just want the number
vcf_sub_MEI[, c("delete", "TE") := tstrsplit(TE, "=", fixed=FALSE)] 
vcf_sub_MEI[, "delete" := NULL] 
vcf_sub_MEI[, c("delete", "START") := tstrsplit(START, "=", fixed=FALSE)] 
vcf_sub_MEI[, "delete" := NULL] 
vcf_sub_MEI[, c("delete", "END") := tstrsplit(END, "=", fixed=FALSE)]
vcf_sub_MEI[, "delete" := NULL]
vcf_sub_MEI[, c("delete", "LENGTH") := tstrsplit(LENGTH, "=", fixed=FALSE)] 
vcf_sub_MEI[, "delete" := NULL] 
vcf_sub_MEI[, c("delete", "AC") := tstrsplit(AC, "=", fixed=FALSE)] 
vcf_sub_MEI[, "delete" := NULL] 
vcf_sub_MEI[, c("order_first_call", "other") := tstrsplit(TE, "/", keep=c(1,2))] 
vcf_sub_MEI[, c("superfamily", "second_call") := tstrsplit(other, "\\|", keep=c(1,2))] 
vcf_sub_MEI[, c("order", "delete") := tstrsplit(order_first_call, "-", keep=c(1,2))] 
vcf_sub_MEI[, "delete" := NULL]
vcf_sub_MEI[, c("fam", "delete") := tstrsplit(superfamily, "-", keep=c(1,2))] 
vcf_sub_MEI[, "delete" := NULL]
vcf_sub_MEI <- vcf_sub_MEI[, order:=as.factor(order)]
vcf_sub_MEI <- vcf_sub_MEI[, fam:=as.factor(fam)]
print(dim(vcf_sub_MEI))
print("the file should have 3522039 rows and 205 columns")

## filter by quality 
vcf_sub_MEI_filter_quality <- vcf_sub_MEI[MEPRED=="MEPRED=PASS",]
print(dim(vcf_sub_MEI_filter_quality))
print("the file should have 3420880 rows and 205 columns")

nc <- c("#CHROM", "POS", "ID", "REF", "INFO", "sum_0_0", 
"sum_1_1", "sum_0_1", "sum_0_missing", "sum_1_missing", 
"SVTYPE", "TE", "MEPRED", "START", "END", "LENGTH", "STRAND", "AC", "order_first_call",
"other", "superfamily", "second_call", "order", "fam")
vcf_sub_filters_lp <- vcf_sub_MEI_filter_quality[, !nc, with=FALSE] ## numerical dataset
vcf_sub_sum_zero_lp <- apply(vcf_sub_filters_lp, 1, function(x) length(which(x=="0"))) ## vector of sums of 0

## print histogram
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-17_MEI_sum_0_in_all_clades.pdf")
hist(vcf_sub_sum_zero_lp, breaks=100) ## most are almost all 0 
dev.off()

vcf_sub_filters_sum_zero <- cbind(vcf_sub_filters_lp,vcf_sub_sum_zero_lp) ## dt num + sum zero
vcf_sub_filtered_no_zero <- vcf_sub_filters_sum_zero[vcf_sub_sum_zero_lp != length(meta_all_names),] 
print(dim(vcf_sub_filtered_no_zero))
print("the file should have 644951 rows and 182 columns")

## just to check
vcf_sub_test <- vcf_sub_filtered_no_zero[1:2,]
write.table(vcf_sub_test, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/test.txt", row.names=F, quote=F, sep="\t")
## had a look at this in excel and it looks good to me :)

## Now, I need to filter out the whole file 
lp_not_zero <- unlist(vcf_sub_filtered_no_zero[,Lp])
vcf_sub_filter_final <- vcf_sub_MEI_filter_quality[Lp %in% lp_not_zero,]
write.table(vcf_sub_filter_final, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-17_vcf_180_filtered_nonzero.txt", 
row.names=F, quote=F, sep="\t")

print(dim(vcf_sub_filter_final))
print("the final file vim_sub_filter_final has 644951 rows and 205 columns")

##### vcf_sub_filter_final is the file I will be using for all downstream analysis

