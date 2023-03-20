#!/usr/bin/env Rscript

### 2023-2-24
### continued 23-2-27
## What I am trying to do is to get a new vcf file 
## Which will contain both insertions and deletions
## which are polymorphic across samples
## and where deletions are encoded such that 0/0 is 1/1; 0/missing is 1/missing, 1/missing is 0/missing and 1/1 is 0/0
## this way, everything is set such that 0/0 is the reference and 1/1 is insertion not present in the reference 
## this should also remove bias from the dataset (which is due to the insertiong being called on the reference)

## libraries
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

## read meta
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

## read the vcf files
vcf_del <- fread("23-2-1_MEGANE_all_cichlids/23-2-16/23-2-27_vcf_subsetted_deletions_not_encoded.txt",
    sep="\t", header=T)
info_col <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
colnames(vcf_del) <- gsub(".mem.crumble.cram", "",colnames(vcf_del))
vcf_del <- vcf_del[,c(info_col, meta_ids), with=FALSE] 
dim(vcf_del)
#### okay all done!

## adjust column names 
vcf_del[, c("SVTYPE", "TE", "START", "END", "LENGTH", "AC") := tstrsplit(INFO, ";", fixed=FALSE)]
vcf_del[, c("delete", "SVTYPE") := tstrsplit(SVTYPE, "=", fixed=FALSE)] ## I just want the number
vcf_del[, "delete" := NULL] ## I just want the number
vcf_del[, c("delete", "TE") := tstrsplit(TE, "=", fixed=FALSE)] 
vcf_del[, "delete" := NULL] 
vcf_del[, c("delete", "START") := tstrsplit(START, "=", fixed=FALSE)] 
vcf_del[, "delete" := NULL] 
vcf_del[, c("delete", "END") := tstrsplit(END, "=", fixed=FALSE)]
vcf_del[, "delete" := NULL]
vcf_del[, c("delete", "LENGTH") := tstrsplit(LENGTH, "=", fixed=FALSE)] 
vcf_del[, "delete" := NULL] 
vcf_del[, c("delete", "AC") := tstrsplit(AC, "=", fixed=FALSE)] 
vcf_del[, "delete" := NULL] 
vcf_del[, c("order_first_call", "other") := tstrsplit(TE, "/", keep=c(1,2))] 
vcf_del[, c("order", "delete") := tstrsplit(order_first_call, "-", keep=c(1,2))] ### clean up the orders
vcf_del[, "delete" := NULL]
vcf_del[, c("superfamily_info", "add") := tstrsplit(other, "\\|", keep=c(1,2))] 
vcf_del[, c("family", "delete") := tstrsplit(superfamily_info, "_AstCal")] ### only leaves the superfamily column
vcf_del[, c("superfamily", "add") := tstrsplit(family, "-", keep=c(1,2))]
vcf_del <- vcf_del[, order:=as.factor(order)]
vcf_del <- vcf_del[, family:=as.factor(family)]
vcf_del <- vcf_del[, superfamily:=as.factor(superfamily)]
print(dim(vcf_del))

### do encoding already keeping in mind that 0/0 is 2 and 1/1 is 0
write.table(vcf_del, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/26-2-27_vcf_subsetted_deletions_not_encoded.txt", 
row.names=F, quote=F, sep="\t")

#### This is in bash
#### the encoding is such that 0/0 == 2, 1/1 == 0 (this makes it equivalent to insertions) etc
"sed -e 's#0/0#2#g' -e 's#0/1#1#g' -e 's#1/1#0#g' -e 's#0/.#1.5#g' -e 's#./0#1.5#g' -e 's#./1#0.5#g' -e 's#./1#0.5#g' \
< 23-2-1_MEGANE_all_cichlids/23-2-16_subset/26-2-27_vcf_subsetted_deletions_not_encoded.txt > \
23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_subsetted_deletions_encoded.txt"

######## Now checking if everything is okay with the file (ie encoding worked)
vcf_del_encoded <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_subsetted_deletions_encoded.txt",
    sep="\t", header=T)
vcf_del <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_subsetted_deletions_not_encoded.txt",
    sep="\t", header=T)
head(vcf_del)
head(vcf_del_encoded)
## so all good 

vcf_ins <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-22_vcf_filtered_evidence_for_ins.txt", 
    sep = "\t", header=T)
vcf_del_filter <- vcf_del_encoded[FILTER=="PASS",]
dim(vcf_del_filter) # 73024 row, 201 columns
vcf_ins <- vcf_ins[,FILTER := "PASS"]

#### keep these columns in the vcf for both
col_info <- c("#CHROM", "POS","ID", "REF", "FILTER", "INFO", "SVTYPE", "TE", "START", "END",
"LENGTH", "AC", "order_first_call", "other")
vcf_ins_sub <- vcf_ins[,c(col_info, meta_ids), with=FALSE]
vcf_del_sub <- vcf_del_filter[,c(col_info, meta_ids), with=FALSE]

###### okay so both vcf files have the right fish, right metadata and they are filtered 
###### such that only the insertions with filter=="PASS" are present
###### Now filtering: I want to remove the insertions which are 0 from the deletion file 
vcf_del_num <- vcf_del_sub[,meta_ids, with=FALSE]
del_sum_deletions <- apply(vcf_del_num, 1, function(x) length(which(x=="2")))
vcf_del_sub <- cbind(vcf_del_sub, del_sum_deletions)
vcf_del_sub_filt <- vcf_del_sub[del_sum_deletions != 180,]
dim(vcf_del_sub_filt)
vcf_del_sub <- vcf_del_sub_filt[,del_sum_deletions := NULL]

vcf_ins_sub <- vcf_ins_sub[,type := "MEI"]
vcf_del_sub <- vcf_del_sub[,type := "MEA"]

vcf_final <- rbind(vcf_ins_sub, vcf_del_sub)

##### Now sort out metadata 
vcf_final[, c("delete", "TE") := tstrsplit(TE, "=", fixed=FALSE)] 
vcf_final[, "delete" := NULL] 
vcf_final[, c("order", "other") := tstrsplit(TE, "/", keep=c(1,2))] 
vcf_final[, c("order_info", "delete") := tstrsplit(order, "-", keep=c(1,2))] ### clean up the orders
vcf_final[, "delete" := NULL]
vcf_final[, c("superfamily_info", "add") := tstrsplit(other, "\\|", keep=c(1,2))] 
vcf_final[, c("family", "delete") := tstrsplit(superfamily_info, "_AstCal")] ### only leaves the superfamily column
vcf_final[, c("superfamily", "add") := tstrsplit(superfamily_info, "-", keep=c(1,2))]
vcf_final[, superfamily_info := NULL]
vcf_final[, delete := NULL]
vcf_final[, add := NULL]
vcf_final[, order := NULL]
vcf_final <- vcf_final[, order_info:=as.factor(order_info)]
vcf_final <- vcf_final[, family:=as.factor(family)]
vcf_final <- vcf_final[, superfamily:=as.factor(superfamily)]
setnames(vcf_final, "order_info", "order")
print(dim(vcf_final))

write.table(vcf_final, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_ins_del_all_final.txt", 
row.names=F, quote=F, sep="\t") #### I don't know how to describe it to be fair 

