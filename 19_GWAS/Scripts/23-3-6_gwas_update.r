#### [2023-03-06]
##### Trying to maybe do this smarter
##### Repeat sums but now do this no the merged file

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
library(RColorBrewer)
library(ggpmisc)
library(HardyWeinberg)
# library(umap) ## UMAP
# library(M3C)
print("loaded libraries successfully")

getwd()
# "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara"

### see how the full vcf was subsetted in script 23-2-15_subset_vcf_new_meta.r
####### Modified 23-2-22 to include the filtered ones (good evidence for at least one insertion only)
## Read the subsetted vcf
vcf <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_ins_del_all_final.txt", 
    sep = "\t", header=T, fill=TRUE)

## factors  
vcf <- vcf[, order:=as.factor(order)]
vcf <- vcf[, family:=as.factor(family)]
vcf <- vcf[, superfamily:=as.factor(superfamily)]
upd.cols.vcf = sapply(vcf, is.factor)
vcf <- vcf[, names(vcf)[upd.cols.vcf] := lapply(.SD, factor), .SDcols = upd.cols.vcf] ## this gets rid of levels which are absent in the dt 

# trait	insID	dist2Gene	VEP	geneID	geneSymbol
# colour	chr7:41166580-41166583_O	0	intron	ENSACLG00000014805	ntrk3a #### did not find this - may be more complicated because it has an inversion 
### cichlid_all_ins_3159429
# colour	chr3:8249800-8249807_O	0	intron	ENSACLG00000018516	nlgn2a #### is in the dataset ### hAT_AC-17
# Caudal peduncle depth	chr22:13852670-13852671_O	0	intron	ENSACLG00000013131	trps1 #### unknown 19

####  GWAS update
## cichlid_all_ins_2638112 - Caudal peduncle depth	chr22:13852670-13852671_O	0	intron	ENSACLG00000013131	trps1
## cichlid_all_ins_3159429 - colour	chr3:8249800-8249807_O	0	intron	ENSACLG00000018516	nlgn2a #### is in the dataset ### hAT_AC-17

vcf_trps <- vcf[ID=="cichlid_all_ins_2638112",]
vcf_ngln <- vcf[ID=="cichlid_all_ins_3159429",]
vcf_trps_num <- vcf_trps[,c(meta_ids),with=FALSE]
vcf_ngln_num <- vcf_ngln[,c(meta_ids),with=FALSE]
meta_gwas_hits <- cbind(meta[,c("pop_id", "primary_id", "clade"), with=FALSE], t(vcf_trps_num))
setnames(meta_gwas_hits, colnames(meta_gwas_hits), c("population", "primary_id", "clade", "trps1"))
meta_gwas_hits <- cbind(meta_gwas_hits, t(vcf_ngln_num))
setnames(meta_gwas_hits, "V1", "ngln2")
write.table(meta_gwas_hits, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-7_gwas_hits_updated.txt", 
sep="\t", row.names=F)

#### What to do with these next
## get the data for all samples (2764) - check in which samples it is present? in bash 
# grep -e "cichlid_all_ins_2638112" -e "cichlid_all_ins_3159429" \
# 23-2-1_MEGANE_all_cichlids/vcf_for_phasing/cichlid_all_MEI_biallelic_no_meta_encoded.vcf \
# > 23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-7_gwas_hits_all_cichlids.txt

###############################################################
###############################################################
### doing this again (because I need to include cichlid number per population)
### meta
meta_full <- fread("data/cichlid_callset_metadata.txt", header=T)
meta_full <- meta_full[, clade:=as.factor(clade)] 
meta_full <- meta_full[, species:=as.factor(species)] 
meta_full <- meta_full[, location:=as.factor(location)]
meta_full <- meta_full[, name:=paste(genus, species, sep=" ")] 
meta_full <- meta_full[, name_loc:=paste(name, location, sep=".")]
meta_full <- meta_full[, pop_id := name_loc][(genus %in% c("Rhamphochromis", "Chilotilapia")) | is.na(name_loc), pop_id := name][]
meta_full <- meta_full[, pop_id:=as.factor(pop_id)] 

vcf_gwas <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-7_gwas_hits_all_cichlids.txt", header=F)
cichlid_names <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-6_cichlid_names.txt", header=T)
setnames(vcf_gwas, colnames(vcf_gwas), colnames(cichlid_names))
colnames(vcf_gwas) <-  sub(".mem.crumble.cram", "", colnames(vcf_gwas))
colnames(vcf_gwas) <-  sub(".mem.cram", "", colnames(vcf_gwas))

vcf_num <- vcf_gwas[,10:2678]
vcf_trans <- t(vcf_num) %>% as.data.frame() %>% rownames_to_column("primary_id")
vcf_trans <- as.data.table(vcf_trans)
meta_gwas <- merge(meta_full, vcf_trans, by="primary_id")
setnames(meta_gwas, c("V1", "V2"), c("trps1", "nlgn2"))

#### okay now we can analyse presence in each population
sum_hets_by_pop_id_trps1 <- meta_gwas[, .(sum_het_trps=length(which(trps1=="1"))), by=pop_id][order(sum_het_trps)]
sum_hom_by_pop_id_trps1 <- meta_gwas[, .(sum_hom_trps=length(which(trps1=="2"))), by=pop_id][order(sum_hom_trps)]
sum_by_pop_id_trps1 <- meta_gwas[, .(sum_trps=sum(trps1)), by=pop_id][order(sum_trps)]

sum_hets_by_pop_id_nlgn2 <- meta_gwas[, .(sum_het_nlgn=length(which(nlgn2=="1"))), by=pop_id][order(sum_het_nlgn)]
sum_hom_by_pop_id_nlgn2 <- meta_gwas[, .(sum_hom_nlgn=length(which(nlgn2=="2"))), by=pop_id][order(sum_hom_nlgn)]
sum_by_pop_id_nlgn2 <- meta_gwas[, .(sum_nlgn=sum(nlgn2)), by=pop_id][order(sum_nlgn)]

sum_gwas  <- list(sum_by_pop_id_trps1, sum_hets_by_pop_id_trps1, sum_hom_by_pop_id_trps1,
sum_by_pop_id_nlgn2, sum_hets_by_pop_id_nlgn2, sum_hom_by_pop_id_nlgn2)
merge_func <- function(...) merge(..., all = TRUE, by='pop_id')
gwas_merged <- Reduce(merge_func, sum_gwas)[order(pop_id)]
sum_pop_id <- meta_gwas[, .N, by=pop_id][order(pop_id)]

gwas_merged2 <- merge(gwas_merged, sum_pop_id, by='pop_id')[order(pop_id)]
gwas_merged2 <- gwas_merged2[,maf_trps := ((sum_het_trps + (2*sum_hom_trps)) / (2*N)) * 100]
gwas_merged2 <- gwas_merged2[,maf_nlgn := ((sum_het_nlgn + (2*sum_hom_nlgn)) / (2*N)) * 100]

gwas_merged2 <- gwas_merged2[,sum_hom_abs_trps := N - (sum_het_trps + sum_hom_trps)]
gwas_merged2 <- gwas_merged2[,sum_hom_abs_nlgn := N - (sum_het_nlgn + sum_hom_nlgn)]

write.table(gwas_merged2, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-7_gwas_by_pop_all_cichlids.txt", 
sep="\t")

### okay now try to do the HWE test
gwas_merged2_hw <- gwas_merged2[N>12 & (maf_trps>5 | maf_nlgn>5),][order(-maf_nlgn)]
gwas_hwe_nlgn <- gwas_merged2_hw[,c("pop_id", "sum_hom_abs_nlgn", "sum_het_nlgn", "sum_hom_nlgn"), with=FALSE]
setnames(gwas_hwe_nlgn, c("sum_hom_abs_nlgn", "sum_het_nlgn", "sum_hom_nlgn"), c("NN", "NM", "MM"))
HW.test.pop <- HWChisq(unlist(gwas_hwe_nlgn[9,c(2:4)]), verbose = TRUE)
## for lake Masoko this works but does not really work for other populations (too few sampples)

#### now by name (ie species not including population)
sum_hets_by_name_trps1 <- meta_gwas[, .(sum_het_trps=length(which(trps1=="1"))), by=name][order(sum_het_trps)]
sum_hom_by_name_trps1 <- meta_gwas[, .(sum_hom_trps=length(which(trps1=="2"))), by=name][order(sum_hom_trps)]
sum_by_name_trps1 <- meta_gwas[, .(sum_trps=sum(trps1)), by=name][order(sum_trps)]

sum_hets_by_name_nlgn2 <- meta_gwas[, .(sum_het_nlgn=length(which(nlgn2=="1"))), by=name][order(sum_het_nlgn)]
sum_hom_by_name_nlgn2 <- meta_gwas[, .(sum_hom_nlgn=length(which(nlgn2=="2"))), by=name][order(sum_hom_nlgn)]
sum_by_name_nlgn2 <- meta_gwas[, .(sum_nlgn=sum(nlgn2)), by=name][order(sum_nlgn)]

sum_gwas_name <- list(sum_by_name_trps1, sum_hets_by_name_trps1, sum_hom_by_name_trps1,
sum_by_name_nlgn2, sum_hets_by_name_nlgn2, sum_hom_by_name_nlgn2)
merge_func_name <- function(...) merge(..., all = TRUE, by='name')
gwas_merged_name <- Reduce(merge_func_name, sum_gwas_name)[order(name)]
sum_name_name <- meta_gwas[, .N, by=name][order(name)]

gwas_merged2_name <- merge(gwas_merged_name, sum_name_name, by='name')[order(name)]
gwas_merged2_name <- gwas_merged2_name[,maf_trps := ((sum_het_trps + (2*sum_hom_trps)) / (2*N)) * 100]
gwas_merged2_name <- gwas_merged2_name[,maf_nlgn := ((sum_het_nlgn + (2*sum_hom_nlgn)) / (2*N)) * 100]

gwas_merged2_name <- gwas_merged2_name[,sum_hom_abs_trps := N - (sum_het_trps + sum_hom_trps)]
gwas_merged2_name <- gwas_merged2_name[,sum_hom_abs_nlgn := N - (sum_het_nlgn + sum_hom_nlgn)]

write.table(gwas_merged2_name, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-7_gwas_by_pop_all_cichlids.txt", 
sep="\t")

### okay now try to do the HWE test
gwas_merged2_hw_name <- gwas_merged2_name[N>12 & (maf_trps>5 | maf_nlgn>5),][order(-maf_nlgn)]
gwas_hwe_nlgn_name <- gwas_merged2_hw_name[,c("name", "sum_hom_abs_nlgn", "sum_het_nlgn", "sum_hom_nlgn"), with=FALSE]
setnames(gwas_hwe_nlgn_name, c("sum_hom_abs_nlgn", "sum_het_nlgn", "sum_hom_nlgn"), c("NN", "NM", "MM"))
HW.test.pop_name.1 <- HWChisq(unlist(gwas_hwe_nlgn_name[1,c(2:4)]), cc = 0, verbose = TRUE)
HW.test.pop_name.2 <- HWChisq(unlist(gwas_hwe_nlgn_name[2,c(2:4)]), verbose = TRUE)
HW.test.pop_name.3 <- HWChisq(unlist(gwas_hwe_nlgn_name[3,c(2:4)]), verbose = TRUE)
HW.test.pop_name.4 <- HWChisq(unlist(gwas_hwe_nlgn_name[4,c(2:4)]), verbose = TRUE)
HW.test.pop_name.5 <- HWChisq(unlist(gwas_hwe_nlgn_name[5,c(2:4)]), verbose = TRUE)
HW.test.pop_name.6 <- HWChisq(unlist(gwas_hwe_nlgn_name[6,c(2:4)]), verbose = TRUE)
HW.test.pop_name.7 <- HWChisq(unlist(gwas_hwe_nlgn_name[7,c(2:4)]), verbose = TRUE)
HW.test.pop_name.8 <- HWChisq(unlist(gwas_hwe_nlgn_name[8,c(2:4)]), verbose = TRUE)
HW.test.pop_name.9 <- HWChisq(unlist(gwas_hwe_nlgn_name[9,c(2:4)]), verbose = TRUE)
HW.test.pop_name.10 <- HWChisq(unlist(gwas_hwe_nlgn_name[10,c(2:4)]), verbose = TRUE)
### 1, 3, 10 are significant 

x <- c(NN=300, MN=100, MM=0)
test <- HWChisq(x, verbose = TRUE)
## Permutation test from here: https://cran.r-project.org/web/packages/HardyWeinberg/vignettes/HardyWeinberg.pdf
# Hardy-Weinberg equilibrium refers to the statistical independence of alleles within individuals. 
# This independence can also be assessed by a permutation test, where all 2n alleles of all individuals 
# are written out as a single sequence (E.g. AAAAABABBBAA....). This sequence is then permuted many times, 
# and for each permuted sequence pairs of successive alleles are taken as individuals. 
# For each permutation a test statistic (the pseudo-statistic) for disequilibrium is computed. 
# The test statistic for the original observed sample is compared against the distribution of the 
# pseudo-statistic, where the latter was generated under the null hypothesis. 
# The p value of the test is calculated as the fraction of permuted samples for which the pseudo-statistic 
# is equal to or exceeds the test statistic. Such a test is computer intensive but has the advantage 
# that it does not rely on asymptotic assumptions. Function HWPerm performs this test.

HW.test.perm_name.1 <- HWPerm(unlist(gwas_hwe_nlgn_name[1,c(2:4)]), verbose = TRUE) # 1
HW.test.perm_name.2 <- HWPerm(unlist(gwas_hwe_nlgn_name[2,c(2:4)]), verbose = TRUE) # 1
HW.test.perm_name.3 <- HWPerm(unlist(gwas_hwe_nlgn_name[3,c(2:4)]), verbose = TRUE) # 0
HW.test.perm_name.4 <- HWPerm(unlist(gwas_hwe_nlgn_name[4,c(2:4)]), verbose = TRUE) # 0.122
HW.test.perm_name.5 <- HWPerm(unlist(gwas_hwe_nlgn_name[5,c(2:4)]), verbose = TRUE) # 1
HW.test.perm_name.6 <- HWPerm(unlist(gwas_hwe_nlgn_name[6,c(2:4)]), verbose = TRUE) # 1 
HW.test.perm_name.7 <- HWPerm(unlist(gwas_hwe_nlgn_name[7,c(2:4)]), verbose = TRUE) # 0.1480588
HW.test.perm_name.8 <- HWPerm(unlist(gwas_hwe_nlgn_name[8,c(2:4)]), verbose = TRUE) # 0.3437647 
HW.test.perm_name.9 <- HWPerm(unlist(gwas_hwe_nlgn_name[9,c(2:4)]), verbose = TRUE) # 1
HW.test.perm_name.10 <- HWPerm(unlist(gwas_hwe_nlgn_name[10,c(2:4)]), verbose = TRUE) # 0.006

##### okay so there is evidence for not following HW equilibrium in A. calliptera
##### Still, I have very high frequencies of alleles in Otopharynx and C. chrysonotus and C. virginalis
supp_id_oto <- unlist(meta_full[name=="Otopharynx argyrosoma",supplier_id])
supp_id_cch <- unlist(meta_full[name=="Copadichromis chrysonotus",supplier_id])
supp_id_ccv <- unlist(meta_full[name=="Copadichromis virginalis",supplier_id])

write.table(supp_id_oto, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-7_supp_id_otopharynx.txt", 
sep="\t")
write.table(supp_id_cch, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-7_supp_id_cop_chrysonotus.txt", 
sep="\t")
write.table(supp_id_ccv, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-7_supp_id_cop_virginalis.txt", 
sep="\t")

## where to find pictures: ~/rds/rds-8b3VcZwY7rY/projects/cichlid/photos
# ls -d abc*   # list all files starting with abc---
# ls -d *abc*  # list all files containing --abc--
# ls -d *abc   # list all files ending with --abc



