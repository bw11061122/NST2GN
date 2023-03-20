### 23-3-3 update

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
library(vcfR)
# library(umap) ## UMAP
# library(M3C)
print("loaded libraries successfully")

getwd()
# "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara"

#################################################
#################################################
## read metadata
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
meta_ids <- unlist(meta[,primary_id]) ### this has all the IDs that you are interested in
meta <- meta[, pop_id := name_loc][(genus %in% c("Rhamphochromis", "Chilotilapia")) | is.na(name_loc), pop_id := names][]
meta <- meta[, pop_id:=as.factor(pop_id)] 
upd.cols = sapply(meta, is.factor)
meta <- meta[, names(meta)[upd.cols] := lapply(.SD, factor), .SDcols = upd.cols] ## this gets rid of levels which are absent in the dt 
length(meta_ids) ### expecting 180 samples to be identified 
print("loaded metadata")

#### I run this interactively and it worked I think
# snps <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-3_vcf_snps_set180.sample.full.vcf", sep = "\t", header=F)
# info_col <- c("chromosome", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
# setnames(snps, colnames(snps), c(info_col, meta_ids))

### do a PCA on this again
# snps_num <- snps[,c(meta_ids), with=FALSE]
# new_ids <- colnames(snps_num)
# snps_num <- snps_num[, c(new_ids) := lapply(.SD, function(x)substring(x,1, 3)), .SDcols=ids]
# snps2 <- snps_num[,c(new_ids), with=FALSE]
# snps_info <- snps[,c(info_col), with=FALSE]
# snps_gt <- cbind(snps_info, snps2)
# write.table(snps_gt, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-3_vcf_snps_set180.genotype.full.vcf", row.names=F, quote=F)

#### Get this to run in R
# sbatch -c 8 --time=02:00:00 --wrap="Rscript 23-3-2_SNPs_PCA.r"

### this goes in bash
# "sed -e 's#0/0#0#g' -e 's#./.#0.5#g' -e 's#0/1#1#g' -e 's#1/1#2#g' -e 's#0/.#0.5#g' -e 's#./0#0.5#g' -e 's#./1#1.5#g' -e 's#./1#1.5#g' \
# < 23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-3_vcf_snps_set180.genotype.full.vcf > \
# 23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-3_vcf_snps_set180.genotype.encoded.full.vcf"

##### Back to R :) 
vcf <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-1_vcf_snps_set180.genotype.encoded.full.vcf", header=T, sep=" ", fill=TRUE)
print(dim(vcf))
vcf_num <- vcf[,c(meta_ids), with=FALSE]
print(dim(vcf_num))
vcf_trans <- transpose(vcf_num)
vcf_trans_full <- cbind(vcf_trans, meta)
res.pca <- PCA(vcf_trans,  graph = FALSE) # PCA on JUST numerical
PCA_MEI_vcf_clade <- fviz_pca_ind(res.pca,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full$clade), # color by location
                                     title="PCA - 1000 SNPs"
                                     # addEllipses = TRUE # Concentration ellipses
)

PCA_MEI_vcf_pop <- fviz_pca_ind(res.pca,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full$pop_id), # color by location
                                     title="PCA - 1000 SNPs"
                                     # addEllipses = TRUE # Concentration ellipses
)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_PCA_SNPs_pop_id_full.pdf")
PCA_MEI_vcf_pop
dev.off()

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_PCA_SNPs_pop_clade_full.pdf")
PCA_MEI_vcf_clade
dev.off()

#### sbatch -c 32 --time=24:00:00 -p skylake-himem --wrap="Rscript 23-3-3_snps_pca_update.r"
