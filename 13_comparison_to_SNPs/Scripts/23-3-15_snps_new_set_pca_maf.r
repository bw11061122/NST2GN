#### [2023-3-15]
##### SNP file analysis 
##### This time using SNPs datasets from 5 chromosomes

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
library(R.utils) ## to read gzipped files
# library(umap) ## UMAP
# library(M3C)
print("loaded libraries successfully")

getwd()
# "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara"

## plotting
cbPalette <- c("#EF5350", "#2196F3", "#43A047", "#9575CD", "#BDBDBD", "#D1C4E9", "#757575", "#795548") # orders
popPalette <- c("#9CCC65", "#7CB342", "#7E57C2", "512DA8", "#AB47BC", "#6A1B9A", "#1B5E20", "8BC34A", "388E3C",
"#B71C1C", "#E57373", "#D32F2F", "#1A237E", "#BF360C", "#827717")
cladePalette <- c("#AED581", "#5E35B1", "#1B5E20", "#B71C1C", "#1A237E", "#EF6C00", "#827717")
shapes <- c(3, 4, 8, 15, 16, 17, 18) ## chose to not have "empty" shapes

## read metadata
meta <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-16_cichlid_meta_standardized_dataset.csv",
                              sep=",", header = T) ## this is the full metadata
meta <- meta[, clade:=as.factor(clade)] 
meta <- meta[, sex:=as.factor(sex)] 
meta <- meta[, location:=as.factor(location)] 
meta <- meta[,sublocation:=as.factor(sublocation)] 
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

levels.pop <- c("Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
"Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
"Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
"Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
"Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", "Diplotaxodon limnothrissa.Southwest_arm")
col <- "pop_id"
meta[, (col) := factor(get(col), levels = levels.pop)]

### Read the vcf I created on 23-3-1 (23-3-1_compare_to_vcf.sh) 
vcf1 <- fread("snps/23-3-14_vcf_snps_set180_encoded1.vcf.gz", header=F, sep="\t")
vcf2 <- fread("snps/23-3-14_vcf_snps_set180_encoded2.vcf.gz", header=F, sep="\t")
vcf3 <- fread("snps/23-3-14_vcf_snps_set180_encoded3.vcf.gz", header=F, sep="\t")
vcf4 <- fread("snps/23-3-14_vcf_snps_set180_encoded4.vcf.gz", header=F, sep="\t")
vcf5 <- fread("snps/23-3-14_vcf_snps_set180_encoded5.vcf.gz", header=F, sep="\t")

setnames(vcf1, colnames(vcf1), c("chr", "coord", "v1", "v2", meta_ids))
setnames(vcf2, colnames(vcf2), c("chr", "coord", "v1", "v2", meta_ids))
setnames(vcf3, colnames(vcf3), c("chr", "coord", "v1", "v2", meta_ids))
setnames(vcf4, colnames(vcf4), c("chr", "coord", "v1", "v2", meta_ids))
setnames(vcf5, colnames(vcf5), c("chr", "coord", "v1", "v2", meta_ids))

vcf_nz1 <- vcf1[ rowSums(vcf1[,5:184]) > 0, ] ## polymorphic
vcf_nz2 <- vcf2[ rowSums(vcf2[,5:184]) > 0, ] ## polymorphic
vcf_nz3 <- vcf3[ rowSums(vcf3[,5:184]) > 0, ] ## polymorphic
vcf_nz4 <- vcf4[ rowSums(vcf4[,5:184]) > 0, ] ## polymorphic
vcf_nz5 <- vcf5[ rowSums(vcf5[,5:184]) > 0, ] ## polymorphic

### sample 517,436
set.seed(1)
vcf_sub1 <- vcf_nz1[sample(.N, 103487)]
vcf_sub_num1 <- vcf_sub1[,c(meta_ids), with=FALSE]
vcf_sub2 <- vcf_nz2[sample(.N, 103487)]
vcf_sub_num2 <- vcf_sub2[,c(meta_ids), with=FALSE]
vcf_sub3 <- vcf_nz3[sample(.N, 103487)]
vcf_sub_num3 <- vcf_sub3[,c(meta_ids), with=FALSE]
vcf_sub4 <- vcf_nz4[sample(.N, 103487)]
vcf_sub_num4 <- vcf_sub4[,c(meta_ids), with=FALSE]
vcf_sub5 <- vcf_nz5[sample(.N, 103487)]
vcf_sub_num5 <- vcf_sub5[,c(meta_ids), with=FALSE]

### summed per sampe
vcf_snps <- rbind(vcf_sub_num1, vcf_sub_num2, vcf_sub_num3, vcf_sub_num4, vcf_sub_num5)
write.table(vcf_snps, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-15_snps_subset_for_comparison.txt",
sep="\t", row.names=F, quote=F)

vcf_trans <- transpose(vcf_snps)
vcf_trans_full <- cbind(vcf_trans, meta)
res.pca <- PCA(vcf_trans,  graph = FALSE) # PCA on JUST numerical

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-16_PCA_snps_scree_v3.pdf")
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()

basic_plot_snp <- fviz_pca_ind(res.pca, label="none")
levels.clade <- c("AstCal", "Benthic", "Deep", "Diplo", "Mbuna","Rhampho", "Utaka")
col_clade <- "clade"
vcf_trans_full[, (col_clade) := factor(get(col_clade), levels = levels.clade)]
up_plot_snp <- ggplot(cbind(basic_plot_snp$data,vcf_trans_full[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=4) + 
    theme_bw(base_size=18) +
    labs(title="PCA - polymorphic SNPs",x = "PC1 (2.7%)", y = "PC2 (2.3%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=popPalette) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-15_up_plot_snps.pdf", width=16, height=8)
up_plot_snp
dev.off()

basic_plot_snp_34 <- fviz_pca_ind(res.pca, label="none", axes=c(3,4))
up_plot_snp_34 <- ggplot(cbind(basic_plot_snp_34$data,vcf_trans_full[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=4) + 
    theme_bw(base_size=18) +
    labs(title="PCA - polymorphic SNPs",x = "PC3 (2%)", y = "PC4 (1.9%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=popPalette) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-15_up_plot_snp_34.pdf", width=16, height=8)
up_plot_snp_34
dev.off()

#### comparison to minor allele frequency
sum_ins_all <- vcf_snps[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ids]
vcf_maf <- as.data.table(sum_ins_all)
vcf_maf <- vcf_maf[,total_alleles := 360] ## 2 * 180 samples
vcf_maf <- vcf_maf[,maf := sum_ins_all / total_alleles]
vcf_maf <- vcf_maf[order(maf)]

maf_all_snp <-ggplot(vcf_maf, aes(x=maf)) + 
  geom_histogram(binwidth=0.01)+
  geom_vline(xintercept=c(0.05), color="red")+
  labs(title="Minor allele frequency: SNPs",x = "Minor allele frequency", y = "Frequency")+
  theme_bw(base_size = 26)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-15_maf_snps_v2.pdf", width=12, height=12)
maf_all_snp
dev.off()

vcf_maf_05 <- vcf_maf[,.(maf_05=length(which(maf>0.05)) / length(maf))][order(-maf_05)]

### Rscript 23-3-15_snps_new_set_pca_maf.r

### new samples
### sample 517,436
set.seed(5)
vcf2_sub1 <- vcf_nz1[sample(.N, 103487)]
vcf2_sub_num1 <- vcf2_sub1[,c(meta_ids), with=FALSE]
vcf2_sub2 <- vcf_nz2[sample(.N, 103487)]
vcf2_sub_num2 <- vcf2_sub2[,c(meta_ids), with=FALSE]
vcf2_sub3 <- vcf_nz3[sample(.N, 103487)]
vcf2_sub_num3 <- vcf2_sub3[,c(meta_ids), with=FALSE]
vcf2_sub4 <- vcf_nz4[sample(.N, 103487)]
vcf2_sub_num4 <- vcf2_sub4[,c(meta_ids), with=FALSE]
vcf2_sub5 <- vcf_nz5[sample(.N, 103487)]
vcf2_sub_num5 <- vcf2_sub5[,c(meta_ids), with=FALSE]

vcf2_snps <- rbind(vcf2_sub_num1, vcf2_sub_num2, vcf2_sub_num3, vcf2_sub_num4, vcf2_sub_num5)
vcf2_trans <- transpose(vcf2_snps)
vcf2_trans_full <- cbind(vcf2_trans, meta)
res.pca <- PCA(vcf2_trans,  graph = FALSE) # PCA on JUST numerical

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-15_PCA_snps_scree_v2_s2.pdf")
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()

basic_plot_snp <- fviz_pca_ind(res.pca, label="none")
levels.clade <- c("AstCal", "Benthic", "Deep", "Diplo", "Mbuna","Rhampho", "Utaka")
col_clade <- "clade"
vcf2_trans_full[, (col_clade) := factor(get(col_clade), levels = levels.clade)]
up_plot_snp <- ggplot(cbind(basic_plot_snp$data,vcf2_trans_full[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=4) + 
    theme_bw(base_size=18) +
    labs(title="PCA - polymorphic SNPs",x = "PC1 (2.7%)", y = "PC2 (2.3%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=popPalette) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-15_up_plot_snps_v2_s2.pdf", width=16, height=8)
up_plot_snp
dev.off()

basic_plot_snp_34 <- fviz_pca_ind(res.pca, label="none", axes=c(3,4))
up_plot_snp_34 <- ggplot(cbind(basic_plot_snp_34$data,vcf_trans_full[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=4) + 
    theme_bw(base_size=18) +
    labs(title="PCA - polymorphic SNPs",x = "PC3 (2%)", y = "PC4 (1.9%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=popPalette) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-15_up_plot_snp_34_v2.pdf", width=16, height=8)
up_plot_snp_34
dev.off()

#### comparison to minor allele frequency
sum_ins_all <- vcf_snps[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ids]
vcf_maf <- as.data.table(sum_ins_all)
vcf_maf <- vcf_maf[,total_alleles := 360] ## 2 * 180 samples
vcf_maf <- vcf_maf[,maf := sum_ins_all / total_alleles]
vcf_maf <- vcf_maf[order(maf)]

maf_all_snp <-ggplot(vcf_maf, aes(x=maf)) + 
  geom_histogram(binwidth=0.01)+
  geom_vline(xintercept=c(0.05), color="red")+
  labs(title="Minor allele frequency: SNPs",x = "Minor allele frequency", y = "Frequency")+
  theme_bw(base_size = 26)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-15_maf_snps_v2.pdf", width=12, height=12)
maf_all_snp
dev.off()

vcf_maf_05 <- vcf_maf[,.(maf_05=length(which(maf>0.05)) / length(maf))][order(-maf_05)]

