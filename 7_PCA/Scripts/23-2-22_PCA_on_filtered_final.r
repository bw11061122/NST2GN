#!/usr/bin/env Rscript

### 2023-2-22
## In this script, I am going to try to deal with the new set of cichlids that we want to look at 
## first, subset the dataset for MEI cichlids
## do encoding and filtering on the correct dataset 
## then do PCA and analyse data structure
## modified 23-3-9 to sort out plot graphics 

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

## set colors for plotting
cbPalette <- c("#EF5350", "#2196F3", "#43A047", "#880E4F", "#FFC107", "#D1C4E9", "#757575", "#795548") # orders
popPalette <- c("#9CCC65", "#7CB342", "#7E57C2", "512DA8", "#AB47BC", "#6A1B9A", "#1B5E20", "8BC34A", "388E3C",
"#B71C1C", "#E57373", "#D32F2F", "#1A237E", "#BF360C", "#827717")
cladePalette <- c("#AED581", "#5E35B1", "#1B5E20", "#B71C1C", "#1A237E", "#EF6C00", "#827717")

#################################################
#################################################
## to see how I subsetted cichlids such that they match what Richard suggested on 23-2-14 
## see script 23-2-14_standardized dataset (locally and also in this folder)

#################################################
#################################################
## read metadata
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
setnames(meta, "name_loc", "population")
meta <- meta[, pop_id := population][(genus %in% c("Rhamphochromis", "Chilotilapia")) | is.na(population), pop_id := names][]
meta <- meta[, pop_id:=as.factor(pop_id)]
upd.cols = sapply(meta, is.factor)
meta <- meta[, names(meta)[upd.cols] := lapply(.SD, factor), .SDcols = upd.cols] ## this gets rid of levels which are absent in the dt 
print("loaded metadata")
levels.pop <- c("Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
"Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
"Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
"Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
"Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", "Diplotaxodon limnothrissa.Southwest_arm")
col <- "pop_id"
meta[, (col) := factor(get(col), levels = levels.pop)]

### see how the full vcf was subsetted in script 23-2-16_subset_vcf_new_meta.r
## Read the subsetted vcf
vcf <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_ins_del_all_final.txt", 
    sep = "\t", header=T, fill=TRUE)

## split information about transposons 
vcf <- vcf[, order:=as.factor(order)]
vcf <- vcf[, fam:=as.factor(fam)]
vcf <- vcf[, superfamily:=as.factor(superfamily)]
levels.order <- c("DNA", "LINE", "LTR", "SINE", "Unknown", "RC", "Retrotransposon", "rRNA")
col <- "order"
vcf[, (col) := factor(get(col), levels = levels.order)]

## PCA for all transposons?
no_col <- c("#CHROM", "POS", "ID", "REF", "INFO", 
"SVTYPE", "TE", "START", "END", "LENGTH", "AC", "order",
"family", "superfamily", "type") # columns to remove 
vcf_col_info <- vcf[,c(no_col), with=FALSE]
vcf_num <- vcf[,c(meta_ids), with=FALSE]
vcf_full <- cbind(vcf, vcf_col_info)

res.pca.MEI <- PCA(vcf_num,  graph = FALSE) # PCA on JUST numerical
PCA_MEI_order <- fviz_pca_ind(res.pca.MEI,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_full$order),# color by location
                                     title="PCA on all polymorphic insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

PCA_MEI_fam <- fviz_pca_ind(res.pca.MEI,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_full$fam), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

vcf_trans <- transpose(vcf_num)
vcf_trans_full <- cbind(vcf_trans, meta)

res.pca.MEI_trans <- PCA(vcf_trans,  graph = FALSE) # PCA on JUST numerical

PCA_MEI_order_cichlid <- fviz_pca_ind(res.pca.MEI_trans,
                                     label = "none", # hide individual labels
                                     title="PCA on all polymorphic insertions",
                                     habillage=as.factor(vcf_trans_full$pop_id))+
                                    # scale_fill_manual(values=coul_pop)) # color by location
                                     # addEllipses = TRUE # Concentration ellipses


PCA_MEI_fam_cichlid <- fviz_pca_ind(res.pca.MEI_trans,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full$pop_id), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-24_PCA_final_scree_plot.pdf")
fviz_eig(res.pca.MEI_trans, addlabels = TRUE, ylim = c(0, 50))
dev.off()

PCA_MEI_fam_cichlid34 <- fviz_pca_ind(res.pca.MEI_trans, axes = c(3,4),
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full$pop_id), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)
PCA_MEI_order_cichlid34 <- fviz_pca_ind(res.pca.MEI_trans, axes = c(3,4),
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full$clade), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-22_PCA_MEI_by_TE_order.pdf")
PCA_MEI_order
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-22_PCA_MEI_by_TE_fam.pdf")
PCA_MEI_fam
dev.off()

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-22_PCA_final_clade.pdf")
PCA_MEI_order_cichlid
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-22_PCA_final_name_loc.pdf")
PCA_MEI_fam_cichlid
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-24_PCA_final_group.pdf")
PCA_MEI_fam_cichlid
dev.off()

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-24_PCA_final_clade_34.pdf")
PCA_MEI_order_cichlid34
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-24_PCA_final_group_34.pdf")
PCA_MEI_fam_cichlid34
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-24_PCA_final_group.pdf")
PCA_MEI_fam_cichlid
dev.off()

#####################################################################
#####################################################################
### added 23-3-9
### redo scree plot

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_PCA_final_scree_plot.pdf")
fviz_eig(res.pca.MEI_trans, addlabels = TRUE, ylim = c(0, 50))
dev.off()

# set shapes
shapes <- c(3, 4, 8, 15, 16, 17, 18) ## chose to not have "empty" shapes

### basic plot
basic_plot_cichlid <- fviz_pca_ind(res.pca.MEI_trans, label="none")
### upgraded plot
up_plot_cichlid <- ggplot(cbind(basic_plot_cichlid$data,vcf_trans_full[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=4) + 
    theme_bw(base_size=18) +
    labs(title="PCA on all polymorphic insertions",x = "PC1 (3.4%)", y = "PC2 (3.4%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=popPalette) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-10_up_plot_cichlid.pdf", width=16, height=8)
up_plot_cichlid
dev.off()

basic_plot_cichlid_34 <- fviz_pca_ind(res.pca.MEI_trans, label="none", axes=c(3,4))
### upgraded plot
up_plot_cichlid_34 <- ggplot(cbind(basic_plot_cichlid_34$data,vcf_trans_full[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=4) + 
    theme_bw(base_size=18) +
    labs(title="PCA on all polymorphic insertions",x = "PC3 (2.9%)", y = "PC4 (2.1%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=popPalette) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-10_up_plot_cichlid_34.pdf", width=16, height=8)
up_plot_cichlid_34
dev.off()

## Rscript 23-2-22_PCA_on_filtered_final.r
