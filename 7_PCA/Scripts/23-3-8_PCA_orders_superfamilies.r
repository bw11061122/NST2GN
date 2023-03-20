###### 2023-03-08 update 
###### redo with different TE orders 
##### TE orders to look at: DNA, LINE, LTR, SINE
##### TE superfamilies: Gypsy, Pao 
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

### see how the full vcf was subsetted in script 23-2-16_subset_vcf_new_meta.r
## Read the subsetted vcf
vcf <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_ins_del_all_final.txt", 
    sep = "\t", header=T, fill=TRUE)

## split information about transposons 
vcf <- vcf[, order:=as.factor(order)]
vcf <- vcf[, family:=as.factor(family)]
vcf <- vcf[, superfamily:=as.factor(superfamily)]

## PCA across different orders 
info_col <- c("order", "family", "superfamily", "type") # columns to remove 
vcf_DNA <- vcf[order=="DNA",]
vcf_col_info_DNA <- vcf_DNA[,c(info_col), with=FALSE]
vcf_num_DNA <- vcf_DNA[,c(meta_ids), with=FALSE]
vcf_trans_DNA <- transpose(vcf_num_DNA)
vcf_trans_full_DNA <- cbind(vcf_trans_DNA, meta)
res.pca.DNA <- PCA(vcf_trans_DNA,  graph = FALSE) # PCA on JUST numerical
PCA_cichlid_DNA_clade <- fviz_pca_ind(res.pca.DNA,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_DNA$clade), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

PCA_cichlid_DNA_pop <- fviz_pca_ind(res.pca.DNA,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_DNA$pop_id), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_DNA_clade.pdf", width=14, height=14)
PCA_cichlid_DNA_clade
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_DNA_pop.pdf", width=14, height=14)
PCA_cichlid_DNA_pop
dev.off()

vcf_LINE <- vcf[order=="LINE",]
vcf_col_info_LINE <- vcf_LINE[,c(info_col), with=FALSE]
vcf_num_LINE <- vcf_LINE[,c(meta_ids), with=FALSE]
vcf_trans_LINE <- transpose(vcf_num_LINE)
vcf_trans_full_LINE <- cbind(vcf_trans_LINE, meta)
res.pca.LINE <- PCA(vcf_trans_LINE,  graph = FALSE) # PCA on JUST numerical
PCA_cichlid_LINE_clade <- fviz_pca_ind(res.pca.LINE,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_LINE$clade), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

PCA_cichlid_LINE_pop <- fviz_pca_ind(res.pca.LINE,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_LINE$pop_id), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_LINE_clade.pdf", width=14, height=14)
PCA_cichlid_LINE_clade
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_LINE_pop.pdf", width=14, height=14)
PCA_cichlid_LINE_pop
dev.off()

vcf_LTR <- vcf[order=="LTR",]
vcf_col_info_LTR <- vcf_LTR[,c(info_col), with=FALSE]
vcf_num_LTR <- vcf_LTR[,c(meta_ids), with=FALSE]
vcf_trans_LTR <- transpose(vcf_num_LTR)
vcf_trans_full_LTR <- cbind(vcf_trans_LTR, meta)
res.pca.LTR <- PCA(vcf_trans_LTR,  graph = FALSE) # PCA on JUST numerical
PCA_cichlid_LTR_clade <- fviz_pca_ind(res.pca.LTR,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_LTR$clade), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

PCA_cichlid_LTR_pop <- fviz_pca_ind(res.pca.LTR,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_LTR$pop_id), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_LTR_clade.pdf", width=14, height=14)
PCA_cichlid_LTR_clade
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_LTR_pop.pdf", width=14, height=14)
PCA_cichlid_LTR_pop
dev.off()

vcf_SINE <- vcf[order=="SINE",]
vcf_col_info_SINE <- vcf_SINE[,c(info_col), with=FALSE]
vcf_num_SINE <- vcf_SINE[,c(meta_ids), with=FALSE]
vcf_trans_SINE <- transpose(vcf_num_SINE)
vcf_trans_full_SINE <- cbind(vcf_trans_SINE, meta)
res.pca.SINE <- PCA(vcf_trans_SINE,  graph = FALSE) # PCA on JUST numerical
PCA_cichlid_SINE_clade <- fviz_pca_ind(res.pca.SINE,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_SINE$clade), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

PCA_cichlid_SINE_pop <- fviz_pca_ind(res.pca.SINE,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_SINE$pop_id), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_SINE_clade.pdf", width=14, height=14)
PCA_cichlid_SINE_clade
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_SINE_pop.pdf", width=14, height=14)
PCA_cichlid_SINE_pop
dev.off()

vcf_Gypsy <- vcf[superfamily=="Gypsy",]
vcf_col_info_Gypsy <- vcf_Gypsy[,c(no_col), with=FALSE]
vcf_num_Gypsy <- vcf_Gypsy[,c(meta_ids), with=FALSE]
vcf_trans_Gypsy <- transpose(vcf_num_Gypsy)
vcf_trans_full_Gypsy <- cbind(vcf_trans_Gypsy, meta)
res.pca.Gypsy <- PCA(vcf_trans_Gypsy,  graph = FALSE) # PCA on JUST numerical
PCA_cichlid_Gypsy_clade <- fviz_pca_ind(res.pca.Gypsy,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_Gypsy$clade), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

PCA_cichlid_Gypsy_pop <- fviz_pca_ind(res.pca.Gypsy,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_Gypsy$pop_id), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_Gypsy_clade.pdf", width=14, height=14)
PCA_cichlid_Gypsy_clade
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_Gypsy_pop.pdf", width=14, height=14)
PCA_cichlid_Gypsy_pop
dev.off()

vcf_Pao <- vcf[superfamily=="Pao",]
vcf_col_info_Pao <- vcf_Pao[,c(no_col), with=FALSE]
vcf_num_Pao <- vcf_Pao[,c(meta_ids), with=FALSE]
vcf_trans_Pao <- transpose(vcf_num_Pao)
vcf_trans_full_Pao <- cbind(vcf_trans_Pao, meta)
res.pca.Pao <- PCA(vcf_trans_Pao,  graph = FALSE) # PCA on JUST numerical
PCA_cichlid_Pao_pop <- fviz_pca_ind(res.pca.Pao,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_Pao$pop_id), pointsize = 3, # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

PCA_cichlid_Pao_clade <- fviz_pca_ind(res.pca.Pao,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_Pao$clade),
                                     pointsize = 3, # color by location, set size
                                     title="PCA - TE insertions")
                                     # addEllipses = TRUE # Concentration ellipses


pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_Pao_clade.pdf", width=14, height=14)
PCA_cichlid_Pao_clade + theme(legend.text = element_text(size = 20),
                                     axis.text = element_text(size = 15))

dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_Pao_pop.pdf", width=14, height=14)
PCA_cichlid_Pao_pop
dev.off()

vcf_PiggyBac <- vcf[superfamily=="PiggyBac",]
vcf_col_info_PiggyBac <- vcf_PiggyBac[,c(no_col), with=FALSE]
vcf_num_PiggyBac <- vcf_PiggyBac[,c(meta_ids), with=FALSE]
vcf_trans_PiggyBac <- transpose(vcf_num_PiggyBac)
vcf_trans_full_PiggyBac <- cbind(vcf_trans_PiggyBac, meta)
res.pca.PiggyBac <- PCA(vcf_trans_PiggyBac,  graph = FALSE) # PCA on JUST numerical
PCA_cichlid_PiggyBac_clade <- fviz_pca_ind(res.pca.PiggyBac,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_PiggyBac$clade), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

PCA_cichlid_PiggyBac_pop <- fviz_pca_ind(res.pca.PiggyBac,
                                     label = "none", # hide individual labels
                                     habillage = as.factor(vcf_trans_full_PiggyBac$pop_id), # color by location
                                     title="PCA - TE insertions"
                                     # addEllipses = TRUE # Concentration ellipses
)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_PiggyBac_clade.pdf", width=14, height=14)
PCA_cichlid_PiggyBac_clade
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PCA_PiggyBac_pop.pdf", width=14, height=14)
PCA_cichlid_PiggyBac_pop
dev.off()

### added 23-3-9 - graphics
# set colors and shapes
library(RColorBrewer)
coul_pop <- brewer.pal(8, "Set1") 
coul_pop <- colorRampPalette(coul_pop)(15)
shapes <- c(3, 4, 8, 15, 16, 17, 18) ## chose to not have "empty" shapes
#### DNA
## scree plot
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_PCA_final_scree_plot_DNA.pdf")
fviz_eig(res.pca.DNA, addlabels = TRUE, ylim = c(0, 50))
dev.off()
### basic plot
basic_plot_DNA <- fviz_pca_ind(res.pca.DNA, label="none")
### upgraded plot
up_plot_DNA <- ggplot(cbind(basic_plot_DNA$data,vcf_trans_full_DNA[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=3) + 
    theme_bw(base_size=18) +
    labs(title="PCA on polymorphic insertions: DNA TEs",x = "PC1 (3.3%)", y = "PC2 (3.2%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=coul_pop) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_up_plot_DNA.pdf", width=16, height=8)
up_plot_DNA
dev.off()

#### LINE
## scree plot
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_PCA_final_scree_plot_LINE.pdf")
fviz_eig(res.pca.LINE, addlabels = TRUE, ylim = c(0, 50))
dev.off()
### basic plot
basic_plot_LINE <- fviz_pca_ind(res.pca.LINE, label="none")
### upgraded plot
up_plot_LINE <- ggplot(cbind(basic_plot_LINE$data,vcf_trans_full_LINE[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=3) + 
    theme_bw(base_size=18) +
    labs(title="PCA on polymorphic insertions: LINE TEs",x = "PC1 (4.2%)", y = "PC2 (4%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=coul_pop) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_up_plot_LINE.pdf", width=16, height=8)
up_plot_LINE
dev.off()

#### LTR
## scree plot
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_PCA_final_scree_plot_LTR.pdf")
fviz_eig(res.pca.LTR, addlabels = TRUE, ylim = c(0, 50))
dev.off()
### basic plot
basic_plot_LTR <- fviz_pca_ind(res.pca.LTR, label="none")
### upgraded plot
up_plot_LTR <- ggplot(cbind(basic_plot_LTR$data,vcf_trans_full_LTR[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=3) + 
    theme_bw(base_size=18) +
    labs(title="PCA on polymorphic insertions: LTR TEs",x = "PC1 (2.7%)", y = "PC2 (2.6%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=coul_pop) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_up_plot_LTR.pdf", width=16, height=8)
up_plot_LTR
dev.off()

#### SINE
## scree plot
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_PCA_final_scree_plot_SINE.pdf")
fviz_eig(res.pca.SINE, addlabels = TRUE, ylim = c(0, 50))
dev.off()
### basic plot
basic_plot_SINE <- fviz_pca_ind(res.pca.SINE, label="none")
### upgraded plot
up_plot_SINE <- ggplot(cbind(basic_plot_SINE$data,vcf_trans_full_SINE[,c("clade","pop_id")]),
    aes(x=x,y=y,col=factor(pop_id), shape=factor(clade))) + 
    geom_point(size=3) + 
    theme_bw(base_size=18) +
    labs(title="PCA on polymorphic insertions: SINE TEs",x = "PC1 (5.8%)", y = "PC2 (5.4%)")+
    guides(shape = guide_legend(title = "clade"), col = guide_legend(title = "population"))+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=coul_pop) 
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_up_plot_SINE.pdf", width=16, height=8)
up_plot_SINE
dev.off()
### sbatch -c 8 --time=01:00:00 --wrap="Rscript 23-3-8_PCA_orders_superfamilies.r"