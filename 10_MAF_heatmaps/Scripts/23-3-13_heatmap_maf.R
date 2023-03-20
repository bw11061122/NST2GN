### 23-3-13
## heatmap on minor allele frequency 
library(ggplot2)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(reshape)
library(tidyr)
library(tidyverse)
library(gplots)

meta <- fread("~/Desktop/KEY DATA/23-2-16_cichlid_meta_standardized_dataset.csv",
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

### subset meta for specific species
meta_ap <- unlist(meta[pop_id == "Alticorpus peterdaviesi.Cape_Maclear", primary_id])
meta_ask <- unlist(meta[pop_id == "Astatotilapia calliptera.Lake_Kingiri", primary_id])
meta_asm <- unlist(meta[pop_id == "Astatotilapia calliptera.Lake_Masoko", primary_id])
meta_crhoa <- unlist(meta[pop_id == "Chilotilapia rhoadesii", primary_id])
meta_cch <- unlist(meta[pop_id == "Copadichromis chrysonotus.Lake_Malombe", primary_id])
meta_cvm <- unlist(meta[pop_id == "Copadichromis virginalis.Lake_Malombe", primary_id])
meta_cvs <- unlist(meta[pop_id == "Copadichromis virginalis.Southwest_arm", primary_id])
meta_czeb <- unlist(meta[pop_id == "Cynotilapia zebroides.Cape_Maclear", primary_id])
meta_diplo <- unlist(meta[pop_id == "Diplotaxodon limnothrissa.Southwest_arm", primary_id])
meta_fr <- unlist(meta[pop_id == "Fossorochromis rostratus.Lake_Malombe", primary_id])
meta_fuel <- unlist(meta[pop_id == "Labeotropheus fuelleborni.Chilumba", primary_id])
meta_trew <- unlist(meta[pop_id == "Labeotropheus trewavasae.Chilumba", primary_id])
meta_mzeb <- unlist(meta[pop_id == "Maylandia zebra.Cape_Maclear", primary_id])
meta_oto <- unlist(meta[pop_id == "Otopharynx argyrosoma.Southeast_arm", primary_id])
meta_rham <- unlist(meta[pop_id == "Rhamphochromis longiceps", primary_id])

### read the maf dataset
### this dataset contains minor allele frequency (frequency of each insertion in each population)
maf <- fread("~/Desktop/FIGURES/23-3-13_maf_by_pop.txt", 
             sep = "\t", header=T, fill=TRUE)
## factors  
maf <- maf[type=="MEI",] ## only interested in new insertions
maf <- maf[, order:=as.factor(order)]
maf <- maf[, family:=as.factor(family)]
maf <- maf[, superfamily:=as.factor(superfamily)]
upd.cols.maf = sapply(maf, is.factor)
maf <- maf[, names(maf)[upd.cols.maf] := lapply(.SD, factor), .SDcols = upd.cols.maf] ## this gets rid of levels which are absent in the dt 

## analysis at the family level
maf_sub_fam <- maf[,c(2,6:20)]
maf_agg_fam <- aggregate(.~family, as.data.frame(maf_sub_fam), mean) ## average maf by family 
fam_names <- maf_agg_fam$family
maf_ann_fam <- maf %>% select(c(order, family)) %>% filter(family %in% fam_names) 
maf_ann_fam  <- unique(maf_ann_fam) 
maf_merge_fam <- merge(maf_ann_fam, as.data.table(maf_agg_fam), by="family")
maf_merge_fam <- maf_merge_fam[order(order, family)]
maf_merge_fam$combo <- paste(maf_merge_fam$order, maf_merge_fam$family, sep="_")
c_rem <- c("DNA_Rex-Babar-6", "DNA_Rex-Babar-13", 
           "DNA_Rex-Babar-16", "DNA_Rex-Babar-12",
           "DNA_Dong-R4-1", "DNA_Dong-R4-2",
           "DNA_L2-13", "DNA_L2-17", "DNA_L2-6",
           "DNA_Ngaro-1", "DNA_Gypsy-35",
           "DNA_Pao-36", "DNA_L1-3", "DNA_L1-18",
           "DNA_Unknown-12", "DNA_Unknown-63",
           "DNA_Unknown-24", "DNA_L2-32", "DNA_Rex-Babar-17")
maf_merge_fam <- maf_merge_fam %>% filter(!(combo %in% c_rem))
maf_merge_fam <- maf_merge_fam[order(order, family)]
maf_num_fam <- maf_merge_fam[,3:17]
maf_num_fam <- maf_num_fam %>% as.matrix()
rownames(maf_num_fam) <- maf_merge_fam$family
colnames(maf_num_fam) <- c("Alticorpus peterdaviesi.Cape_Maclear",
                                  "Astatotilapia calliptera.Lake_Kingiri", "Astatotilapia calliptera.Lake_Masoko", "Chilotilapia rhoadesii",
                                  "Copadichromis chrysonotus.Lake_Malombe", "Copadichromis virginalis.Lake_Malombe", "Copadichromis virginalis.Southwest_arm", 
                                  "Cynotilapia zebroides.Cape_Maclear", "Diplotaxodon limnothrissa.Southwest_arm", "Fossorochromis rostratus.Lake_Malombe", 
                                  "Labeotropheus fuelleborni.Chilumba", "Labeotropheus trewavasae.Chilumba", "Maylandia zebra.Cape_Maclear", 
                                  "Otopharynx argyrosoma.Southeast_arm", "Rhamphochromis longiceps")  
## pheatmap package
col_ann <- meta %>% select(clade, pop_id)
col_ann <- unique(col_ann)
col_ann <- col_ann %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "pop_id")

## row annotation 
row_ann_maf_fam <- maf_merge_fam %>% select(order, family) %>% column_to_rownames(var="family")

## set colors 
annotation_colors = list(
  clade = c(AstCal="#AED581", Mbuna="#5E35B1", Utaka="#1B5E20", Benthic="#B71C1C", 
            Deep="#1A237E", Rhampho="#EF6C00", Diplo="#827717"),
  order = c(DNA="#EF5350", LINE="#2196F3", LTR="#43A047", SINE="#9575CD",
            RC="#D1C4E9"))

## modify the column names function in pheatmap to get rotated labels
## courtesy of https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps 
library(grid) 
draw_colnames_90 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                  ns=asNamespace("pheatmap"))
## heat map 
pheatmap(maf_num_fam, scale="row", annotation_row = row_ann_maf_fam,
         annotation_col = col_ann,
         cellwidth=20,
         annotation_colors = annotation_colors,
         main="Families: MAF", cluster_rows=F, 
         cluster_cols=T, show_rownames = F,
         fontsize=22, cexCol=8)

## analysis at the superfamily level
maf_sub_sfam <- maf[,c(4,6:20)]
maf_agg_sfam <- aggregate(.~superfamily, as.data.frame(maf_sub_sfam), mean) ## average by superfamily
sfam_names <- maf_agg_sfam$superfamily
maf_ann_sfam <- maf %>% select(c(order, superfamily)) %>% filter(superfamily %in% sfam_names) 
maf_ann_sfam  <- unique(maf_ann_sfam) 
maf_merge_sfam <- merge(maf_ann_sfam, as.data.table(maf_agg_sfam), by="superfamily")
maf_merge_sfam <- maf_merge_sfam[order(order, superfamily)]
maf_merge_sfam$combo <- paste(maf_merge_sfam$order, maf_merge_sfam$superfamily, sep="_")
c_rem_sfam <- c("DNA_Dong", "DNA_Gypsy", "DNA_L1", "DNA_L2",
           "DNA_Pao", "DNA_Ngaro", "DNA_Rex", "DNA_Unknown")
maf_merge_sfam <- maf_merge_sfam %>% filter(!(combo %in% c_rem_sfam))
maf_merge_sfam <- maf_merge_sfam[order(order, superfamily)]
maf_num_sfam <- maf_merge_sfam[,3:17]
maf_num_sfam <- maf_num_sfam %>% as.matrix()
rownames(maf_num_sfam) <- maf_merge_sfam$superfamily
colnames(maf_num_sfam) <- c("Alticorpus peterdaviesi.Cape_Maclear",
                            "Astatotilapia calliptera.Lake_Kingiri", "Astatotilapia calliptera.Lake_Masoko", "Chilotilapia rhoadesii",
                            "Copadichromis chrysonotus.Lake_Malombe", "Copadichromis virginalis.Lake_Malombe", "Copadichromis virginalis.Southwest_arm", 
                            "Cynotilapia zebroides.Cape_Maclear", "Diplotaxodon limnothrissa.Southwest_arm", "Fossorochromis rostratus.Lake_Malombe", 
                            "Labeotropheus fuelleborni.Chilumba", "Labeotropheus trewavasae.Chilumba", "Maylandia zebra.Cape_Maclear", 
                            "Otopharynx argyrosoma.Southeast_arm", "Rhamphochromis longiceps")  

## row annotation 
row_ann_maf_sfam <- maf_merge_sfam %>% select(order, superfamily) %>% column_to_rownames(var="superfamily")

## heat map 
pheatmap(maf_num_sfam, annotation_row = row_ann_maf_sfam,
         annotation_col = col_ann,
         cellwidth=16,
         cellheight=16,
         annotation_colors = annotation_colors,
         main="MAF in superfamilies, not scaled", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=18, cexCol=8) ### no scaling 
pheatmap(maf_num_sfam, scale="row", annotation_row = row_ann_maf_sfam,
         annotation_col = col_ann,
         cellwidth=16,
         cellheight=16,
         annotation_colors = annotation_colors,
         main="MAF in superfamilies, scaled by superfamily (row)", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=18, cexCol=8) ### scale by row 
pheatmap(maf_num_sfam, scale="column", annotation_row = row_ann_maf_sfam,
         annotation_col = col_ann,
         cellwidth=16,
         cellheight=16,
         annotation_colors = annotation_colors,
         main="MAF in superfamilies, scaled by population (col)", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=18, cexCol=8) ### scale by column

cor(c(maf_num_sfam), c(vcf_df_num_sfam_ins)) # 0.2323
cor.test(c(maf_num_sfam), c(vcf_df_num_sfam_ins)) # 2.298e-06

### get everything to one column
maf_values <- c(maf_num_sfam)
ins_values <- c(vcf_df_num_sfam_ins)
plot(ins_values, maf_values)

par(mfrow=c(1,2))
heatmap(maf_num_sfam)
heatmap(vcf_df_num_sfam_ins)

matrix_diff <- vcf_df_num_sfam_ins / maf_num_sfam
pheatmap(matrix_diff, scale="row")


