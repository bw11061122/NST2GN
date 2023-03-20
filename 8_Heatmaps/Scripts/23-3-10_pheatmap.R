## 23-3-10 
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

### read the vcf dataset
vcf <- fread("~/Desktop/KEY DATA/23-2-27_vcf_ins_del_all_final.txt", 
             sep = "\t", header=T, fill=TRUE)

## factors  
vcf <- vcf[, order:=as.factor(order)]
vcf <- vcf[, family:=as.factor(family)]
vcf <- vcf[, superfamily:=as.factor(superfamily)]
upd.cols.vcf = sapply(vcf, is.factor)
vcf <- vcf[, names(vcf)[upd.cols.vcf] := lapply(.SD, factor), .SDcols = upd.cols.vcf] ## this gets rid of levels which are absent in the dt 

### how i want the heatmap to look like:
### columns: individuals 
### rows: superfamilies 
### I will split by insertions and deletions because I think 
### this way it will be more informative
sum_ins_ap <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ap]
sum_ins_ask <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ask]
sum_ins_asm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_asm]
sum_ins_cch <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cch]
sum_ins_crhoa <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_crhoa]
sum_ins_cvm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cvs]
sum_ins_cvs <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cvm]
sum_ins_czeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_czeb]
sum_ins_diplo <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_diplo]
sum_ins_fr <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_fr]
sum_ins_fuel <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_fuel]
sum_ins_trew <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_trew]
sum_ins_mzeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_mzeb]
sum_ins_oto <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_oto]
sum_ins_rham <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_rham]

vcf_sub_fam_ins <- vcf[type=="MEI","family",with=FALSE]
vcf_dt_sums_fam_ins <- as.data.table(cbind(vcf_sub_fam_ins, sum_ins_ap, sum_ins_ask, sum_ins_asm, 
                                       sum_ins_crhoa, sum_ins_cch, sum_ins_cvm, sum_ins_cvs, sum_ins_czeb, 
                                       sum_ins_diplo, sum_ins_fr, sum_ins_fuel,  sum_ins_trew, sum_ins_mzeb, 
                                       sum_ins_oto, sum_ins_rham))

vcf_sub_sfam_ins <- vcf[type=="MEI","superfamily",with=FALSE]
vcf_dt_sums_sfam_ins <- as.data.table(cbind(vcf_sub_sfam_ins, sum_ins_ap, sum_ins_ask, sum_ins_asm, 
                                       sum_ins_crhoa, sum_ins_cch, sum_ins_cvm, sum_ins_cvs, sum_ins_czeb, 
                                       sum_ins_diplo, sum_ins_fr, sum_ins_fuel,  sum_ins_trew, sum_ins_mzeb, 
                                       sum_ins_oto, sum_ins_rham))
sums <- c("sum_ins_ap", "sum_ins_ask", "sum_ins_asm", 
"sum_ins_crhoa", "sum_ins_cch", "sum_ins_cvm", "sum_ins_cvs", "sum_ins_czeb", 
"sum_ins_diplo", "sum_ins_fr", "sum_ins_fuel",  "sum_ins_trew", "sum_ins_mzeb", 
"sum_ins_oto", "sum_ins_rham")

vcf_df_agg_fam_ins <- aggregate(.~family, as.data.frame(vcf_dt_sums_fam_ins), sum) ## sum by family 
fam_names <- vcf_df_agg_fam_ins$family
df_ins_ann_fam <- vcf %>% filter(type=="MEI") %>% select(c(order, family)) %>% filter(family %in% fam_names) 
df_ins_ann_fam  <- unique(df_ins_ann_fam) 
vcf_df_ins_fam_merge <- merge(df_ins_ann_fam, as.data.table(vcf_df_agg_fam_ins), by="family")
vcf_df_ins_fam_merge <-vcf_df_ins_merge[order(order, family)]
vcf_df_num_fam_ins <- vcf_df_ins_merge %>% select(sums) 
vcf_df_num_fam_ins <- vcf_df_num_fam_ins%>% as.matrix()
rownames(vcf_df_num_fam_ins) <- vcf_df_ins_merge$family
colnames(vcf_df_num_fam_ins) <- c("Alticorpus peterdaviesi.Cape_Maclear",
                                  "Astatotilapia calliptera.Lake_Kingiri", "Astatotilapia calliptera.Lake_Masoko", "Chilotilapia rhoadesii",
                                  "Copadichromis chrysonotus.Lake_Malombe", "Copadichromis virginalis.Lake_Malombe", "Copadichromis virginalis.Southwest_arm", 
                                  "Cynotilapia zebroides.Cape_Maclear", "Diplotaxodon limnothrissa.Southwest_arm", "Fossorochromis rostratus.Lake_Malombe", 
                                  "Labeotropheus fuelleborni.Chilumba", "Labeotropheus trewavasae.Chilumba", "Maylandia zebra.Cape_Maclear", 
                                  "Otopharynx argyrosoma.Southeast_arm", "Rhamphochromis longiceps")  

vcf_df_agg_sfam_ins <- aggregate(.~superfamily, as.data.frame(vcf_dt_sums_sfam_ins), sum) ## sum by family 
sfam_names <- vcf_df_agg_sfam_ins$superfamily
df_ins_ann_sfam <- vcf %>% filter(type=="MEI") %>% select(c(order, superfamily)) %>% filter(superfamily %in% sfam_names) 
df_ins_ann_sfam  <- unique(df_ins_ann_sfam) 
vcf_df_ins_sfam_merge <- merge(df_ins_ann_sfam, as.data.table(vcf_df_agg_sfam_ins), by="superfamily")
vcf_df_ins_sfam_merge <-vcf_df_ins_sfam_merge[order(order, superfamily)]
vcf_df_num_sfam_ins <- vcf_df_ins_sfam_merge %>% select(sums) 
vcf_df_num_sfam_ins <- vcf_df_num_sfam_ins%>% as.matrix()
rownames(vcf_df_num_sfam_ins) <- vcf_df_ins_sfam_merge$superfamily
colnames(vcf_df_num_sfam_ins) <- c("Alticorpus peterdaviesi.Cape_Maclear",
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
row_ann_sfam <- vcf_df_ins_sfam_merge %>% select(order, superfamily) %>% column_to_rownames(var="superfamily")
row_ann_fam <- vcf_df_ins_fam_merge %>% select(order, family) %>% column_to_rownames(var="family")

## set colors 
annotation_colors = list(
  clade = c(AstCal="#AED581", Mbuna="#5E35B1", Utaka="#1B5E20", Benthic="#B71C1C", 
            Deep="#1A237E", Rhampho="#EF6C00", Diplo="#827717"),
  order = c(DNA="#EF5350", LINE="#2196F3", LTR="#43A047", SINE="#9575CD",
            RC="#D1C4E9"))

## modify the stupid column names function in pheatmap to get rotated labels
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
pheatmap(vcf_df_num_sfam_ins, scale="row", annotation_row = row_ann_sfam,
         annotation_col = col_ann,
         cellwidth=40,
         annotation_colors = annotation_colors,
         main="Superfamilies (insertions absent from reference)", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=22, cexCol=8)
pheatmap(vcf_df_num_fam_ins, scale="row", annotation_row = row_ann_fam,
         annotation_col = col_ann,
         cellwidth=f40,
         annotation_colors = annotation_colors,
         main="Families (insertions absent from reference)", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=22, cexCol=8)
pheatmap(vcf_df_num_sfam_del, scale="row", annotation_row = row_ann_sfam,
         annotation_col = col_ann,
         cellwidth=40,
         annotation_colors = annotation_colors,
         main="Superfamilies (insertions shared with reference)", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=22, cexCol=8)

## scale by columns
#### superfamilies 
pheatmap(vcf_df_num_sfam_del, scale="column", annotation_row = row_ann_sfam,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="Superfamilies (insertions shared with the reference)", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)
pheatmap(vcf_df_num_sfam_ins, scale="column", annotation_row = row_ann_sfam,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="Superfamilies (insertions absent from the reference)", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)

#### Added 23-3-13 to showcase specific families
#### looking at Maverick families 
maverick_fam <- vcf_df_num_fam_ins %>% as.data.frame() %>% rownames_to_column("family") %>% filter(grepl('Maverick', family))
maverick_fam_mt <- maverick_fam %>% column_to_rownames(var="family") %>% as.matrix()
row_ann_maverick <- vcf_df_ins_fam_merge %>% filter(grepl('Maverick', family)) %>% select(order, family) %>% column_to_rownames(var="family")
pheatmap(maverick_fam_mt, scale="row", annotation_row = row_ann_maverick,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="Maverick insertions", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)
pheatmap(maverick_fam_mt, #scale="column", 
         annotation_row = row_ann_maverick,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="Maverick insertions", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)
### DIRS
dirs_fam <- vcf_df_num_fam_ins %>% as.data.frame() %>% rownames_to_column("family") %>% filter(grepl('DIRS', family))
dirs_fam_mt <- dirs_fam %>% column_to_rownames(var="family") %>% as.matrix()
row_ann_dirs <- vcf_df_ins_fam_merge %>% filter(grepl('DIRS', family)) %>% select(order, family) %>% column_to_rownames(var="family")
pheatmap(dirs_fam_mt, scale="row", annotation_row = row_ann_dirs,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="DIRS insertions", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)
pheatmap(dirs_fam_mt, #scale="column", 
         annotation_row = row_ann_dirs,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="DIRS insertions", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)

dirs_fam <- vcf_df_num_fam_ins %>% as.data.frame() %>% rownames_to_column("family") %>% filter(grepl('dirs', family))
dirs_fam_mt <- dirs_fam %>% column_to_rownames(var="family") %>% as.matrix()
row_ann_dirs <- vcf_df_ins_fam_merge %>% filter(grepl('dirs', family)) %>% select(order, family) %>% column_to_rownames(var="family")
pheatmap(dirs_fam_mt, scale="row", annotation_row = row_ann_dirs,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="dirs insertions", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)
pheatmap(dirs_fam_mt, #scale="column", 
         annotation_row = row_ann_dirs,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="dirs insertions", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)

piggyback_fam <- vcf_df_num_fam_ins %>% as.data.frame() %>% rownames_to_column("family") %>% filter(grepl('PiggyBac', family))
piggyback_fam_mt <- piggyback_fam %>% column_to_rownames(var="family") %>% as.matrix()
row_ann_piggyback <- vcf_df_ins_fam_merge %>% filter(grepl('PiggyBac', family)) %>% select(order, family) %>% column_to_rownames(var="family")
pheatmap(piggyback_fam_mt, scale="row", annotation_row = row_ann_piggyback,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="PiggyBac insertions", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)
pheatmap(piggyback_fam_mt, #scale="column", 
         annotation_row = row_ann_piggyback,
         annotation_col = col_ann,
         cellwidth=10,
         cellheight=10,
         annotation_colors = annotation_colors,
         main="PiggyBac insertions", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=10, cexCol=8)

