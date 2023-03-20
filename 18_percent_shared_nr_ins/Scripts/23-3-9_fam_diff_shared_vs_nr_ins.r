#### [2023-03-09]
##### For each family, look at what % of insertions contributes to the total number of insertions?

.libPaths( c( "~/R/x86_64-redhat-linux-gnu-library/4.2" , .libPaths() ) ) 
.libPaths()
library(data.table)
library(gridExtra)
library(ggplot2)
library(plyr)
library(dplyr)
library(tibble)
library(ggrepel)
library("FactoMineR") ## PCA
library("factoextra") ## PCA 
library(RColorBrewer)
library(ggpmisc)
# library(umap) ## UMAP
# library(M3C)
print("loaded libraries successfully")

getwd()
# "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara"

### set for plottting
library(RColorBrewer)
coul_pop <- brewer.pal(8, "Set1") 
coul_pop <- colorRampPalette(coul_pop)(15)
cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999")

### read the metadata
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
vcf <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_ins_del_all_final.txt", 
    sep = "\t", header=T, fill=TRUE)

## factors  
vcf <- vcf[, order:=as.factor(order)]
vcf <- vcf[, family:=as.factor(family)]
vcf <- vcf[, superfamily:=as.factor(superfamily)]
upd.cols.vcf = sapply(vcf, is.factor)
vcf <- vcf[, names(vcf)[upd.cols.vcf] := lapply(.SD, factor), .SDcols = upd.cols.vcf] ## this gets rid of levels which are absent in the dt 
vcf[is.na(order==TRUE), order := "Unknown"]
levels.order <- c("DNA", "LINE", "LTR", "SINE", "Unknown", "RC", "Retrotransposon", "rRNA")
col <- "order"
vcf[, (col) := factor(get(col), levels = levels.order)]

######################################################################################
######################################################################################
### IDEA:
### you can divide the families into 4 classes
### percentage of insertions which are shared for each family 
### number of insertions identified in this family between ALL populations
### plot number of insertions ~ percentage shared insertions for each family
### Expected four classes: 
### 1 few and not shared - active very recently
### 2 many and not shared - highly active recently 
### 3 few and shared - not active recently, degraded / silenced?
### 4 many and shared - has been active in the ancestral population

## percentage of shared insertions for each family
## 1 how many samples have any insertion at all? (may also redo with !(x %in% c(0, 0.5))
sum_nz_ap <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_ap]
sum_nz_ask <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_ask]
sum_nz_asm <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_asm]
sum_nz_cch <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_cch]
sum_nz_crhoa <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_crhoa]
sum_nz_cvm <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_cvs]
sum_nz_cvs <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_cvm]
sum_nz_czeb <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_czeb]
sum_nz_diplo <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_diplo]
sum_nz_fr <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_fr]
sum_nz_fuel <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_fuel]
sum_nz_trew <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_trew]
sum_nz_mzeb <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_mzeb]
sum_nz_oto <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_oto]
sum_nz_rham <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_rham]

vcf_sub <- vcf[,c("ID", "family", "order", "superfamily", "type")]
vcf_sums_nz <- cbind(vcf_sub, sum_nz_ap, sum_nz_ask, sum_nz_asm, 
sum_nz_crhoa, sum_nz_cch, sum_nz_cvm, sum_nz_cvs, sum_nz_czeb, 
sum_nz_diplo, sum_nz_fr, sum_nz_fuel,  sum_nz_trew, sum_nz_mzeb, 
sum_nz_oto, sum_nz_rham) 
sums_nz <- c("sum_nz_ap", "sum_nz_ask", "sum_nz_asm", 
"sum_nz_crhoa", "sum_nz_cch", "sum_nz_cvm", "sum_nz_cvs", "sum_nz_czeb", 
"sum_nz_diplo", "sum_nz_fr", "sum_nz_fuel", "sum_nz_trew", "sum_nz_mzeb", 
"sum_nz_oto", "sum_nz_rham") 
sum_nz_all <- vcf_sums_nz[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=sums_nz]
vcf_sums_nz <- cbind(vcf_sums_nz, sum_nz_all)
vcf_sums_nz_melt <- melt(vcf_sums_nz)

### value !0 means that the population has at least one insertion identified
vcf_sum_nz_all <- vcf_sums_nz[sum_nz_all==15,c("family", "sum_nz_all"), with=FALSE]
sum_shared <- vcf_sum_nz_all[,.N, by=family]
setnames(sum_shared, "N", "N_shared_ins")
sum_total <- vcf[,.N, by=family]
setnames(sum_total, "N", "N_total_ins")
shared_merge <- merge(sum_shared, sum_total, by="family")
shared_merge <- shared_merge[, percent_shared := (N_shared_ins / N_total_ins) * 100][order(-percent_shared)]

### how do I get the percetage of insertions shared between all populations
nr_ins_fam <- vcf[,.N, by=family]
fam_data_merge <- merge(nr_ins_fam, shared_merge, by="family")

fam_data <- vcf[,c("family", "superfamily", "order"), with=FALSE]
fam_data_merge <- merge(fam_data_merge, fam_data, by="family") # color by order, superfamily, family
fam_data_merge[is.na(order==TRUE), order := "Unknown"]
levels.order <- c("DNA", "LINE", "LTR", "SINE", "Unknown", "RC", "Retrotransposon", "rRNA")
col <- "order"
fam_data_merge[, (col) := factor(get(col), levels = levels.order)]

corr_nr_ins_p_shared_order <- ggplot(data = fam_data_merge, aes(x = log(N), y = percent_shared, color=order))+
  geom_point(size=3)+
  labs(title="Correlation between number of insertions in a family and % insertions shared between all populations",x = "log(number of insertions in the family)", y = "% insertions shared in the family")+
  theme_bw(base_size = 24)+
  scale_color_manual(values=cbPalette)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_corr_nr_ins_p_shared_order.pdf", width=20, height=12)
corr_nr_ins_p_shared_order
dev.off()

### Add the total number of insertions and color by it
### For each family, what is the total number of 1 + 2 * 2
sum_ins_in_fam <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ids, by=family]
sum_ins_in_fam <- sum_ins_in_fam[,.(sum_ins=sum(V1)), by=.(family)]
fam_merge2 <- merge(fam_data_merge, sum_ins_in_fam, by="family") # color by order, superfamily, family
fam_merge2[is.na(order==TRUE), order := "Unknown"]
fam_merge2[, (col) := factor(get(col), levels = levels.order)]

corr_nr_ins_shared_sum_ins <- ggplot(data = fam_merge2, aes(x = log(N), y = percent_shared, color=log(sum_ins)))+
  geom_point(size=3)+
  scale_color_gradient(low = "black", high = "lightblue")+
  labs(title="Correlation between number of insertions and % insertions shared - DNA TEs",x = "log(number of insertions in the family)", y = "% insertions shared in the family")+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_corr_nr_ins_p_shared_sum_ins.pdf", width=20, height=12)
corr_nr_ins_shared_sum_ins
dev.off()

### split for individual orders
fam_data_merge_dna <- fam_data_merge[order=="DNA",] # color by order, superfamily, family
corr_nr_ins_shared_dna <- ggplot(data = fam_data_merge_dna, aes(x = log(N), y = percent_shared, color=superfamily))+
  geom_point(size=3)+
  labs(title="Correlation between number of insertions in a family and % insertions shared in all populations - DNA TEs",x = "log(number of insertions in the family)", y = "% insertions shared in the family")+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_corr_nr_ins_p_shared_dna.pdf", width=20, height=12)
corr_nr_ins_shared_dna
dev.off()

fam_data_merge_line <- fam_data_merge[order=="LINE",] # color by order, superfamily, family
corr_nr_ins_shared_line <- ggplot(data = fam_data_merge_line, aes(x = log(N), y = percent_shared, color=superfamily))+
  geom_point(size=3)+
  labs(title="Correlation between number of insertions and % insertions shared - - LINE TEs",x = "log(number of insertions in the family)", y = "% insertions shared in the family")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_corr_nr_ins_p_shared_line.pdf", width=20, height=12)
corr_nr_ins_shared_line
dev.off()

fam_data_merge_sine <- fam_data_merge[order=="SINE",] # color by order, superfamily, family
corr_nr_ins_shared_sine <- ggplot(data = fam_data_merge_sine, aes(x = log(N), y = percent_shared, color=superfamily))+
  geom_point(size=3)+
  labs(title="Correlation between number of insertions and % insertions shared - SINE TEs",x = "log(number of insertions in the family)", y = "% insertions shared in the family")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_corr_nr_ins_p_shared_sine.pdf", width=20, height=12)
corr_nr_ins_shared_sine
dev.off()

fam_data_merge_ltr <- fam_data_merge[order=="LTR",] # color by order, superfamily, family
corr_nr_ins_shared_ltr <- ggplot(data = fam_data_merge_ltr, aes(x = log(N), y = percent_shared, color=superfamily))+
  geom_point(size=3)+
  labs(title="Correlation between number of insertions and % insertions shared - LTR TEs",x = "log(number of insertions in a family)", y = "% insertions shared in the family")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_corr_nr_ins_p_shared_ltr.pdf", width=20, height=12)
corr_nr_ins_shared_ltr
dev.off()