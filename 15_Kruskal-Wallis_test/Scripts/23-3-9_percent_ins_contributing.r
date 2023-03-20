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

##########################################################################################################
##########################################################################################################
#### Another way of looking at this 
#### for each insertion, sum it up for each population
sum_ins_ap <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ap]
sum_ins_ask <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ask]
sum_ins_asm <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_asm]
sum_ins_cch <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cch]
sum_ins_crhoa <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_crhoa]
sum_ins_cvm <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cvs]
sum_ins_cvs <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cvm]
sum_ins_czeb <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_czeb]
sum_ins_diplo <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_diplo]
sum_ins_fr <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_fr]
sum_ins_fuel <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_fuel]
sum_ins_trew <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_trew]
sum_ins_mzeb <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_mzeb]
sum_ins_oto <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_oto]
sum_ins_rham <- vcf[, apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_rham]

vcf_sub <- vcf[,c("ID", "family", "order", "superfamily", "type")]
vcf_sums_ins <- cbind(vcf_sub, sum_ins_ap, sum_ins_ask, sum_ins_asm, 
sum_ins_crhoa, sum_ins_cch, sum_ins_cvm, sum_ins_cvs, sum_ins_czeb, 
sum_ins_diplo, sum_ins_fr, sum_ins_fuel,  sum_ins_trew, sum_ins_mzeb, 
sum_ins_oto, sum_ins_rham) 
vcf_sums_ins_melt <- melt(vcf_sums_ins)

percent_ins_zero_len <- vcf_sums_ins_melt[, .(percent_zero=length(which((value==0))) / length(value) * 100), by=.(variable, family)][order(-percent_zero)]
percent_ins_non_zero_len <- vcf_sums_ins_melt[, .(percent_non_zero=100-(length(which((value==0))) / length(value) * 100)), by=.(variable, family)][order(percent_non_zero)]

# The arguments to spread():
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
library(tidyr)
ins_zero_sum_wide <- spread(percent_ins_zero_len, variable, percent_zero)
ins_non_zero_sum_wide <- spread(percent_ins_non_zero_len, variable, percent_non_zero)

#### Kruskall-Wallis test to correlate the percentage with the p-value
## okay so Kruskall on (sum_hets + 2*sum_homs) / nr fish ~ population for each family
## Bonferroni is 1/(nr tests)
vcf_sub <- vcf[,c("ID", "family", "order", "superfamily", "type")]
vcf_sums <- cbind(vcf_sub, sum_ins_ap, sum_ins_ask, sum_ins_asm, 
sum_ins_crhoa, sum_ins_cch, sum_ins_cvm, sum_ins_cvs, sum_ins_czeb, 
sum_ins_diplo, sum_ins_fr, sum_ins_fuel,  sum_ins_trew, sum_ins_mzeb, 
sum_ins_oto, sum_ins_rham) 
vcf_sums_melt <- melt(vcf_sums)[order(ID)]
vcf_sums_melt[, c("del", "del2", "pop_id") := tstrsplit(variable, "_", fixed=TRUE)]
vcf_sums_melt[,c("del", "del2") := NULL]
vcf_sums_melt_piggyback <- vcf_sums_melt[family=="PiggyBac-5",]
vcf_sums_melt <- vcf_sums_melt[,family := as.factor(family)]
vcf_sums_melt <- vcf_sums_melt[,order := as.factor(order)]
vcf_sums_melt <- vcf_sums_melt[,superfamily := as.factor(superfamily)]

### QUESTION: is there a significant difference in the number of insertions of each family between populations,
### normalised to the total number of insertions present in a given population?

sum_ins_by_fam <- vcf_sums_melt[, .(sum_ins_fam=sum(value)), by=pop_id]

vcf_sums_merged <- merge(vcf_sums_melt, sum_ins_by_fam, by="pop_id")
vcf_sums_merged <- vcf_sums_merged[,value_norm := value / sum_ins_fam] ## divide by the total number of insertions

vcf_merg_df <- as.data.frame(vcf_sums_merged)
vcf_merg_df$family <- as.factor(vcf_merg_df$family)
results_merged <- list()
for(i in levels(vcf_merg_df$family)){  
  results_merged[[i]] <- kruskal.test(formula(value_norm ~ pop_id), data = vcf_merg_df %>% filter(family==i))$p.value
} 
res_df_merge <- ldply (results_merged, data.frame)
colnames(res_df_merge) <- c("family", "p_value")
res_df_filter_merge <- res_df_merge %>% filter(p_value < (0.05 / length(res_df_merge$family)))
res_df_filter_merge <- as.data.table(res_df_filter_merge)[order(p_value)] 
res_p_zero_merge <- res_df_filter_merge[p_value==0,] 

res_dt_merge <- as.data.table(res_df_merge)
ins_zero_sum_dt <- as.data.table(ins_zero_sum_wide)
ins_zero_sum_merge <- merge(ins_zero_sum_dt, res_dt_merge, by="family")  ### what percentage of insertions in a given family are 0
ins_non_zero_sum_dt <- as.data.table(ins_non_zero_sum_wide)
ins_non_zero_sum_merge <- merge(ins_non_zero_sum_dt, res_dt_merge, by="family") ## for each family, this is the % of insertisons which DO contribute 

### now find out the total number of insertions for a given family and merge (number of ins per family for context)
nr_ins_fam <- vcf[,.N, by=family]
ins_zero_sum_merge <- merge(ins_zero_sum_merge, nr_ins_fam, by="family")
ins_non_zero_sum_merge <- merge(ins_non_zero_sum_merge, nr_ins_fam, by="family")

percent_ins_sign <- ins_zero_sum_merge[p_value==0,] ## 28 families 
percent_ins_non_sign <- ins_zero_sum_merge[p_value>0.994,] ## p_value = 0.994 gives the same number of families (28) 

percent_ins_nz_sign <- ins_non_zero_sum_merge[p_value==0,]
percent_ins_nz_non_sign <- ins_non_zero_sum_merge[p_value>0.994,]

### extra stuff: is there a correlation between the significance and the number of insertions identified
res_df_nr <- merge(res_df_merge, nr_ins_fam, by="family")
res_df_plot <- res_df_nr %>% filter(p_value!=0)
corr_pval_nr <- ggplot(data = res_df_plot, aes(x = log(N), y = -log(p_value)))+
  geom_point()+
  ylim(c(0,1000))+
  labs(title="Correlation between number of insertions in a family and significance in Kruskal test",x = "log(number of insertions in a family)", y = "-log(adj P value)")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_corr_nr_ins_fam_p_val.pdf", width=20, height=12)
corr_pval_nr
dev.off()
