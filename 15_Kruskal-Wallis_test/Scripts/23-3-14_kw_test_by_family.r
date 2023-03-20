#### [2023-03-08]
##### 

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

## plotting
cbPalette <- c("#EF5350", "#2196F3", "#43A047", "#9575CD", "#BDBDBD", "#D1C4E9", "#757575", "#795548") # orders
popPalette <- c("#9CCC65", "#7CB342", "#7E57C2", "512DA8", "#AB47BC", "#6A1B9A", "#1B5E20", "8BC34A", "388E3C",
"#B71C1C", "#E57373", "#D32F2F", "#1A237E", "#BF360C", "#827717")
cladePalette <- c("#AED581", "#5E35B1", "#1B5E20", "#B71C1C", "#1A237E", "#EF6C00", "#827717")

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

### order factor variables
levels.pop <- c("Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
"Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
"Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
"Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
"Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", "Diplotaxodon limnothrissa.Southwest_arm")
col <- "pop_id"
meta[, (col) := factor(get(col), levels = levels.pop)]
levels.clade <- c("AstCal", "Mbuna", "Utaka", "Benthic", "Deep","Rhampho", "Diplo")
col_clade <- "clade"
meta[, (col_clade) := factor(get(col_clade), levels = levels.clade)]

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

#####################################################################
#####################################################################
## in a given family, how many insertions are present in each population?
## how many insertions = 1 * het + 2 * hom
## plan: do KW test
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

### so this is summing up each insertion for each of the 15 populations that I have
### then, I will want to match this to families and do separate plots by families 

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

### KW test for each family
box_PB5 <- ggplot(data = vcf_sums_melt_piggyback, aes(x = pop_id, y = value, color=pop_id))+
  geom_boxplot(outlier.shape=NA)+
  ylim(c(0,24))+
  #geom_jitter()+
  labs(title="Distirbution of insertions in the PiggyBac-5 family",x = "Population", y = "Sum of polymorphisms per sample")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_PiggyBac5_trial_comaprison.pdf", width=12, height=12)
box_PB5
dev.off()

#### QUESTION: is there a significant difference in the number of inseertions of each family between populations?
vcf_melt_df <- as.data.frame(vcf_sums_melt)
vcf_melt_df$family <- as.factor(vcf_melt_df$family)
results <- list()
for(i in levels(vcf_melt_df$family)){  
  results[[i]] <- kruskal.test(formula(value ~ pop_id), data = vcf_melt_df %>% filter(family==i))$p.value
} 

res_df <- ldply (results, data.frame)
colnames(res_df) <- c("family", "p_value")
res_df_filter <- res_df %>% filter(p_value < (0.05 / length(res_df$family)))
res_df_filter <- as.data.table(res_df_filter)[order(p_value)]
res_p_zero <- res_df_filter[p_value==0,]

### on values normalised to population average
sum_ins_by_fam <- vcf_sums_melt[, .(sum_ins_fam=sum(value)), by=pop_id]
### this is the same if you do sum(unlist(sum_ins_pop for any given population))
### so basically gives you the total number of insertions in the population

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
res_df_filter_merge <- as.data.table(res_df_filter_merge)[order(p_value)] ## worryingly this has even more significant ones
res_p_zero_merge <- res_df_filter_merge[p_value==0,] 

vcf_merg_df$superfamily <- as.factor(vcf_merg_df$superfamily)
results_sfam <- list()
for(i in levels(vcf_merg_df$superfamily)){  
  results_sfam[[i]] <- kruskal.test(formula(value_norm ~ pop_id), data = vcf_merg_df %>% filter(superfamily==i))$p.value
} 
res_df_merge <- ldply (results_merged, data.frame)
colnames(res_df_merge) <- c("family", "p_value")
res_df_filter_merge <- res_df_merge %>% filter(p_value < (0.05 / length(res_df_merge$family)))
res_df_filter_merge <- as.data.table(res_df_filter_merge)[order(p_value)] ## worryingly this has even more significant ones
res_p_zero_merge <- res_df_filter_merge[p_value==0,] 

hist_p_value <-ggplot(data = res_df_merge, aes(x = -log(p_value))) + 
  geom_histogram(binwidth=10)+
  geom_vline(xintercept=(-log(0.05)), color="red")+
  labs(title="Significance in Kruskal-Wallis test",x = "-log(adjusted p-value)", y = "Frequency")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")  
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-14_hist_pval_kruskal.pdf", width=12, height=12)
hist_p_value
dev.off()



