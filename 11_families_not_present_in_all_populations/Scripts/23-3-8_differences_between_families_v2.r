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

## separate significant and non-significant to check later on
res_df_sign <- res_df %>% filter(p_value < (0.05 / length(res_df$family)))
res_df_non_sign <- res_df %>% filter(p_value >= (0.05 / length(res_df$family)))

### QUESTION: is there a significant difference in the number of insertions of each family between populations,
### normalised to the total number of insertions present in a given population?

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


#### actually, on the global scale, are there differences between groups?
kruskal.test(sum_ins_fam ~ pop_id, data = vcf_sums_merged)

#### what is the distribution of p-values?
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_hist_pvalue_kruskal_norm.pdf")
hist(-log(res_df_merge$p_value), breaks=100)
dev.off()

box_ins_per_pop_norm <- ggplot(data = vcf_sums_merged, aes(x = pop_id, y = value_norm, color=pop_id))+
  geom_boxplot(outlier.shape=NA)+
  #geom_jitter()+
  ylim(c(0, 0.00001))+
  labs(title="Distribution of normalised values",x = "Population", y = "Sum of ins for a given site / total ins")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_ins_pop_id.pdf", width=12, height=12)
box_ins_per_pop_norm
dev.off()
 
hist_p_value <-ggplot(data = res_df_merge, aes(x = -log(p_value))) + 
  geom_histogram(binwidth=10)+
  geom_vline(xintercept=(-log(0.05)), color="red")+
  labs(title="Significance in Kruskal-Wallis test",x = "Adjusted p-value in the Kruskall-Wallis test", y = "Frequency")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")  
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-12_hist_pval_kruskal2.pdf", width=12, height=12)
hist_p_value
dev.off()

fam_most_sign_save <- res_df_merge %>% filter(p_value==0) %>% select(family)
write.table(fam_most_sign_save, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-12_top_fam_kw.txt",
quote=F, sep="\t")
## comparison of distirbution of normalised values for some families to visualise this
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_hist_ap.pdf")
hist(vcf_sums_merged[pop_id=="ap",value_norm])
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_hist_ask.pdf")
hist(vcf_sums_merged[pop_id=="ask",value_norm])
dev.off()
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-8_hist_asm.pdf")
hist(vcf_sums_merged[pop_id=="asm",value_norm])
dev.off()

res_df_sign_merged <- res_df_merge %>% filter(p_value < (0.05 / length(res_df_merge$family)))
res_df_non_sign_merged <- res_df_merge %>% filter(p_value >= (0.05 / length(res_df_merge$family)))
family_sign <- res_df_sign_merged %>% select(family) %>% unlist()
family_non_sign <- res_df_non_sign_merged %>% select(family) %>% unlist()

### okay, now it would be interesting to compare the significant and non-significant families
### and see whether they contributed the same or different insertions
### how exactly can I do that?

### calculate how many insertions are present in ALL samples
### what % of insertions of all insertions from a given family is present in all samples
### does this match whether the family is significantly different between populations or not?
### Expectiation: families with a higher % of insertions present in all samples will not be significant 
### because these insertions are ancestral (were present in the ancestor as well)
sum_ins_non_zero <- vcf[, apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=meta_ids]
vcf <- cbind(vcf, sum_ins_non_zero)

ins_in_all <- vcf[sum_ins_non_zero==180,] ### insertions present in all samples
ins_non_all <- vcf[sum_ins_non_zero<180,] ### insertions absent from at least one sample

nr_ins_by_fam <-vcf[, .N, by=family] ## how many insertions are there in a given family?

#### percentage of insertions which are present in all samples, out of all insertions in a given population
vcf_percent_all <- vcf[, .(percent_in_all=sum(sum_ins_non_zero==180) / .N * 100), by=family][order(-percent_in_all)]

### which families have at least one insertion whcih is present in all populations
family_ins_in_all <- vcf_percent_all[percent_in_all>0,family]
family_ins_non_all <- vcf_percent_all[percent_in_all==0,family]

common_fam_sign <- intersect(family_sign, family_ins_non_all)
common_fam_non_sign <- intersect(family_non_sign, family_ins_in_all)

### maybe determine what percentage of significant families have an insertion whcih is present in all samples?
percent_fam_sign_ins_non_all <- length(intersect(family_sign, family_ins_non_all)) / length(family_sign) * 100
percent_fam_sign_ins_in_all <- length(intersect(family_sign, family_ins_in_all)) / length(family_sign) * 100
percent_fam_non_sign_ins_non_all <- length(intersect(family_non_sign, family_ins_non_all)) / length(family_non_sign) * 100
percent_fam_non_sign_ins_in_all <- length(intersect(family_non_sign, family_ins_in_all)) / length(family_non_sign) * 100

### okay so the above is looking at TEs whcih are present in all samples
### what I want is TEs which are in at least one sample in each population
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

pop_sums <- c("sum_nz_ap", "sum_nz_ask", "sum_nz_asm", 
"sum_nz_crhoa", "sum_nz_cch", "sum_nz_cvm", "sum_nz_cvs", "sum_nz_czeb", 
"sum_nz_diplo", "sum_nz_fr", "sum_nz_fuel",  "sum_nz_trew", "sum_nz_mzeb", 
"sum_nz_oto", "sum_nz_rham") 

nr_pop_with_ins <- vcf_sums_nz[,apply(.SD, 1, function(x) length(which(x!="0"))), .SDcols=pop_sums]
### sum which are not zero in either population == there is at least one insertion in each population
vcf_sums_nz <- cbind(vcf_sums_nz, nr_pop_with_ins)

#### how many insertions in all insertions of a given family are present in each population?
vcf_percent_all_nz <- vcf_sums_nz[, .(percent_in_all=sum(nr_pop_with_ins==15) / .N * 100), by=family][order(-percent_in_all)]

### which families have at least one insertion whcih is present in all populations
family_ins_in_all_nz <- vcf_percent_all_nz[percent_in_all>0,family]
family_ins_non_all_nz <- vcf_percent_all_nz[percent_in_all==0,family]

### maybe determine what percentage of significant families have an insertion whcih is present in all samples?
### expected: significant families are less likely to be present in each population
### non-significant families are more likely to be present in each population
percent_fam_sign_ins_in_all_nz <- length(intersect(family_sign, family_ins_in_all_nz)) / length(family_sign) * 100
### what percentage of significant families have insertions present in each population?
percent_fam_sign_ins_non_all_nz <- length(intersect(family_sign, family_ins_non_all_nz)) / length(family_sign) * 100

percent_fam_non_sign_ins_in_all_nz <- length(intersect(family_non_sign, family_ins_in_all_nz)) / length(family_non_sign) * 100
### what percentage of non-significant families have insertions present in each population?
percent_fam_non_sign_ins_non_all_nz <- length(intersect(family_non_sign, family_ins_non_all_nz)) / length(family_non_sign) * 100

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

percent_ins_zero <- vcf_sums_ins_melt[, .(percent_zero=sum(value==0) / sum(value) * 100), by=.(variable, family)]
### what is this asking: what percentage of insertions in a given family are 0

percent_ins_non_zero <- vcf_sums_ins_melt[, .(percent_non_zero= 100 - (sum(value==0) / sum(value) * 100)), by=.(variable, family)][order(family)]
### what 

######### Ranking families by what % of insertions they contributed to in each population
### maybe I would like to know what % is contributed by each family in each pop to total insertions?
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

## combine sums for each population
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

vcf_sums_melt_fam <- vcf_sums_melt[, .(sum_fam=sum(value)), by=.(family, pop_id)]
vcf_sums_melt_fam <- merge(vcf_sums_melt_fam, sum_ins_by_fam, by="pop_id")
max(vcf_sums_melt_fam[,kotal]) # 8.99
which(vcf_sums_melt_fam[,percent_total>5])
vcf_sums_melt_fam_5 <- vcf_sums_melt_fam[percent_total>5,]
vcf_sums_melt_fam <- vcf_sums_melt_fam[order(pop_id, percent_total)]

### could I assing a rank to each family and based on this, compare if the ranks between ppulations differ?
### or assing a rank to each and see if there are differences in ranks?
vcf_sums_melt_fam <- vcf_sums_melt_fam[, pop_id := as.factor(pop_id)]
vcf_sums_melt_fam <- vcf_sums_melt_fam[order(pop_id, percent_total)]
vcf_sums_melt_fam <- vcf_sums_melt_fam[,rank_fam := rep(c(1:447),times=15)]
vcf_sums_melt_df <- as.data.frame(vcf_sums_melt_fam)

## add rank based on what % this family contributes to insertions in a given population
setDT(vcf_sums_melt_fam[, myrank := rank(percent_total), by = pop_id])
vcf_sums_melt_fam <- vcf_sums_melt_fam[,.SD[order(myrank)], by=pop_id]
vcf_sums_melt_fam[, .(diff_rank = max(myrank) - min(myrank)), by=family][order(-diff_rank)]
vcf_rank_diff <- vcf_sums_melt_fam[, .(diff_rank = max(myrank) - min(myrank)), by=family][order(-diff_rank)]
vcf_rank_diff_mean <- vcf_sums_melt_fam[, .(mean_diff_rank = (max(myrank) - min(myrank)) / mean(myrank)), by=family][order(-mean_diff_rank)]
vcf_rank_diff_filt <- vcf_rank_diff[diff_rank > 50,]

### okay so let's see how that looks like 
vcf_sums_melt_fam[family=="ERV1-8",]
vcf_sums_melt_fam[family=="Unknown-25",]
vcf_sums_melt_fam[family=="Unknown-58",]
vcf_sums_melt_fam[family=="Gypsy-27",]

### okay I guess I can also check which family is absent from specific populations
fam_absent_in_pop <- vcf_sums_melt_fam[,.(sum_absent = length(which(sum_fam=="0")), sum_tot=sum(sum_fam)), by=family][order(-sum_absent)]

#### there are 16 families which are absent in at least one population
#### one family is absent in all but 1 and two in all but two
vcf_sums_melt_fam[family=="Unknown-25",]
vcf_sums_melt_fam[family=="Unknown-27",] ## not very interesting - very rare insertions
vcf_sums_melt_fam[family=="Pao-15",] ## not very interesting - very rare insertions

### can filter to find those of interest w/ reasonably high frequencies
fam_absent_in_pop[sum_absent>0 & sum_tot>11,]

##### It does not make sense to do the below 
##### because you have only one datapoint (%) for each population in a given comparison
##### an alternative would be to get the percentage contributed by each insertion 
##### and then do this for each family 
vcf_sums_melt_df <- as.data.frame(vcf_sums_melt_fam)
results <- list()
for(i in levels(vcf_sums_melt_df$family)){  
  results[[i]] <- kruskal.test(formula(myrank ~ pop_id), data = vcf_sums_melt_df %>% filter(family==i))$p.value
} 

vcf_sums_melt_fam[family=="5S-Deu-L2-1",]
vcf_sums_melt_fam[family=="PiggyBac-5",]
res_df %>% filter(family=="PiggyBac-5")
kruskal.test(formula(myrank ~ pop_id), data = vcf_sums_melt_df %>% filter(family=="PiggyBac-5"))

max(res_df$p_value)
res_df <- ldply (results, data.frame)
colnames(res_df) <- c("family", "p_value")
res_df_filter <- res_df %>% filter(p_value < (0.05 / length(res_df$family)))
res_df_filter <- as.data.table(res_df_filter)[order(p_value)]
res_p_zero <- res_df_filter[p_value==0,]

### another question: in how many populations is each family present?
# Rscript 23-3-8_differences_between_families_v2.r

#### for KW test, see here: https://stackoverflow.com/questions/50194609/r-kruskal-wallis-test-in-loop-over-specified-columns-in-data-frame

#### update 23-3-11: do this for each insertion and look at distribution along the chromosome 