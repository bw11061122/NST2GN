#### [2023-03-06]
##### Trying to maybe do this smarter
##### Repeat sums but now do this no the merged file

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

#### subset meta for specific species
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

### see how the full vcf was subsetted in script 23-2-15_subset_vcf_new_meta.r
####### Modified 23-2-22 to include the filtered ones (good evidence for at least one insertion only)
## Read the subsetted vcf
vcf <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_ins_del_all_final.txt", 
    sep = "\t", header=T, fill=TRUE)

## factors  
vcf <- vcf[, order:=as.factor(order)]
vcf <- vcf[, family:=as.factor(family)]
vcf <- vcf[, superfamily:=as.factor(superfamily)]
upd.cols.vcf = sapply(vcf, is.factor)
vcf <- vcf[, names(vcf)[upd.cols.vcf] := lapply(.SD, factor), .SDcols = upd.cols.vcf] ## this gets rid of levels which are absent in the dt 

################################################################################
################################################################################
##### sum up for different orders
##### this is all operating on the sum of heterozygous sites 

## first check if there are families whcih are not present AT ALL in a given population
## how many 
sum_0_ap <- vcf[, apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_ap]
sum_1_ap <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_ap]
sum_2_ap <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_ap]
sum_0_ask <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_ask]
sum_1_ask <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_ask]
sum_2_ask <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_ask]
sum_0_asm <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_asm]
sum_1_asm <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_asm]
sum_2_asm <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_asm]
sum_0_crhoa <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_crhoa]
sum_1_crhoa <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_crhoa]
sum_2_crhoa <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_crhoa]
sum_0_cch <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_cch]
sum_1_cch <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_cch]
sum_2_cch <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_cch]
sum_0_cvm <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_cvm]
sum_1_cvm <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_cvm]
sum_2_cvm <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_cvm]
sum_0_cvs <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_cvs]
sum_1_cvs <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_cvs]
sum_2_cvs <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_cvs]
sum_0_czeb <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_czeb]
sum_1_czeb <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_czeb]
sum_2_czeb <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_czeb]
sum_0_diplo <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_diplo]
sum_1_diplo <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_diplo]
sum_2_diplo <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_diplo]
sum_0_fr <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_fr]
sum_1_fr <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_fr]
sum_2_fr <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_fr]
sum_0_fuel <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_fuel]
sum_1_fuel <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_fuel]
sum_2_fuel <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_fuel]
sum_0_trew <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_trew]
sum_1_trew <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_trew]
sum_2_trew <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_trew]
sum_0_mzeb <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_mzeb]
sum_1_mzeb <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_mzeb]
sum_2_mzeb <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_mzeb]
sum_0_oto <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_oto]
sum_1_oto <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_oto]
sum_2_oto <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_oto]
sum_0_rham <- vcf[, apply(.SD,1,function(x) length(which(x=="0"))), .SDcols=meta_rham]
sum_1_rham <- vcf[, apply(.SD,1,function(x) length(which(x=="1"))), .SDcols=meta_rham]
sum_2_rham <- vcf[, apply(.SD,1,function(x) length(which(x=="2"))), .SDcols=meta_rham]

vcf_sub <- vcf[,c("ID", "family", "order", "superfamily", "type")]
vcf_sums <- cbind(vcf_sub, sum_0_ap, sum_1_ap, sum_2_ap, 
sum_0_ask, sum_1_ask, sum_2_ask, 
sum_0_asm, sum_1_asm, sum_2_asm, 
sum_0_crhoa, sum_1_crhoa, sum_2_crhoa, 
sum_0_cch, sum_1_cch, sum_2_cch, 
sum_0_cvm, sum_1_cvm, sum_2_cvm, 
sum_0_cvs, sum_1_cvs, sum_2_cvs, 
sum_0_czeb, sum_1_czeb, sum_2_czeb, 
sum_0_diplo, sum_1_diplo, sum_2_diplo, 
sum_0_fr, sum_1_fr, sum_2_fr, 
sum_0_fuel, sum_1_fuel, sum_2_fuel, 
sum_0_trew, sum_1_trew, sum_2_trew, 
sum_0_mzeb, sum_1_mzeb, sum_2_mzeb, 
sum_0_oto, sum_1_oto, sum_2_oto, 
sum_0_rham, sum_1_rham, sum_2_rham) 
vcf_sums_melt <- melt(vcf_sums)[order(ID)]

vcf_sums_melt[, c("del", "state", "pop_id") := tstrsplit(variable, "_", fixed=TRUE)]
vcf_sums_melt[,del := NULL]

vcf_sums_melt_fam <- vcf_sums_melt[, .(mean_fam=mean(value)), by=.(family, state)]
vcf_sums_melt_pop <- vcf_sums_melt[, .(mean_pop=mean(value)), by=.(pop_id, state)]
vcf_sums_melt_fam_pop <- vcf_sums_melt[, .(mean_fam_pop = mean(value)), by=.(pop_id, family, state)][order(family)]
vcf_sums_group_state <- vcf_sums_melt[, .(mean_fam=mean(value)), by=.(family, state)]

melt_ins <- vcf_sums_melt[type=="MEI",]
melt_del <- vcf_sums_melt[type=="MEA",]
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_hist_by_genotype_ins.pdf")
ggplot(melt_ins,aes(x=value))+geom_histogram()+facet_grid(~state)+theme_bw()
dev.off()

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_hist_by_genotype_del.pdf")
ggplot(melt_del,aes(x=value))+geom_histogram()+facet_grid(~state)+theme_bw()
dev.off()

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_hist_by_genotype.pdf")
ggplot(vcf_sums_melt,aes(x=value))+geom_histogram()+facet_grid(~state)+theme_bw()
dev.off()

#### Now onto Kruskall-Wallis test
#### this will be to normalise the values 
vcf_sums_melt_pop <- vcf_sums_melt[, .(mean_pop=mean(value)), by=.(pop_id, state)]
kruskal.test(mean_pop ~ pop_id, data = vcf_sums_melt_pop) 
## this says there are no stat significant differences actually between values in these families

vcf_sums_merged <- merge(vcf_sums_melt, vcf_sums_melt_pop, by=c("state", "pop_id"))[order(family)]

kruskal.test(mean_pop ~ pop_id, data = vcf_sums_melt_pop) 

library(dplyr)
group_by((vcf_sums_melt), family) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
  )

kruskal.test(value ~ pop_id, data = vcf_sums_melt)

vcf_sums_merged <- vcf_sums_merged[, value_norm := value / mean_pop]
kruskal.test(value_norm ~ pop_id, data = vcf_sums_merged)
kruskal.test(value_norm ~ family, data = vcf_sums_merged)

# data:  value_norm by pop_id
# Kruskal-Wallis chi-squared = 275479, df = 14, p-value < 2.2e-16

# data:  value_norm by family
# Kruskal-Wallis chi-squared = 40056, df = 446, p-value < 2.2e-16
#### okay so by this it is significant although I don't know if I'm actually doing it right tbh

## okay can I now find families which are significantly differen?
sum_by_fam_pop <- vcf_sums_melt[,.(sum_state=sum(value)), by = .(family, pop_id, state)]
kruskal.test(sum_state ~ pop_id, data = sum_by_fam_pop) ## no difference between populations
kruskal.test(sum_state ~ family, data = sum_by_fam_pop)

### filter out where value is 0 for state 1 and value is 0 for state 2
sum_by_fam_pop_absent_hom <- sum_by_fam_pop[(state=="2" & sum_state=="0"),]
sum_by_fam_pop_absent_het <- sum_by_fam_pop[(state=="1" & sum_state=="0"),]

### how does this look like for PiggyBac
sum_by_fam_pb <- sum_by_fam_pop %>% as.data.frame %>% filter(grepl('PiggyBac', family))

#### This is what I've done 2023-03-06
###############################################################################
###############################################################################
###############################################################################
## so I dont think there are any families which are private to a specific population
## there still may be insertions which are private to some 

box_TEs_per_sample <- ggplot(data = vcf_sums_melt, aes(x = value, y = sum_pol_per_sample, color=order))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  ylim(c(0,65000))+
  labs(title="Sum of polymorphisms per sample",x = "Population", y = "Sum of polymorphisms per sample")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-2_MEI_sum_polymorphisms_per_sample.pdf", width=18, height=14)
box_TEs_per_sample
dev.off()

######################################################################################################
######################################################################################################
####### 2023-03-06 and 2023-03-08
####### This is what I was trying to get done previously
####### But I am not sure if this made any sense 
## for each family, superfamily and order, I am calculating how many het sites are there in each sample
sum_ins_by_order <- vcf[, lapply(.SD, function(x) length(which(x=="1"))+2*length(which(x=="2"))), by=order, .SDcols=meta_ids]
sum_ins_by_superfamily <- vcf[, lapply(.SD, function(x) length(which(x=="1"))+2*length(which(x=="2"))), by=superfamily, .SDcols=meta_ids]
sum_ins_by_family <- vcf[, lapply(.SD,function(x) length(which(x=="1"))+2*length(which(x=="2"))), by=family, .SDcols=meta_ids]

## now I am converting it to a data frame
t_sum_by_order <- as.data.frame(t(sum_ins_by_order[,-1]))
colnames(t_sum_by_order) <- sum_ins_by_order$order
t_sum_by_superfamily <- as.data.frame(t(sum_ins_by_superfamily[,-1]))
colnames(t_sum_by_superfamily) <- sum_ins_by_superfamily$superfamily
t_sum_by_family <- as.data.frame(t(sum_ins_by_family[,-1]))
colnames(t_sum_by_family) <- sum_ins_by_family$family

#### this gives you for each cichlid, how many insertions of each order there are 
meta_sub <- meta[,c("pop_id", "clade"), with=FALSE]
meta_sums_family <- cbind(meta_sub, t_sum_by_family)
meta_sums_superfamily <- cbind(meta_sub, t_sum_by_superfamily)
meta_sums_order <- cbind(meta_sub, t_sum_by_order)
meta_sums_order_melt <- melt(meta_sums_order)
box_TEs_per_sample_per_order <- ggplot(data = meta_sums_order_melt, aes(x = variable, y = value))+
  geom_boxplot()+
  geom_jitter(aes(color=pop_id, size=5))+
  theme_bw()+
  # facet_wrap(~population)+
  labs(title="Number of polymorphic sites per sample",x = "TE order", y = "TE sites per sample")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_MEI_sum_ins_per_sample_per_order.pdf", width=22, height=16)
box_TEs_per_sample_per_order
dev.off()

box_TEs_per_sample <- ggplot(data = meta_sums_order_melt, aes(x = pop_id, y = value))+
  geom_boxplot()+
  geom_jitter(aes(color=variable, size=5))+
  theme_bw()+
  # facet_wrap(~population)+
  labs(title="Number of polymorphic sites per sample",x = "TE order", y = "TE sites per sample")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_MEI_sum_ins_per_sample.pdf", width=22, height=16)
box_TEs_per_sample
dev.off()

meta_sums_family <- cbind(meta_sub, t_sum_by_family)
meta_sums_family_melt <- melt(meta_sums_family)
meta_family <- meta_sums_family_melt[, mean_pop := mean(value), by=pop_id] ### how many insertions are there 
meta_family <- meta_sums_family_melt[, mean_fam := mean(value), by=variable]
meta_family <- meta_sums_family_melt[, mean_fam_pop := mean(value), by=.(variable, pop_id)]
meta_family_mean <- unique(meta_family[,c("clade", "pop_id", "variable", "mean_fam", "mean_pop", "mean_fam_pop"), with=FALSE])

### calculate absolute difference 
mean_diff_fam <- meta_family_mean[, diff := mean_fam_pop - mean_fam]

### calculate percentage difference 
mean_diff_fam <- meta_family_mean[, diff_percent := (mean_fam_pop - mean_fam) / mean_fam * 100]

mean_diff_fam[, c("superfamily", "number") := tstrsplit(variable, "-", fixed=FALSE, keep=c(1,2))]
mean_diff_fam <- mean_diff_fam[,superfamily := as.factor(superfamily)]
mean_diff_fam[order(-diff_percent)]
mean_piggy <- mean_diff_fam[superfamily == "PiggyBac"]

box_TEs_per_family <- ggplot(data = mean_piggy, aes(x = pop_id, y = diff_percent))+
  geom_point(aes(color=variable, size=5))+
  theme_bw()+
  facet_wrap(~variable)+
  labs(title="Distribution of PiggyBack families",x = "Population", y = "% number of mean insertions in the family relative to average number for the family")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_MEI_percent_diff_family_piggyback.pdf", width=22, height=16)
box_TEs_per_family
dev.off()

### Is there a way I could identify families which are only present in a subset of families?
mean_diff_fam[mean_fam_pop==0,] ## this will return all families which are zero in the specific population they are absent from

### based on this, can I find families with a large percentage change?
mean_diff_pc <- mean_diff_fam[diff_percent > 100 & mean_fam > 5,]

#####################################################################
#####################################################################
## in a given family, how many insertions are present in each population?
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
### then, I will want to match this to familie and do separate plots by families 

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
vcf_sums_melt_fam <- vcf_sums_melt[, .(mean_fam=mean(value)), by=.(family, state)]
vcf_sums_melt_pop <- vcf_sums_melt[, .(mean_pop=mean(value)), by=.(pop_id, state)]
vcf_sums_melt_fam_pop <- vcf_sums_melt[, .(mean_fam_pop = mean(value)), by=.(pop_id, family, state)][order(family)]
vcf_sums_group_state <- vcf_sums_melt[, .(mean_fam=mean(value)), by=.(family, state)]

kruskal.test(sum_state ~ pop_id, data = sum_by_fam_pop, p.adj='bonferroni')

#####################################################################
#####################################################################
### 2023-03-07 New Kruskall-Wallis
### we eventually decided it did not work 
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
vcf_sums_melt_fam <- vcf_sums_melt[, .(mean_fam=mean(value)), by=.(family, state)]
vcf_sums_melt_pop <- vcf_sums_melt[, .(mean_pop=mean(value)), by=.(pop_id, state)]
vcf_sums_melt_fam_pop <- vcf_sums_melt[, .(mean_fam_pop = mean(value)), by=.(pop_id, family, state)][order(family)]
vcf_sums_group_state <- vcf_sums_melt[, .(mean_fam=mean(value)), by=.(family, state)]

kruskal.test(sum_state ~ pop_id, data = sum_by_fam_pop, p.adj='bonferroni')

# Rscript 23-3-6_differences_between_families_v2.r
