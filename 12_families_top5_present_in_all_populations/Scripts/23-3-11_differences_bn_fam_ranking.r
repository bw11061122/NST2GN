#### [2023-03-11]
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
levels.order <- c("DNA", "LINE", "LTR", "SINE", "Unknown", "RC", "Retrotransposon", "rRNA")
col_order <- "order"
vcf[, (col_order) := factor(get(col_order), levels = levels.order)]
vcf[is.na(order==TRUE), order := "Unknown"]
upd.cols = sapply(vcf, is.factor)

## in a given family, how many insertions are present in each population?
## how many insertions = 1 * het + 2 * hom
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

### calculate how many insertions are present in ALL samples
### what % of insertions of all insertions from a given family is present in all samples
### does this match whether the family is significantly different between populations or not?
### Expectiation: families with a higher % of insertions present in all samples will not be significant 
### because these insertions are ancestral (were present in the ancestor as well)
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
vcf_sums_melt <- vcf_sums_melt[,family := as.factor(family)]
vcf_sums_melt <- vcf_sums_melt[,order := as.factor(order)]
vcf_sums_melt <- vcf_sums_melt[,superfamily := as.factor(superfamily)]

vcf_sums_melt_fam <- vcf_sums_melt[, .(sum_fam=sum(value)), by=.(family, pop_id)]
vcf_sums_melt_fam <- merge(vcf_sums_melt_fam, sum_ins_by_fam, by="pop_id")
nr_ins_pop <- vcf_sums_melt_fam[, .(sum_pop=sum(sum_fam)), by=pop_id]
vcf_sums_melt_merge <- merge(vcf_sums_melt_fam, nr_ins_pop, by="pop_id")
vcf_sums_melt_merge <- vcf_sums_melt_merge[,percent_total := sum_fam / sum_pop * 100][order(-percent_total)] # 8.99
vcf_sums_melt_merge_non <- vcf_sums_melt_merge[is.na(family==TRUE), family := "drop"]
vcf_sums_melt_merge_non <- vcf_sums_melt_merge_non[family!="drop"]

### could I assing a rank to each family and based on this, compare if the ranks between ppulations differ?
### or assing a rank to each and see if there are differences in ranks?
vcf_sums_melt_merge_non <- vcf_sums_melt_merge_non[, pop_id := as.factor(pop_id)]
vcf_sums_melt_merge_non <- vcf_sums_melt_merge_non[order(pop_id, percent_total)]
setDT(vcf_sums_melt_merge_non[, myrank := rank(-percent_total), by = pop_id])
vcf_sums_melt_merge_non <- vcf_sums_melt_merge_non[,.SD[order(myrank)], by=pop_id]

### include a plot to show the top 5 families 
vcf_sums_melt_merge_non <- vcf_sums_melt_merge_non[,pop_id := as.factor(pop_id)]
levels(vcf_sums_melt_merge_non$pop_id) <- c("Alticorpus peterdaviesi.Cape_Maclear",
"Astatotilapia calliptera.Lake_Kingiri", "Astatotilapia calliptera.Lake_Masoko", "Chilotilapia rhoadesii",
"Copadichromis chrysonotus.Lake_Malombe", "Copadichromis virginalis.Lake_Malombe", "Copadichromis virginalis.Southwest_arm", 
"Cynotilapia zebroides.Cape_Maclear", "Diplotaxodon limnothrissa.Southwest_arm", "Fossorochromis rostratus.Lake_Malombe", 
"Labeotropheus fuelleborni.Chilumba", "Labeotropheus trewavasae.Chilumba", "Maylandia zebra.Cape_Maclear", 
"Otopharynx argyrosoma.Southeast_arm", "Rhamphochromis longiceps")
vcf_sums_melt_merge_non <- vcf_sums_melt_merge_non[,pop_id := as.factor(pop_id)]
levels.pop <- c("Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
"Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
"Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
"Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
"Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", "Diplotaxodon limnothrissa.Southwest_arm")
col_meta <- "pop_id"
vcf_sums_melt_merge_non[, (col_meta) := factor(get(col_meta), levels = levels.pop)]

### add orders 
info_fam_order <- unique(vcf[type=="MEI",c("family", "superfamily", "order"), with=FALSE])
vcf_sums_melt_merge_non <- merge(vcf_sums_melt_merge_non, info_fam_order, all.x=TRUE)

top5 <- ggplot(data = vcf_sums_melt_merge_non[myrank %in% c(1:5)], aes(x = pop_id, y = percent_total, fill=order))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=cbPalette)+
  theme_bw()+
  labs(title="Major TE families across populations",x = "Population", y = "% contribution to all insertions")+
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.margin = margin(10, 10, 10, 100))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-14_top5_fam_by_pop.pdf", width=12, height=10)
top5
dev.off()

top5_mean <- vcf_sums_melt_merge_non[myrank %in% c(1:5)]
top5_mean <- top5_mean[,.(mean=mean(percent_total)), by=family]
info_fam_order <- unique(vcf[type=="MEI",c("family", "superfamily", "order"), with=FALSE])
top5_merged <- merge(top5_mean, info_fam_order, all.x=TRUE)

t5 <- ggplot(data = top5_merged, aes(x = reorder(family, c(mean)), y=mean, fill=order))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=cbPalette)+
  theme_bw()+
  labs(title="Major TE families across populations", x = "Family", y = "Average % of insertions derived from this family")+
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.margin = margin(10, 10, 10, 100))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-14_top5_fam.pdf", width=12, height=10)
t5
dev.off()

### top 10 and top 5
vcf_sums_melt_top10 <- vcf_sums_melt_merge_non[myrank %in% c(1:5)]
vcf_sums_melt_top5 <- vcf_sums_melt_merge_non[myrank %in% c(1:5)]
write.table(vcf_sums_melt_top10, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-11_top10_fam_by_pop.txt",
sep="\t", quote=F)
write.table(vcf_sums_melt_top5, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-11_top5_fam_by_pop.txt",
sep="\t", quote=F)

## add rank based on what % this family contributes to insertions in a given population
vcf_sums_melt_merge_non[, .(diff_rank = max(myrank) - min(myrank)), by=family][order(-diff_rank)]
vcf_sums_melt_ranks <- vcf_sums_melt_merge_non[,.SD[order(myrank)], by=pop_id]
vcf_sums_melt_ranks[, .(diff_rank = max(myrank) - min(myrank)), by=family][order(-diff_rank)]
vcf_rank_diff <- vcf_sums_melt_ranks[, .(diff_rank = max(myrank) - min(myrank)), by=family][order(-diff_rank)]
vcf_rank_diff_mean <- vcf_sums_melt_ranks[, .(mean_diff_rank = (max(myrank) - min(myrank)) / mean(myrank)), by=family][order(-mean_diff_rank)]

## we can check how it looks like for the top 1
vcf[family=="ERV1-8"]
vcf_sums_melt_merge_non[family=="ERV1-7"]
#################################################################################
### Families absent in specific populations
fam_absent_in_pop <- vcf_sums_melt_fam[,.(sum_absent = length(which(sum_fam=="0")), sum_tot=sum(sum_fam)), by=family][order(-sum_absent)]

#### there are 16 families which are absent in at least one population
#### one family is absent in all but 1 and two in all but two
vcf_sums_melt_fam[family=="Unknown-25",]
vcf_sums_melt_fam[family=="Unknown-27",] 
vcf_sums_melt_fam[family=="Pao-15",] 

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
