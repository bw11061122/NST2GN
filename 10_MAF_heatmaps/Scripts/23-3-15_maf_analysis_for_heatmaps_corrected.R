#### [2023-03-13]
#### In this script, I am trying to get plots with minimum allele frequency for each population
#### to then plot this as a heatmap for each population 

#### Instructions:
#### minimum allele frequency
#### do a hist for SNPs
#### families with more than 10,000 insertions
#### do also all insertions for each population 
#### compare minimum allele frequency 
#### y axis = frequency 
#### x axis = number of insertions (het + 2 * homozygous)

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

### plotting
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

### sum of insertions in a given population (if doing MAF for each population separately)
## in a given family, how many insertions are present in each population?
## how many insertions = 1 * het + 2 * hom
## MAF just for new insertions
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

#### 2023-03-13
#### Pio suggested presenting MAF as a heatmap 

### average maf for insertions in a given family 
vcf_sub <- vcf[type=="MEI",c("ID", "family", "order", "superfamily", "type")]
sums_pop <- c(sum_ins_ap, sum_ins_ask, sum_ins_asm, 
              sum_ins_crhoa, sum_ins_cch, sum_ins_cvm, sum_ins_cvs, sum_ins_czeb, 
              sum_ins_diplo, sum_ins_fr, sum_ins_fuel,  sum_ins_trew, sum_ins_mzeb, 
              sum_ins_oto, sum_ins_rham)
sums_ids <- c("sum_ins_ap", "sum_ins_ask", "sum_ins_asm", 
              "sum_ins_crhoa", "sum_ins_cch", "sum_ins_cvm", "sum_ins_cvs", "sum_ins_czeb", 
              "sum_ins_diplo", "sum_ins_fr", "sum_ins_fuel", "sum_ins_trew", "sum_ins_mzeb", 
              "sum_ins_oto", "sum_ins_rham") 
maf_pop_num <- as.data.frame(cbind(sum_ins_ap, sum_ins_ask, sum_ins_asm,
                                   sum_ins_cch, sum_ins_crhoa, sum_ins_cvm, sum_ins_cvs, sum_ins_czeb, sum_ins_diplo, 
                                   sum_ins_fr, sum_ins_fuel,  sum_ins_trew, sum_ins_mzeb, 
                                   sum_ins_oto, sum_ins_rham)) / 24 ## 2 * 12 = 24 alleles in each population
maf_pop <- cbind(vcf_sub, maf_pop_num)

### sum of insertions with MAF > 0.05 for each family in each population 
### this is already normalized to the sum of all insertions
### Note: I initially normalized to length instead of sum - corrected 23-3-15
# so I want 
sum_maf_05_ap <- sum(sum_ins_ap/24 >0.05)
sum_maf_05_ask <- sum(sum_ins_ask/24 >0.05)
sum_maf_05_asm <- sum(sum_ins_asm/24 >0.05)
sum_maf_05_cch <- sum(sum_ins_cch/24 >0.05)
sum_maf_05_crhoa <- sum(sum_ins_crhoa/24 >0.05)
sum_maf_05_cvm <- sum(sum_ins_cvm/24 >0.05)
sum_maf_05_cvs <- sum(sum_ins_cvs/24 >0.05)
sum_maf_05_czeb <- sum(sum_ins_czeb/24 >0.05)
sum_maf_05_diplo <- sum(sum_ins_diplo/24 >0.05)
sum_maf_05_fr <- sum(sum_ins_fr/24 >0.05)
sum_maf_05_fuel <- sum(sum_ins_fuel/24 >0.05)
sum_maf_05_trew <- sum(sum_ins_trew/24 >0.05)
sum_maf_05_mzeb <- sum(sum_ins_mzeb/24 >0.05)
sum_maf_05_oto <- sum(sum_ins_oto/24 >0.05)
sum_maf_05_rham <- sum(sum_ins_rham/24 >0.05)

maf_05_fam_ap <- maf_pop[type=="MEI", .(maf_05_ap=sum(sum_ins_ap >0.05) / sum_maf_05_ap),by=family]
maf_05_fam_ask <- maf_pop[type=="MEI", .(maf_05_ask=sum(sum_ins_ask >0.05) / sum_maf_05_ask),by=family]
maf_05_fam_asm <- maf_pop[type=="MEI", .(maf_05_asm=sum(sum_ins_asm >0.05) / sum_maf_05_asm),by=family]
maf_05_fam_cch <- maf_pop[type=="MEI", .(maf_05_cch=sum(sum_ins_cch >0.05) / sum_maf_05_cch),by=family]
maf_05_fam_crhoa <- maf_pop[type=="MEI", .(maf_05_crhoa=sum(sum_ins_crhoa >0.05) / sum_maf_05_crhoa),by=family]
maf_05_fam_cvm <- maf_pop[type=="MEI", .(maf_05_cvm=sum(sum_ins_cvm >0.05) / sum_maf_05_cvm),by=family]
maf_05_fam_cvs <- maf_pop[type=="MEI", .(maf_05_cvs=sum(sum_ins_cvs >0.05) / sum_maf_05_cvs),by=family]
maf_05_fam_czeb <- maf_pop[type=="MEI", .(maf_05_czeb=sum(sum_ins_czeb >0.05) / sum_maf_05_czeb),by=family]
maf_05_fam_diplo <- maf_pop[type=="MEI", .(maf_05_diplo=sum(sum_ins_diplo >0.05) / sum_maf_05_diplo),by=family]
maf_05_fam_fr <- maf_pop[type=="MEI", .(maf_05_fr=sum(sum_ins_fr >0.05) / sum_maf_05_fr),by=family]
maf_05_fam_fuel <- maf_pop[type=="MEI", .(maf_05_fuel=sum(sum_ins_fuel >0.05) / sum_maf_05_fuel),by=family]
maf_05_fam_trew <- maf_pop[type=="MEI", .(maf_05_trew=sum(sum_ins_trew >0.05) / sum_maf_05_trew),by=family]
maf_05_fam_mzeb <- maf_pop[type=="MEI", .(maf_05_mzeb=sum(sum_ins_mzeb >0.05) / sum_maf_05_mzeb),by=family]
maf_05_fam_oto <- maf_pop[type=="MEI", .(maf_05_oto=sum(sum_ins_oto >0.05) / sum_maf_05_oto),by=family]
maf_05_fam_rham <- maf_pop[type=="MEI", .(maf_05_rham=sum(sum_ins_rham >0.05) / sum_maf_05_rham),by=family]

maf_05_by_pop <- cbind(maf_05_fam_ap, maf_05_fam_ask, maf_05_fam_asm, maf_05_fam_cch, maf_05_fam_crhoa,
                       maf_05_fam_cvm, maf_05_fam_cvs, maf_05_fam_czeb, maf_05_fam_diplo, maf_05_fam_fr, maf_05_fam_fuel, maf_05_fam_trew,
                       maf_05_fam_mzeb, maf_05_fam_oto, maf_05_fam_rham)
maf_05_by_pop <- maf_05_by_pop[,c(1, 2, 4, 6, 8, 10, 12, 14, 16, 18,
                                  20, 22, 24, 26, 28, 30)]
fam_names_maf <- unique(maf_pop[type=="MEI",c("family", "superfamily", "order"), with=FALSE])
maf_05_by_pop_merge <- merge(maf_05_by_pop, fam_names_maf, by="family")
maf_05_by_pop_merge <- maf_05_by_pop_merge[is.na(family==TRUE), family := "drop"]
maf_05_by_pop_merge <- maf_05_by_pop_merge[family!="drop"]

write.table(maf_05_by_pop_merge, "~/Desktop/FIGURES/23-3-15_percent_sites_05_by_pop_corrected.txt",
            sep="\t", quote=F, row.names=F)

head(meta)
meta_sum <- meta[,c("seq_depth", "genus", "species", "location", "name_loc", "pop_id"), with=FALSE]
meta_sum <- unique(meta_sum[,.(mean_seq=mean(seq_depth)), by=.(name_loc)])
meta_ord <- unique(meta[,c("location", "pop_id", "name_loc"), with=FALSE])
meta_merge <- merge(meta_sum, meta_ord, by="name_loc")
write.table(meta_merge, "~/Desktop/23-3-15_meta_table_report.txt", quote=F, row.names=F, sep="\t")
