#### [2023-3-11]
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

## set up colors
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
levels.order <- c("DNA", "LINE", "LTR", "SINE", "Unknown", "RC", "Retrotransposon", "rRNA")
col <- "order"
vcf_sum_order[, (col) := factor(get(col), levels = levels.order)]

##### meta for specific cichlids
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

#####################################################################################################
#####################################################################################################
#####################################################################################################
#### How many sites of all identified are 0, 0.5, 1, 1.5, 2?
sum_ins_0_ap <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_ap]
sum_ins_05_ap <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_ap]
sum_ins_1_ap <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_ap]
sum_ins_15_ap <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_ap]
sum_ins_2_ap <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_ap]
sum_ins_0_ask <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_ask]
sum_ins_05_ask <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_ask]
sum_ins_1_ask <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_ask]
sum_ins_15_ask <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_ask]
sum_ins_2_ask <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_ask]
sum_ins_0_asm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_asm]
sum_ins_05_asm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_asm]
sum_ins_1_asm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_asm]
sum_ins_15_asm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_asm]
sum_ins_2_asm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_asm]
sum_ins_0_cch <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_cch]
sum_ins_05_cch <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_cch]
sum_ins_1_cch <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_cch]
sum_ins_15_cch <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_cch]
sum_ins_2_cch <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_cch]
sum_ins_0_crhoa <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_crhoa]
sum_ins_05_crhoa <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_crhoa]
sum_ins_1_crhoa <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_crhoa]
sum_ins_15_crhoa <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_crhoa]
sum_ins_2_crhoa <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_crhoa]
sum_ins_0_cvm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_cvm]
sum_ins_05_cvm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_cvm]
sum_ins_1_cvm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_cvm]
sum_ins_15_cvm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_cvm]
sum_ins_2_cvm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_cvm]
sum_ins_0_cvs <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_cvs]
sum_ins_05_cvs <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_cvs]
sum_ins_1_cvs <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_cvs]
sum_ins_15_cvs <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_cvs]
sum_ins_2_cvs <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_cvs]
sum_ins_0_czeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_czeb]
sum_ins_05_czeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_czeb]
sum_ins_1_czeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_czeb]
sum_ins_15_czeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_czeb]
sum_ins_2_czeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_czeb]
sum_ins_0_diplo <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_diplo]
sum_ins_05_diplo <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_diplo]
sum_ins_1_diplo <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_diplo]
sum_ins_15_diplo <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_diplo]
sum_ins_2_diplo <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_diplo]
sum_ins_0_fr <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_fr]
sum_ins_05_fr <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_fr]
sum_ins_1_fr <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_fr]
sum_ins_15_fr <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_fr]
sum_ins_2_fr <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_fr]
sum_ins_0_fuel <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_fuel]
sum_ins_05_fuel <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_fuel]
sum_ins_1_fuel <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_fuel]
sum_ins_15_fuel <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_fuel]
sum_ins_2_fuel <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_fuel]
sum_ins_0_trew <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_trew]
sum_ins_05_trew <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_trew]
sum_ins_1_trew <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_trew]
sum_ins_15_trew <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_trew]
sum_ins_2_trew <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_trew]
sum_ins_0_mzeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_mzeb]
sum_ins_05_mzeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_mzeb]
sum_ins_1_mzeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_mzeb]
sum_ins_15_mzeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_mzeb]
sum_ins_2_mzeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_mzeb]
sum_ins_0_oto <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_oto]
sum_ins_05_oto <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_oto]
sum_ins_1_oto <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_oto]
sum_ins_15_oto <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_oto]
sum_ins_2_oto <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_oto]
sum_ins_0_rham <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_rham]
sum_ins_05_rham <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_rham]
sum_ins_1_rham <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_rham]
sum_ins_15_rham <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_rham]
sum_ins_2_rham <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_rham]

sum_ap <- cbind(sum_ins_0_ap, sum_ins_05_ap, sum_ins_1_ap, sum_ins_15_ap, sum_ins_2_ap)
sum_ask <- cbind(sum_ins_0_ask, sum_ins_05_ask, sum_ins_1_ask, sum_ins_15_ask, sum_ins_2_ask)
sum_asm <- cbind(sum_ins_0_asm, sum_ins_05_asm, sum_ins_1_asm, sum_ins_15_asm, sum_ins_2_asm)
sum_cch <- cbind(sum_ins_0_cch, sum_ins_05_cch, sum_ins_1_cch, sum_ins_15_cch, sum_ins_2_cch)
sum_crhoa <- cbind(sum_ins_0_crhoa, sum_ins_05_crhoa, sum_ins_1_crhoa, sum_ins_15_crhoa, sum_ins_2_crhoa)
sum_cvm <- cbind(sum_ins_0_cvm, sum_ins_05_cvm, sum_ins_1_cvm, sum_ins_15_cvm, sum_ins_2_cvm)
sum_cvs <- cbind(sum_ins_0_cvs, sum_ins_05_cvs, sum_ins_1_cvs, sum_ins_15_cvs, sum_ins_2_cvs)
sum_czeb <- cbind(sum_ins_0_czeb, sum_ins_05_czeb, sum_ins_1_czeb, sum_ins_15_czeb, sum_ins_2_czeb)
sum_diplo <- cbind(sum_ins_0_diplo, sum_ins_05_diplo, sum_ins_1_diplo, sum_ins_15_diplo, sum_ins_2_diplo)
sum_fr <- cbind(sum_ins_0_fr, sum_ins_05_fr, sum_ins_1_fr, sum_ins_15_fr, sum_ins_2_fr)
sum_fuel <- cbind(sum_ins_0_fuel, sum_ins_05_fuel, sum_ins_1_fuel, sum_ins_15_fuel, sum_ins_2_fuel)
sum_trew <- cbind(sum_ins_0_trew, sum_ins_05_trew, sum_ins_1_trew, sum_ins_15_trew, sum_ins_2_trew)
sum_mzeb <- cbind(sum_ins_0_mzeb, sum_ins_05_mzeb, sum_ins_1_mzeb, sum_ins_15_mzeb, sum_ins_2_mzeb)
sum_oto <- cbind(sum_ins_0_oto, sum_ins_05_oto, sum_ins_1_oto, sum_ins_15_oto, sum_ins_2_oto)
sum_rham <- cbind(sum_ins_0_rham, sum_ins_05_rham, sum_ins_1_rham, sum_ins_15_rham, sum_ins_2_rham)
sums <- cbind(sum_ap, sum_ask, sum_asm, sum_cch, sum_crhoa, sum_cvm, sum_cvs, sum_czeb, sum_diplo, sum_fr, sum_fuel,
sum_trew, sum_mzeb, sum_oto, sum_rham)
count_sum_states <- colSums(sums) ## this will be for all insertions
count_sum_states <- as.data.frame(count_sum_states) %>% rownames_to_column("variable")
count_sum_dt <- as.data.table(count_sum_states)
count_sum_dt[, c("del", "del2","state", "pop_id") := tstrsplit(variable, "_", fixed=TRUE)]
count_sum_dt[,c("del", "del2") := NULL]
count_sum_dt[,state := as.factor(state)]
count_sum_dt[,pop_id := as.factor(pop_id)]

count_sum_dt[, state:=`levels<-`(state, c("0", "0.5", "1", "1.5", "2"))]
count_sum_dt[, pop_id:=`levels<-`(pop_id, c("Alticorpus peterdaviesi.Cape_Maclear", "Astatotilapia calliptera.Lake_Kingiri", 
"Astatotilapia calliptera.Lake_Masoko", "Copadichromis chrysonotus.Lake_Malombe", "Chilotilapia rhoadesii",
"Copadichromis virginalis.Lake_Malombe","Copadichromis virginalis.Southwest_arm", "Cynotilapia zebroides.Cape_Maclear",
"Diplotaxodon limnothrissa.Southwest_arm", "Fossorochromis rostratus.Lake_Malombe", "Labeotropheus fuelleborni.Chilumba",
"Labeotropheus trewavasae.Chilumba", "Maylandia zebra.Cape_Maclear", "Otopharynx argyrosoma.Southeast_arm",
"Rhamphochromis longiceps"))]
levels.pop <- c("Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
"Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
"Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
"Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
"Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", "Diplotaxodon limnothrissa.Southwest_arm")
col <- "pop_id"
count_sum_dt[, (col) := factor(get(col), levels = levels.pop)]

meta_clade <- unique(meta[,c("pop_id", "clade"), with=FALSE])
count_sum_merged <- merge(count_sum_dt, meta_clade, by="pop_id")
levels.clade <- c("AstCal", "Mbuna", "Utaka", "Benthic", "Deep","Rhampho", "Diplo")
col_clade <- "clade"
count_sum_merged[, (col_clade) := factor(get(col_clade), levels = levels.clade)]
setnames(count_sum_merged, "pop_id", "population")
ins_by_state <- ggplot(data = count_sum_merged, aes(x = population, y = count_sum_states, fill=state))+
  geom_col(aes(fill=state))+
  scale_fill_manual(values=c("#FFAB91", "#FF7043", "#F4511E", "#D84315", "#BF360C"))+
  labs(title="Distribution of insertions in each population",x = "Population", y = "Number of insertions")+
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.margin = margin(10, 10, 10, 100))
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_nr_ins_states_by_pop.pdf", width=22, height=10)
ins_by_state
dev.off()

## do I need to care about 1.5?
total_ins_by_pop <- count_sum_merged[, .(sum_total=sum(count_sum_states)), by=population]
count_0 <- count_sum_merged[state==0,]
count_0 <- count_0[,percent := count_sum_states / 5969184 * 100]
count_1 <- count_sum_merged[state==1,]
count_1 <- count_1[,percent := count_sum_states / 5969184 * 100]
count_2 <- count_sum_merged[state==2,]
count_2 <- count_2[,percent := count_sum_states / 5969184 * 100]
count_05 <- count_sum_merged[state==0.5,]
count_05 <- count_05[,percent := count_sum_states / 5969184 * 100]
count_15 <- count_sum_merged[state==1.5,]
count_15 <- count_15[,percent := count_sum_states / 5969184 * 100][order(percent)]
counts_percent <- rbind(count_0, count_05, count_1, count_15, count_2)[order(percent)]
## I would say there were so few I did not account for this further on

################################################################################################
#### Same for deletions
sum_del_0_ap <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_ap]
sum_del_05_ap <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_ap]
sum_del_1_ap <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_ap]
sum_del_15_ap <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_ap]
sum_del_2_ap <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_ap]
sum_del_0_ask <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_ask]
sum_del_05_ask <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_ask]
sum_del_1_ask <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_ask]
sum_del_15_ask <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_ask]
sum_del_2_ask <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_ask]
sum_del_0_asm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_asm]
sum_del_05_asm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_asm]
sum_del_1_asm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_asm]
sum_del_15_asm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_asm]
sum_del_2_asm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_asm]
sum_del_0_cch <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_cch]
sum_del_05_cch <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_cch]
sum_del_1_cch <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_cch]
sum_del_15_cch <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_cch]
sum_del_2_cch <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_cch]
sum_del_0_crhoa <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_crhoa]
sum_del_05_crhoa <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_crhoa]
sum_del_1_crhoa <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_crhoa]
sum_del_15_crhoa <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_crhoa]
sum_del_2_crhoa <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_crhoa]
sum_del_0_cvm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_cvm]
sum_del_05_cvm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_cvm]
sum_del_1_cvm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_cvm]
sum_del_15_cvm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_cvm]
sum_del_2_cvm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_cvm]
sum_del_0_cvs <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_cvs]
sum_del_05_cvs <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_cvs]
sum_del_1_cvs <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_cvs]
sum_del_15_cvs <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_cvs]
sum_del_2_cvs <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_cvs]
sum_del_0_czeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_czeb]
sum_del_05_czeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_czeb]
sum_del_1_czeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_czeb]
sum_del_15_czeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_czeb]
sum_del_2_czeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_czeb]
sum_del_0_diplo <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_diplo]
sum_del_05_diplo <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_diplo]
sum_del_1_diplo <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_diplo]
sum_del_15_diplo <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_diplo]
sum_del_2_diplo <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_diplo]
sum_del_0_fr <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_fr]
sum_del_05_fr <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_fr]
sum_del_1_fr <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_fr]
sum_del_15_fr <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_fr]
sum_del_2_fr <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_fr]
sum_del_0_fuel <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_fuel]
sum_del_05_fuel <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_fuel]
sum_del_1_fuel <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_fuel]
sum_del_15_fuel <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_fuel]
sum_del_2_fuel <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_fuel]
sum_del_0_trew <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_trew]
sum_del_05_trew <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_trew]
sum_del_1_trew <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_trew]
sum_del_15_trew <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_trew]
sum_del_2_trew <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_trew]
sum_del_0_mzeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_mzeb]
sum_del_05_mzeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_mzeb]
sum_del_1_mzeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_mzeb]
sum_del_15_mzeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_mzeb]
sum_del_2_mzeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_mzeb]
sum_del_0_oto <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_oto]
sum_del_05_oto <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_oto]
sum_del_1_oto <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_oto]
sum_del_15_oto <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_oto]
sum_del_2_oto <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_oto]
sum_del_0_rham <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_rham]
sum_del_05_rham <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_rham]
sum_del_1_rham <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_rham]
sum_del_15_rham <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_rham]
sum_del_2_rham <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_rham]

sum_del <- cbind(sum_del_0_ap, sum_del_05_ap, sum_del_1_ap, sum_del_15_ap, sum_del_2_ap,
sum_del_0_ask, sum_del_05_ask, sum_del_1_ask, sum_del_15_ask, sum_del_2_ask,
sum_del_0_asm, sum_del_05_asm, sum_del_1_asm, sum_del_15_asm, sum_del_2_asm,
sum_del_0_cch, sum_del_05_cch, sum_del_1_cch, sum_del_15_cch, sum_del_2_cch,
sum_del_0_crhoa, sum_del_05_crhoa, sum_del_1_crhoa, sum_del_15_crhoa, sum_del_2_crhoa,
sum_del_0_cvm, sum_del_05_cvm, sum_del_1_cvm, sum_del_15_cvm, sum_del_2_cvm,
sum_del_0_cvs, sum_del_05_cvs, sum_del_1_cvs, sum_del_15_cvs, sum_del_2_cvs,
sum_del_0_czeb, sum_del_05_czeb, sum_del_1_czeb, sum_del_15_czeb, sum_del_2_czeb,
sum_del_0_diplo, sum_del_05_diplo, sum_del_1_diplo, sum_del_15_diplo, sum_del_2_diplo,
sum_del_0_fr, sum_del_05_fr, sum_del_1_fr, sum_del_15_fr, sum_del_2_fr,
sum_del_0_fuel, sum_del_05_fuel, sum_del_1_fuel, sum_del_15_fuel, sum_del_2_fuel,
sum_del_0_trew, sum_del_05_trew, sum_del_1_trew, sum_del_15_trew, sum_del_2_trew,
sum_del_0_mzeb, sum_del_05_mzeb, sum_del_1_mzeb, sum_del_15_mzeb, sum_del_2_mzeb,
sum_del_0_oto, sum_del_05_oto, sum_del_1_oto, sum_del_15_oto, sum_del_2_oto,
sum_del_0_rham, sum_del_05_rham, sum_del_1_rham, sum_del_15_rham, sum_del_2_rham)

count_sum_states_del <- colSums(sum_del) ## this will be for all delertions
count_sum_states_del <- as.data.frame(count_sum_states_del) %>% rownames_to_column("variable")
count_sum_dt_del <- as.data.table(count_sum_states_del)
count_sum_dt_del[, c("del", "del2","state", "pop_id") := tstrsplit(variable, "_", fixed=TRUE)]
count_sum_dt_del[,c("del", "del2") := NULL]
count_sum_dt_del[,state := as.factor(state)]
count_sum_dt_del[,pop_id := as.factor(pop_id)]

count_sum_dt_del[, state:=`levels<-`(state, c("0", "0.5", "1", "1.5", "2"))]
count_sum_dt_del[, pop_id:=`levels<-`(pop_id, c("Alticorpus peterdaviesi.Cape_Maclear", "Astatotilapia calliptera.Lake_Kingiri", 
"Astatotilapia calliptera.Lake_Masoko", "Copadichromis chrysonotus.Lake_Malombe", "Chilotilapia rhoadesii",
"Copadichromis virginalis.Lake_Malombe","Copadichromis virginalis.Southwest_arm", "Cynotilapia zebroides.Cape_Maclear",
"Diplotaxodon limnothrissa.Southwest_arm", "Fossorochromis rostratus.Lake_Malombe", "Labeotropheus fuelleborni.Chilumba", "Maylandia zebra.Cape_Maclear", "Otopharynx argyrosoma.Southeast_arm",
"Rhamphochromis longiceps", "Labeotropheus trewavasae.Chilumba"))]
levels.pop <- c("Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
"Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
"Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
"Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
"Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", "Diplotaxodon limnothrissa.Southwest_arm")
col <- "pop_id"
count_sum_dt_del[, (col) := factor(get(col), levels = levels.pop)]

count_sum_merged_del <- merge(count_sum_dt_del, meta_clade, by="pop_id")
count_sum_merged_del[, (col_clade) := factor(get(col_clade), levels = levels.clade)]
setnames(count_sum_merged_del, "pop_id", "population")
del_by_state <- ggplot(data = count_sum_merged_del, aes(x = population, y = count_sum_states_del, fill=state))+
  geom_col(aes(fill=state))+
  scale_fill_manual(values=c("#FFAB91", "#FF7043", "#F4511E", "#D84315", "#B71C1C"))+
  labs(title="Distribution of deletions in each population",x = "Population", y = "Number of deletions")+
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.margin = margin(10, 10, 10, 100))
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEA_nr_del_states_by_pop.pdf", width=22, height=10)
del_by_state
dev.off()

sum_del_0 <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_ids]
sum_del_05 <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_ids]
sum_del_1 <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_ids]
sum_del_15 <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_ids]
sum_del_2 <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_ids]
sum_ins_0 <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0"))), .SDcols=meta_ids]
sum_ins_05 <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="0.5"))), .SDcols=meta_ids]
sum_ins_1 <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))), .SDcols=meta_ids]
sum_ins_15 <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1.5"))), .SDcols=meta_ids]
sum_ins_2 <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="2"))), .SDcols=meta_ids]

sum_del_total <- sum(c(sum_del_0, sum_del_05, sum_del_1, sum_del_15, sum_del_2))
sum_ins_total <- sum(c(sum_ins_0, sum_ins_05, sum_ins_1, sum_ins_15, sum_ins_2))

sum_del <- as.data.table(cbind(sum_del_0, sum_del_05, sum_del_1, sum_del_15, sum_del_2))
sum_ins <- as.data.table(cbind(sum_ins_0, sum_ins_05, sum_ins_1, sum_ins_15, sum_ins_2))
sum_del_state <- as.data.frame(colSums(sum_del)) %>% rownames_to_column("variable")
sum_ins_state <- as.data.frame(colSums(sum_ins)) %>% rownames_to_column("variable")
sum_del_state$total <- sum_del_total
sum_ins_state$total <- sum_ins_total

sum_del_state <- as.data.table(sum_del_state)
sum_del_state <- sum_del_state[, "state" := tstrsplit(variable, "_", keep=3)]
setnames(sum_del_state, colnames(sum_del_state), c("variable", "sum_ins", "total", "state"))
sum_del_state <- sum_del_state[, percent := sum_ins / total]
sum_del_state[,state := as.factor(state)]
sum_del_state[, state:=`levels<-`(state, c("0", "0.5", "1", "1.5", "2"))]

sum_ins_state <- as.data.table(sum_ins_state)
sum_ins_state <- sum_ins_state[, "state" := tstrsplit(variable, "_", keep=3)]
setnames(sum_ins_state, colnames(sum_ins_state), c("variable", "sum_ins", "total", "state"))
sum_ins_state <- sum_ins_state[, percent := sum_ins / total]
sum_ins_state[,state := as.factor(state)]
sum_ins_state[, state:=`levels<-`(state, c("0", "0.5", "1", "1.5", "2"))]

pie_state_ins <- ggplot(sum_ins_state, aes(x="", y=percent, fill=state))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y")+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values=c("#FFAB91", "#FF7043", "#BF360C", "#D81B60", "#4A148C"))+
  guides(fill = guide_legend(title = "state")) +
  ggtitle("Sum of insertions (total 89,537,760)")

pie_state_del <- ggplot(sum_del_state, aes(x="", y=percent, fill=state))+
 geom_bar(width = 1, stat = "identity")+
 coord_polar("y")+
 theme_minimal()+
 theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
       axis.text.x=element_blank(),
       panel.border = element_blank(),
       panel.grid=element_blank(),
       axis.ticks = element_blank())+
 scale_fill_manual(values=c("#FFAB91", "#FF7043", "#BF360C", "#D81B60", "#4A148C"))+
 guides(fill = guide_legend(title = "state")) +
 ggtitle("Sum of insertions (shared with reference, total 3,600,720)")

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-11_states_pie.pdf", width=16, height=8)
grid.arrange(pie_state_ins, pie_state_del, ncol=2)
dev.off()

