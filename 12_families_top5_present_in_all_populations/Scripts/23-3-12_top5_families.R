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

cbPalette <- c("#EF5350", "#2196F3", "#43A047", "#9575CD", "#BDBDBD", "#D1C4E9", "#757575", "#795548") # orders
popPalette <- c("#9CCC65", "#7CB342", "#7E57C2", "512DA8", "#AB47BC", "#6A1B9A", "#1B5E20", "8BC34A", "388E3C",
                "#B71C1C", "#E57373", "#D32F2F", "#1A237E", "#BF360C", "#827717")
cladePalette <- c("#AED581", "#5E35B1", "#1B5E20", "#B71C1C", "#1A237E", "#EF6C00", "#827717")

meta <- fread("~/Desktop/23-2-16_cichlid_meta_standardized_dataset.csv",
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

dt5 <- fread("~/Desktop/23-3-11_top5_fam_by_pop.txt", sep="\t")
dt10 <- fread("~/Desktop/23-3-11_top10_fam_by_pop.txt", sep="\t")
dt5_filt <- dt5[,c("pop_id", "family", "myrank"), with=FALSE]
dt5_wide <- spread(dt5_filt, pop_id, family)
dt5_wide[, myrank:= as.factor(myrank)]
dt5_wide[, myrank:=`levels<-`(myrank, c("5", "4", "3", "2", "1"))]
setcolorder(dt5_wide, c("myrank", "asm", "ask", "fuel","trew", "czeb", "mzeb", "cvs", "cvm", "cch",
                        "crhoa", "fr", "oto", "ap", "rham", "diplo"))
setnames(dt5_wide, colnames(dt5_wide), c("rank", "Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
                                         "Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
                                         "Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
                                         "Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
                                         "Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", 
                                         "Diplotaxodon limnothrissa.Southwest_arm"))
dt5_wide <- dt5_wide[order(-rank)]
write.table(dt5_wide, "~/Desktop/23-3-12_top5_ranks.txt", sep="\t", quote=T, row.names = F)

