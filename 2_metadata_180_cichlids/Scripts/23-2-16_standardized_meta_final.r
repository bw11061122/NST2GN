## 23-2-16
#####################################################
#####################################################
## Finalised metadata (version 3)

### This script is to make sure that the clades are alright
library(circlize)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(EnvStats)
library(dplyr)
library(tibble)
library(ggrepel)
library("FactoMineR") ## PCA
library("factoextra") ## PCA 
library(stringr)

#########################################
#########################################
# select species with at least 12 homogenous individuals
# species that I want to look at
# deep - Alticorpus peterdaviesi (1)
# Diplo - Diplotaxodon limnothrissa (1)
# Rhampho - Rhamphochromic longiceps (1)
# AstCal - Astatotilapia calliptera (Lake Masoko, Kingiri lake) (2)
# Mbuna - M. zebra, L. trewavasae, L. fuelleborni, C. zebroides (Thumbi west) (4)
# Utaka - C. virginalis (Southwest Arm, Lake Malombe), C. chrysonotus (Lake Malombe) (3)
# Benthic - Fossorochromis rostratus (Southwest Arm), Otopharynx argyrosoma, C. rhoadesii (3)

#########################################
### cichlid metadata
cichlid_all_metadata <- fread("~/Desktop/Genetics_project/23-1-23/Data/meta/cichlid_callset_metadata.csv",
                              sep=",", header = T) ## this is the full metadat
cichlid_all_metadata <- cichlid_all_metadata[, clade:=as.factor(clade)] 
cichlid_all_metadata <- cichlid_all_metadata[, sex:=as.factor(sex)] 
cichlid_all_metadata <- cichlid_all_metadata[, location:=as.factor(location)] 
cichlid_all_metadata <- cichlid_all_metadata[, sublocation:=as.factor(sublocation)] 
cichlid_all_metadata <- cichlid_all_metadata[, species:=as.factor(species)] 
cichlid_all_metadata <- cichlid_all_metadata[, genus:=as.factor(genus)] 
cichlid_all_metadata <- cichlid_all_metadata[, names:=paste(genus, species)] 
cichlid_all_metadata <- cichlid_all_metadata[, name_loc:=paste0(names,".",location)] 
cichlid_all_metadata <- cichlid_all_metadata[, names:=as.factor(names)] 
cichlid_all_metadata <- cichlid_all_metadata[, name_loc:=as.factor(name_loc)] 

######## choice of sequencing depth 
## plot histogram of sequence depth 
hist(cichlid_all_metadata[,seq_depth], breaks=100)
### This shows that the threshold at 10 is very reasonable 
## threshold = 10 (see the histogram)

###### choice of species 
## how many cichlids do I have for each 
cichlid_number_samples_name <- cichlid_all_metadata[seq_depth > 10, .N, keyby = names][order (-N)]
## how many of these have species name without sp
## I don't want to have the ones with .sp bc these are ones where we are not sure what the species is
cichlid_number_samples_name_12ind <- cichlid_number_samples_name[N >= 12,]
cichlid_number_samples_name_12ind_filtered <-cichlid_number_samples_name_12ind[- grep("sp. ", names),]

###### 23-2-16 meeting with Richard and Pio - decided to focus on 15 cichlid populations
## filter by sequence
cichlid_all_metadata_seq <- cichlid_all_metadata[seq_depth > 10,]

### cichlids of interest 
populations_of_interest <- c("Alticorpus peterdaviesi.Cape_Maclear", "Diplotaxodon limnothrissa.Southwest_arm",
                             "Rhamphochromis longiceps", "Astatotilapia calliptera.Lake_Masoko", 
                             "Astatotilapia calliptera.Lake_Kingiri", "Maylandia zebra.Cape_Maclear",
                             "Labeotropheus fuelleborni.Chilumba", "Labeotropheus trewavasae.Chilumba",
                             "Cynotilapia zebroides.Cape_Maclear", "Copadichromis virginalis.Lake_Malombe",
                             "Copadichromis virginalis.Southwest_arm", "Copadichromis chrysonotus.Lake_Malombe",
                             "Fossorochromis rostratus.Lake_Malombe", "Chilotilapia rhoadesii",
                             "Otopharynx argyrosoma.Southeast_arm")
meta_pop_of_interest <- cichlid_all_metadata_seq[name_loc %in% c(populations_of_interest) | names %in% c(populations_of_interest),]

## subset for populations you are interested in 
meta_filtered_1 <- meta_pop_of_interest[name_loc=="Alticorpus peterdaviesi.Cape_Maclear"]
meta_filtered_2 <- meta_pop_of_interest[name_loc=="Astatotilapia calliptera.Lake_Masoko"]
meta_filtered_3 <- meta_pop_of_interest[name_loc=="Astatotilapia calliptera.Lake_Kingiri"]
meta_filtered_4 <- meta_pop_of_interest[name_loc=="Copadichromis chrysonotus.Lake_Malombe"]
meta_filtered_5 <- meta_pop_of_interest[name_loc=="Copadichromis virginalis.Southwest_arm"]
meta_filtered_6 <- meta_pop_of_interest[name_loc=="Copadichromis virginalis.Lake_Malombe"]
meta_filtered_7 <- meta_pop_of_interest[name_loc=="Cynotilapia zebroides.Cape_Maclear"]
meta_filtered_8 <- meta_pop_of_interest[name_loc=="Diplotaxodon limnothrissa.Southwest_arm"]
meta_filtered_9 <- meta_pop_of_interest[name_loc=="Fossorochromis rostratus.Lake_Malombe"]
meta_filtered_10 <- meta_pop_of_interest[name_loc=="Labeotropheus fuelleborni.Chilumba"]
meta_filtered_11 <- meta_pop_of_interest[name_loc=="Labeotropheus trewavasae.Chilumba"]
meta_filtered_12 <- meta_pop_of_interest[name_loc=="Maylandia zebra.Cape_Maclear"]
meta_filtered_13 <- meta_pop_of_interest[name_loc=="Otopharynx argyrosoma.Southeast_arm"]
meta_filtered_14 <- meta_pop_of_interest[names=="Chilotilapia rhoadesii"]
meta_filtered_15 <- meta_pop_of_interest[names=="Rhamphochromis longiceps"]

# choose 12 homogenous individuals from each 
### for the Masoko samples, I don't know - do I care about sublocation too?
### if no sublocation, sample randomly
set.seed(1)
meta_filtered_1_final <- meta_filtered_1[sample(.N, 12)]

## for AstCal, I want littoral in the notes
## there are no cichlids marked littoral
meta_filtered_2_sub <- meta_filtered_2[sublocation=="5-20m",]
meta_filtered_2_final <- meta_filtered_2_sub[sample(.N, 12)]

## choose by sublocation
meta_filtered_3_sub <- meta_filtered_3[sublocation=="shallow",]
meta_filtered_3_final <- meta_filtered_3_sub[sample(.N, 12)]

# choose by sublocation
meta_filtered_4_sub <- meta_filtered_4[sublocation=="Chimwala",]
meta_filtered_4_final <- meta_filtered_4_sub[sample(.N, 12)]

# make sure clade is okay
meta_filtered_5_sub <- meta_filtered_5[clade=="Utaka" & sublocation=="Msaka",]
meta_filtered_5_final <- meta_filtered_5_sub[sample(.N, 12)]

# idk exclude notes
meta_filtered_6_sub <- meta_filtered_6[notes == "",]
meta_filtered_6_final <- meta_filtered_6_sub[sample(.N, 12)]

# choose by sublocation
meta_filtered_7_sub <- meta_filtered_7[sublocation == "Thumbi_West",]
meta_filtered_7_final <- meta_filtered_7_sub[sample(.N, 12)]

### if no sublocation, sample randomly
meta_filtered_8_final <- meta_filtered_8[sample(.N, 12)]
meta_filtered_9_final <- meta_filtered_9[sample(.N, 12)]

meta_filtered_10_sub <- meta_filtered_10[sublocation == "Luwino_reef",]
meta_filtered_10_final <- meta_filtered_10_sub[sample(.N, 12)]

meta_filtered_11_sub <- meta_filtered_11[sublocation == "Luwino_reef",]
meta_filtered_11_final <- meta_filtered_11[sample(.N, 12)]

meta_filtered_12_sub <- meta_filtered_12[sublocation == "Thumbi_West",]
meta_filtered_12_final <- meta_filtered_12_sub[sample(.N, 12)]

meta_filtered_13_final <- meta_filtered_13[sample(.N, 12)]

meta_filtered_14_sub <- meta_filtered_14[sublocation != "Makanjila",]
meta_filtered_14_final <- meta_filtered_14_sub[sample(.N, 12)]

meta_filtered_15_final <- meta_filtered_15[sample(.N, 12)]

meta_all15_final <- rbind(meta_filtered_1_final, meta_filtered_2_final, meta_filtered_3_final,
                          meta_filtered_4_final, meta_filtered_5_final, meta_filtered_6_final,
                          meta_filtered_7_final, meta_filtered_8_final, meta_filtered_9_final,
                          meta_filtered_10_final, meta_filtered_11_final, meta_filtered_12_final,
                          meta_filtered_13_final, meta_filtered_14_final, meta_filtered_15_final)
hist(meta_all15_final[,seq_depth], breaks=100)
summary(meta_all15_final)
#### NOTE: 23-2-14 I am saving a suggested callset (not consulted completely with Richard and Pio)
#### just to get some more stuff done but can over-write with new file
write.csv(meta_all15_final, "~/Desktop/Genetics_project/23-2-16/23-2-16_cichlid_meta_standardized_dataset.csv", 
          row.names=F, quote=F)

#########################
##########################
## NOTE: for clade 16
### there are 2 out of 12 fish which are hard to ID 
### I am not removing them but need to keep this in mind

dt <- data.table(col1 = c(1:9), col2=3, col3 = "blue")
dtx <- data.table(col1 = c(1:9), col2=3)
write.table(dt, "~/Desktop/mock.txt", row.names=F, col.names=T)
dt2 <- fread("~/Desktop/mock.txt", header=T)
dt_sums1 <- colSums(dt)  
dt_sum1 <- sapply(dtx, function(x) sum(x))
dt_sum2 <- sapply(dt, function(x) sum(x != 1))
dt_sum3 <- sapply(dt, function(x) sum(x == 1))
dt_sum3 <- sapply(dt, function(x) sum(x %in% c(1:5)))
dt_sum4 <- sapply(round(dt,1), function(x) sum(x!=1))

dt_sums2 <- sapply(round(dt[-1],1), function(x) sum(!(x %in% c(2, 4))))
dt_sums3 <- sapply(dt, function(x) sum(!(x %in% c(2,4))))
dt_sums4 <- colSums(dt[,1:2] != 4) 

library(dplyr)
round(dt, 1) %>% summarise(across(-1, ~sum(!(. %in% c(2,4) | is.na(.)))))