#### [2023-03-09]
#### In this script, I am analysing the general metadata set 
#### Questions: how many clades, species there are? how many locations have been sampled?

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
"/rds/project/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara"

### read the metadata
meta <- fread("data/cichlid_callset_metadata.txt",
                              sep="\t", header = T) ## this is the full metadata
meta <- meta[, clade:=as.factor(clade)] 
meta <- meta[, sex:=as.factor(sex)] 
meta <- meta[, location:=as.factor(location)] 
meta <- meta[, sublocation:=as.factor(sublocation)] 
meta <- meta[, species:=as.factor(species)] 
meta <- meta[, genus:=as.factor(genus)] 
meta <- meta[,names := paste(genus, species, sep="_")]
meta <- meta[, names:=as.factor(names)] 
meta <- meta[, name_loc:=as.factor(name_loc)] 
upd.cols = sapply(meta, is.factor)
meta <- meta[, names(meta)[upd.cols] := lapply(.SD, factor), .SDcols = upd.cols] ## this gets rid of levels which are absent in the dt 
meta_ids <- unlist(meta[,primary_id]) 
length(meta_ids) 

### how many clades are there?
levels(meta[,clade])
### how many species are there?
levels(meta[,names]) ## I can see that in many cases, the species is not clear
species_sp <- meta %>% as.data.frame %>% select(names) %>% filter(grepl('_sp. ', names))  
species_cf <- meta %>% as.data.frame %>% select(names) %>% filter(grepl('_cf. ', names)) 
species_sp$names <- as.factor(species_sp$names) 
species_cf$names <- as.factor(species_cf$names) 
length(unique(species_sp$names)) # 91
length(unique(species_cf$names)) # 24

### from how many locations were the samples taken?
levels(meta[,location])

