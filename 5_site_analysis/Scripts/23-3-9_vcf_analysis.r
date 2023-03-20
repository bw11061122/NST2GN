#### [2023-03-09]
#### In this script, I am analysing the general metadata set 
#### Questions: how many insertions are contributed by each order / superfamily / family (in total / per population)?

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

### set up coloring 
cbPalette <- c("#EF5350", "#2196F3", "#43A047", "#9575CD", "#BDBDBD", "#D1C4E9", "#757575", "#795548") # orders
popPalette <- c("#9CCC65", "#7CB342", "#7E57C2", "512DA8", "#AB47BC", "#6A1B9A", "#1B5E20", "8BC34A", "388E3C",
"#B71C1C", "#E57373", "#D32F2F", "#1A237E", "#BF360C", "#827717")
cladePalette <- c("#AED581", "#5E35B1", "#1B5E20", "#B71C1C", "#1A237E", "#EF6C00", "#827717")

### Read the metadata
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

### set order of populations
## I think makes more sense to have it from A. calliptera - Mbuna - Utaka - Benthic - Deep - Rhampo - Diplo
levels.pop <- c("Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
"Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
"Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
"Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
"Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", "Diplotaxodon limnothrissa.Southwest_arm")
col <- "pop_id"
meta[, (col) := factor(get(col), levels = levels.pop)]

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

### Read vcf
vcf <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_ins_del_all_final.txt", 
    sep = "\t", header=T, fill=TRUE)
 
vcf <- vcf[, order:=as.factor(order)]
vcf <- vcf[, family:=as.factor(family)]
vcf <- vcf[, superfamily:=as.factor(superfamily)]
upd.cols.vcf = sapply(vcf, is.factor)
vcf <- vcf[, names(vcf)[upd.cols.vcf] := lapply(.SD, factor), .SDcols = upd.cols.vcf] ## this gets rid of levels which are absent in the dt 

### Contribution of orders, superfamilies, families
vcf_sum_order <- vcf[,.N, by=order]
vcf_sum_superfamily <- vcf[,.N, by=superfamily]
vcf_sum_family <- vcf[,.N, by=family]

vcf_sum_order <- vcf_sum_order[,percent:= N / sum(N) * 100][order(percent)]
vcf_sum_superfamily <- vcf_sum_superfamily[,percent:= N / sum(N) * 100][order(percent)]
vcf_sum_family <- vcf_sum_family[,percent:= N / sum(N) * 100][order(percent)]

## relevel orders
levels.order <- c("DNA", "LINE", "LTR", "SINE", "Unknown", "RC", "Retrotransposon", "rRNA")
col <- "order"
vcf_sum_order[, (col) := factor(get(col), levels = levels.order)]

vcf_sum_order[is.na(order==TRUE), order := "Unknown"]
upd.cols = sapply(vcf_sum_order, is.factor)
box_TEs_order <- ggplot(data = vcf_sum_order, aes(x = order, y = percent, fill=order))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=cbPalette)+
  theme_bw()+
  labs(title="TE orders in 180 samples",x = "Order", y = "Percentage of the total number of sites (%)")+
  theme_bw(base_size = 26)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_nr_ins_by_order.pdf", width=14, height=10)
box_TEs_order
dev.off()

### let's do this for each population - now NUMBER
vcf_sum_order_pop <- vcf[,.N, by=.(order)]
vcf_sum_superfamily_pop <- vcf[,.N, by=.(superfamily)]
vcf_sum_family_pop <- vcf[,.N, by=.(family)]

vcf_sum_order <- vcf_sum_order[,percent:= N / sum(N) * 100][order(percent)]
vcf_sum_superfamily <- vcf_sum_superfamily[,percent:= N / sum(N) * 100][order(percent)]
vcf_sum_family <- vcf_sum_family[,percent:= N / sum(N) * 100][order(percent)]

## get vcf for each family
info_col <- colnames(vcf[,c(1:14, 195:198)])
vcf_ap <- vcf[,c(info_col, meta_ap), with=FALSE]
vcf_ask <- vcf[,c(info_col, meta_ask), with=FALSE]
vcf_asm <- vcf[,c(info_col, meta_asm), with=FALSE]
vcf_cch <- vcf[,c(info_col, meta_cch), with=FALSE]
vcf_crhoa <- vcf[,c(info_col, meta_crhoa), with=FALSE]
vcf_cvm <- vcf[,c(info_col, meta_cvm), with=FALSE]
vcf_cvs <- vcf[,c(info_col, meta_cvs), with=FALSE]
vcf_czeb <- vcf[,c(info_col, meta_czeb), with=FALSE]
vcf_diplo <- vcf[,c(info_col, meta_diplo), with=FALSE]
vcf_fr <- vcf[,c(info_col, meta_fr), with=FALSE]
vcf_fuel <- vcf[,c(info_col, meta_fuel), with=FALSE]
vcf_trew <- vcf[,c(info_col, meta_trew), with=FALSE]
vcf_mzeb <- vcf[,c(info_col, meta_mzeb), with=FALSE]
vcf_oto <- vcf[,c(info_col, meta_oto), with=FALSE]
vcf_rham <- vcf[,c(info_col, meta_rham), with=FALSE]

## get ones which are above 0
vcf_ap_nz <- vcf_ap[ rowSums(vcf_ap[,19:30]) > 0, ]
vcf_ask_nz <- vcf_ap[ rowSums(vcf_ask[,19:30]) > 0, ]
vcf_asm_nz <- vcf_ap[ rowSums(vcf_asm[,19:30]) > 0, ]
vcf_cch_nz <- vcf_ap[ rowSums(vcf_cch[,19:30]) > 0, ]
vcf_crhoa_nz <- vcf_ap[ rowSums(vcf_crhoa[,19:30]) > 0, ]
vcf_cvm_nz <- vcf_ap[ rowSums(vcf_cvm[,19:30]) > 0, ]
vcf_cvs_nz <- vcf_ap[ rowSums(vcf_cvs[,19:30]) > 0, ]
vcf_czeb_nz <- vcf_ap[ rowSums(vcf_czeb[,19:30]) > 0, ]
vcf_diplo_nz <- vcf_ap[ rowSums(vcf_diplo[,19:30]) > 0, ]
vcf_fr_nz <- vcf_ap[ rowSums(vcf_fr[,19:30]) > 0, ]
vcf_fuel_nz <- vcf_ap[ rowSums(vcf_fuel[,19:30]) > 0, ]
vcf_trew_nz <- vcf_ap[ rowSums(vcf_trew[,19:30]) > 0, ]
vcf_mzeb_nz <- vcf_ap[ rowSums(vcf_mzeb[,19:30]) > 0, ]
vcf_oto_nz <- vcf_ap[ rowSums(vcf_oto[,19:30]) > 0, ]
vcf_rham_nz <- vcf_ap[ rowSums(vcf_rham[,19:30]) > 0, ]

num_ap_order <- vcf_ap_nz[, .N, by=.(order)]
total_ap <- sum(num_ap_order[,N])
num_ap_order <- num_ap_order[,percent := N / total_ap * 100]
num_ap_order <- num_ap_order[,population := "Alticorpus peterdaviesi.Cape_Maclear"]
num_ask_order <- vcf_ask_nz[, .N, by=.(order)]
total_ask <- sum(num_ask_order[,N])
num_ask_order <- num_ask_order[,percent := N / total_ask * 100]
num_ask_order <- num_ask_order[,population := "Astatotilapia calliptera.Lake_Kingiri"]
num_asm_order <- vcf_asm_nz[, .N, by=.(order)]
total_asm <- sum(num_asm_order[,N])
num_asm_order <- num_asm_order[,percent := N / total_asm * 100]
num_asm_order <- num_asm_order[,population := "Astatotilapia calliptera.Lake_Masoko"]
num_cch_order <- vcf_cch_nz[, .N, by=.(order)]
total_cch <- sum(num_cch_order[,N])
num_cch_order <- num_cch_order[,percent := N / total_cch * 100]
num_cch_order <- num_cch_order[,population := "Copadichromis chrysonotus.Lake_Malombe"]
num_crhoa_order <- vcf_crhoa_nz[, .N, by=.(order)]
total_crhoa <- sum(num_crhoa_order[,N])
num_crhoa_order <- num_crhoa_order[,percent := N / total_crhoa * 100]
num_crhoa_order <- num_crhoa_order[,population := "Chilotilapia rhoadesii"]
num_cvm_order <- vcf_cvm_nz[, .N, by=.(order)]
total_cvm <- sum(num_cvm_order[,N])
num_cvm_order <- num_cvm_order[,percent := N / total_cvm * 100]
num_cvm_order <- num_cvm_order[,population := "Copadichromis virginalis.Lake_Malombe"]
num_cvs_order <- vcf_cvs_nz[, .N, by=.(order)]
total_cvs <- sum(num_cvs_order[,N])
num_cvs_order <- num_cvs_order[,percent := N / total_cvs * 100]
num_cvs_order <- num_cvs_order[,population := "Copadichromis virginalis.Southwest_arm"]
num_czeb_order <- vcf_czeb_nz[, .N, by=.(order)]
total_czeb <- sum(num_czeb_order[,N])
num_czeb_order <- num_czeb_order[,percent := N / total_czeb * 100]
num_czeb_order <- num_czeb_order[,population := "Cynotilapia zebroides.Cape_Maclear"]
num_diplo_order <- vcf_diplo_nz[, .N, by=.(order)]
total_diplo <- sum(num_diplo_order[,N])
num_diplo_order <- num_diplo_order[,percent := N / total_diplo * 100]
num_diplo_order <- num_diplo_order[,population := "Diplotaxodon limnothrissa.Southwest_arm"]
num_fr_order <- vcf_fr_nz[, .N, by=.(order)]
total_fr <- sum(num_fr_order[,N])
num_fr_order <- num_fr_order[,percent := N / total_fr * 100]
num_fr_order <- num_fr_order[,population := "Fossorochromis rostratus.Lake_Malombe"]
num_fuel_order <- vcf_fuel_nz[, .N, by=.(order)]
total_fuel <- sum(num_fuel_order[,N])
num_fuel_order <- num_fuel_order[,percent := N / total_fuel * 100]
num_fuel_order <- num_fuel_order[,population := "Labeotropheus fuelleborni.Chilumba"]
num_trew_order <- vcf_trew_nz[, .N, by=.(order)]
total_trew <- sum(num_trew_order[,N])
num_trew_order <- num_trew_order[,percent := N / total_trew * 100]
num_trew_order <- num_trew_order[,population := "Labeotropheus trewavasae.Chilumba"]
num_mzeb_order <- vcf_mzeb_nz[, .N, by=.(order)]
total_mzeb <- sum(num_mzeb_order[,N])
num_mzeb_order <- num_mzeb_order[,percent := N / total_mzeb * 100]
num_mzeb_order <- num_mzeb_order[,population := "Maylandia zebra.Cape_Maclear"]
num_oto_order <- vcf_oto_nz[, .N, by=.(order)]
total_oto <- sum(num_oto_order[,N])
num_oto_order <- num_oto_order[,percent := N / total_oto * 100]
num_oto_order <- num_oto_order[,population := "Otopharynx argyrosoma.Southeast_arm"]
num_rham_order <- vcf_rham_nz[, .N, by=.(order)]
total_rham <- sum(num_rham_order[,N])
num_rham_order <- num_rham_order[,percent := N / total_rham * 100]
num_rham_order <- num_rham_order[,population := "Rhamphochromis longiceps"]

nums_order_pop <- rbind(num_ap_order, 
num_ask_order, num_asm_order, num_cch_order,
num_crhoa_order, num_cvm_order, num_cvs_order,
num_czeb_order, num_diplo_order, num_fr_order,
num_fuel_order, num_trew_order, num_mzeb_order,
num_oto_order, num_rham_order)

levels.order <- c("DNA", "LINE", "LTR", "SINE", "Unknown", "RC", "Retrotransposon", "rRNA")
col_order <- "order"
nums_order_pop[, (col_order) := factor(get(col_order), levels = levels.order)]
nums_order_pop[is.na(order==TRUE), order := "Unknown"]
upd.cols = sapply(nums_order_pop, is.factor)
levels.pop <- c("Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
"Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
"Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
"Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
"Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", "Diplotaxodon limnothrissa.Southwest_arm")
col_meta <- "population"
nums_order_pop[, (col_meta) := factor(get(col_meta), levels = levels.pop)]
nums_order_pop <- nums_order_pop[, names(nums_order_pop)[upd.cols] := lapply(.SD, factor), .SDcols = upd.cols]
bar_TEs_order_pop <- ggplot(data = nums_order_pop, aes(x = population, y = percent, fill=order))+
  geom_col(aes(fill=order))+
  scale_fill_manual(values=cbPalette)+
  theme_bw()+
  labs(title="TE orders in 180 samples",x = "Population", y = "Percentage of the total number of sites (%)")+
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.margin = margin(10, 10, 10, 100))
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_MEI_nr_ins_by_order_by_pop.pdf", width=22, height=10)
bar_TEs_order_pop
dev.off()

total_nr_sites <- as.data.table(rbind(total_ap, total_ask, total_asm, total_cch, total_crhoa, total_cvm, total_cvs,
total_czeb, total_diplo, total_fr, total_fuel, total_trew, total_mzeb, total_oto, total_rham))
total_nr_sites <- total_nr_sites[,pop_id := c("Alticorpus peterdaviesi.Cape_Maclear",
"Astatotilapia calliptera.Lake_Kingiri", "Astatotilapia calliptera.Lake_Masoko", "Chilotilapia rhoadesii",
"Copadichromis chrysonotus.Lake_Malombe", "Copadichromis virginalis.Lake_Malombe", "Copadichromis virginalis.Southwest_arm", 
"Cynotilapia zebroides.Cape_Maclear", "Diplotaxodon limnothrissa.Southwest_arm", "Fossorochromis rostratus.Lake_Malombe", 
"Labeotropheus fuelleborni.Chilumba", "Labeotropheus trewavasae.Chilumba", "Maylandia zebra.Cape_Maclear", 
"Otopharynx argyrosoma.Southeast_arm", "Rhamphochromis longiceps")]
total_nr_sites <- total_nr_sites[,pop_id := as.factor(pop_id)]
levels.pop <- c("Astatotilapia calliptera.Lake_Masoko", "Astatotilapia calliptera.Lake_Kingiri", "Labeotropheus fuelleborni.Chilumba", 
"Labeotropheus trewavasae.Chilumba", "Cynotilapia zebroides.Cape_Maclear", "Maylandia zebra.Cape_Maclear", 
"Copadichromis virginalis.Southwest_arm", "Copadichromis virginalis.Lake_Malombe",  "Copadichromis chrysonotus.Lake_Malombe",
"Chilotilapia rhoadesii", "Fossorochromis rostratus.Lake_Malombe","Otopharynx argyrosoma.Southeast_arm", 
"Alticorpus peterdaviesi.Cape_Maclear", "Rhamphochromis longiceps", "Diplotaxodon limnothrissa.Southwest_arm")
col_meta <- "pop_id"
total_nr_sites[, (col_meta) := factor(get(col_meta), levels = levels.pop)]
total_nr_sites <- total_nr_sites[, names(total_nr_sites)[upd.cols] := lapply(.SD, factor), .SDcols = upd.cols]
meta_clade <- meta[,c("pop_id", "clade"), with=FALSE]
total_nr_sites <- merge(meta_clade, total_nr_sites, by="pop_id")
total_nr_sites <- unique(total_nr_sites)
levels.clade <- c("AstCal", "Mbuna", "Utaka", "Benthic", "Deep","Rhampho", "Diplo")
col_clade <- "clade"
total_nr_sites[, (col_clade) := factor(get(col_clade), levels = levels.clade)]
bar_nr_sites <- ggplot(data = total_nr_sites, aes(x = pop_id, y = V1, fill=clade))+
  geom_bar(stat="identity")+
  ylim(c(0,140000))+
  scale_fill_manual(values=cladePalette)+
  theme_bw()+
  labs(title="Distribution of TE sites",x = "Population", y = "Total number of sites")+
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.margin = margin(10, 10, 10, 100))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-9_nr_sites_by_pop.pdf", width=12, height=10)
bar_nr_sites
dev.off()
