#### [2023-2-27]
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

#### QUESTION 1
## How many polymorphic insertions are there in each sample?  by populations and clades
### Firstly, I can determine how many insertions are carried by each sample 
sum_pol_per_sample <- colSums(vcf[,..meta_ids]) ### pol - polymorphisms 
meta <- cbind(meta, sum_pol_per_sample)

box_TEs_per_sample <- ggplot(data = meta, aes(x = pop_id, y = sum_pol_per_sample, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  scale_color_manual(values=cladePalette)+
  ylim(c(0,65000))+
  labs(title="Sum of polymorphisms per sample",x = "Population", y = "Sum of polymorphisms per sample")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_polymorphisms_per_sample.pdf", width=18, height=14)
box_TEs_per_sample
dev.off()

#######################################################
#### insertions (MEI - not present in the reference)
vcf_ins <- vcf[type=="MEI",]
sum_ins_per_sample <- colSums(vcf_ins[,..meta_ids])
meta <- cbind(meta, sum_ins_per_sample)
box_TEs_per_sample_ins <- ggplot(data = meta, aes(x = pop_id, y = sum_ins_per_sample, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  scale_color_manual(values=cladePalette)+
  ylim(c(0,65000))+
  labs(title="Sum of polymorphic insertions per sample",x = "Population", y = "Sum of insertions per sample")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-2_MEI_sum_ins_per_sample.pdf", width=18, height=14)
box_TEs_per_sample_ins
dev.off()

#### deletions (MEA - present in the reference)
vcf_del <- vcf[type=="MEA",]
sum_del_per_sample <- colSums(vcf_del[,..meta_ids])
meta <- cbind(meta, sum_del_per_sample)
box_TEs_per_sample_del <- ggplot(data = meta, aes(x = pop_id, y = sum_del_per_sample, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  ylim(c(0,65000))+
  scale_color_manual(values=cladePalette)+
  labs(title="Sum of polymorphic deletions per sample",x = "Population", y = "Sum of deletions per sample")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-2_MEI_sum_del_per_sample.pdf", width=18, height=14)
box_TEs_per_sample_del
dev.off()

#######################################################
#######################################################
## analyse heterozygous sites not the sum
vcf_num_all <- vcf[,meta_ids, with=FALSE]
sum_het_sites_all <- apply(vcf_num_all, 2, function(x) length(which(x=="1")))
meta <- cbind(meta, sum_het_sites_all)
box_TEs_per_sample_het_all <- ggplot(data = meta, aes(x = pop_id, y = sum_het_sites_all, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  scale_color_manual(values=cladePalette)+
  ylim(c(0,20000))+
  theme_bw()+
  labs(title="Number of heterozygous polymorphisms per sample",x = "Population", y = "Heterozygous sites per sampl")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_het_sites_per_sample_all.pdf", width=18, height=14)
box_TEs_per_sample_het_all
dev.off()

vcf_num_ins <- vcf_ins[,meta_ids, with=FALSE]
sum_het_sites_ins <- apply(vcf_num_ins, 2, function(x) length(which(x=="1")))
meta <- cbind(meta, sum_het_sites_ins)
box_TEs_per_sample_het_ins <- ggplot(data = meta, aes(x = pop_id, y = sum_het_sites_ins, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  ylim(c(0,20000))+
  scale_color_manual(values=cladePalette)+
  theme_bw()+
  labs(title="Number of heterozygous insertions per sample",x = "Population", y = "Heterozygous insertions per sample")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_het_sites_per_sample_ins.pdf", width=18, height=12)
box_TEs_per_sample_het_ins
dev.off()

vcf_num_del <- vcf_del[,meta_ids, with=FALSE]
sum_het_sites_del <- apply(vcf_num_del, 2, function(x) length(which(x=="1")))
meta <- cbind(meta, sum_het_sites_del)
box_TEs_per_sample_het_del <- ggplot(data = meta, aes(x = pop_id, y = sum_het_sites_del, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  ylim(c(0,20000))+
  scale_color_manual(values=cladePalette)+
  labs(title="Number of heterozygous deletions per sample",x = "Population", y = "Heterzoygous deletions per sample")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_het_sites_per_sample_del.pdf", width=18, height=12)
box_TEs_per_sample_het_del
dev.off()

###### add the number of homozygous sites
###### you can also caclculate the ratio 
sum_hom_sites_all <- apply(vcf_num_all, 2, function(x) length(which(x=="2"))) ### homozygotes which HAVE the insertion
het_to_hom_all <- sum_het_sites_all / sum_hom_sites_all
meta <- cbind(meta, het_to_hom_all)
box_TEs_per_sample_het_hom_all <- ggplot(data = meta, aes(x = pop_id, y = het_to_hom_all, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  ylim(c(0,2))+
  theme_bw()+
  scale_color_manual(values=cladePalette)+
  labs(title="Ratio of heterozygous to homozygous polymorphisms per sample",x = "Population", y = "Ratio of heterozygus to homozygous sites")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_het_to_hom_all.pdf", width=18, height=14)
box_TEs_per_sample_het_hom_all
dev.off()

sum_hom_sites_ins <- apply(vcf_num_ins, 2, function(x) length(which(x=="2"))) ### homozygotes which HAVE the insertion
het_to_hom_ins <- sum_het_sites_ins / sum_hom_sites_ins
meta <- cbind(meta, het_to_hom_ins)
box_TEs_per_sample_het_hom_ins <- ggplot(data = meta, aes(x = pop_id, y = het_to_hom_ins, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  ylim(c(0,8))+
  scale_color_manual(values=cladePalette)+
  theme_bw()+
  labs(title="Ratio of heterozygous to homozygous insertions per sample",x = "Population", y = "Ratio of heterozygus to homozygous sites")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_het_to_hom_ins.pdf", width=18, height=14)
box_TEs_per_sample_het_hom_ins
dev.off()

sum_hom_sites_del <- apply(vcf_num_del, 2, function(x) length(which(x=="2"))) ### homozygotes which HAVE the insertion
het_to_hom_del <- sum_het_sites_del / sum_hom_sites_del
meta <- cbind(meta, het_to_hom_del)
box_TEs_per_sample_het_hom_del <- ggplot(data = meta, aes(x = pop_id, y = het_to_hom_del, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  ylim(c(0,2))+
  scale_color_manual(values=cladePalette)+
  theme_bw()+
  labs(title="Ratio of heterozygous to homozygous deletions per sample",x = "Population", y = "Ratio of heterozygus to homozygous sites")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_het_to_hom_del.pdf", width=18, height=14)
box_TEs_per_sample_het_hom_del
dev.off()

#######################################################
#######################################################
## analyse heterozygous sites not the sum
sum_hethom_sites_all <- apply(vcf_num_all, 2, function(x) length(which(x=="1")) + 2*length(which(x=="2")))
meta <- cbind(meta, sum_hethom_sites_all)
box_TEs_per_sample_hethom_all <- ggplot(data = meta, aes(x = pop_id, y = sum_hethom_sites_all, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  scale_color_manual(values=cladePalette)+
  ylim(c(0,60000))+
  theme_bw()+
  labs(title="Total number of polymorphic insertions",x = "Population", y = "Sum of insertions (per sample)")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_hethom_sites_per_sample_all.pdf", width=18, height=14)
box_TEs_per_sample_hethom_all
dev.off()

sum_hethom_sites_ins <- apply(vcf_num_ins, 2, function(x) length(which(x=="1")) + 2*length(which(x=="2")))
meta <- cbind(meta, sum_hethom_sites_ins)
box_TEs_per_sample_hethom_ins <- ggplot(data = meta, aes(x = pop_id, y = sum_hethom_sites_ins, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  ylim(c(0,60000))+
  scale_color_manual(values=cladePalette)+
  theme_bw()+
  labs(title="Total number of polymorphic insertions (absent in reference)",x = "Population", y = "Sum of insertions (per sample)")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_hethom_sites_per_sample_ins.pdf", width=18, height=12)
box_TEs_per_sample_hethom_ins
dev.off()

sum_hethom_sites_del <- apply(vcf_num_del, 2, function(x) length(which(x=="1")) + 2*length(which(x=="2")))
meta <- cbind(meta, sum_hethom_sites_del)
box_TEs_per_sample_hethom_del <- ggplot(data = meta, aes(x = pop_id, y = sum_hethom_sites_del, color=clade))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  ylim(c(0,60000))+
  scale_color_manual(values=cladePalette)+
  labs(title="Total number of polymorphic insertions (shared with reference)",x = "Population", y = "Sum of insertions (per sample)")+
  theme_bw(base_size = 26)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_MEI_sum_hethom_sites_per_sample_del.pdf", width=18, height=12)
box_TEs_per_sample_hethom_del
dev.off()

#######################################################
#######################################################
#### ok now on total by sequencing depth
seq_depth_by_id <- ggplot(data = meta, aes(x = seq_depth, y = sum_pol_per_sample, color=pop_id))+
  geom_point(size=5)+
  theme_bw()+
  scale_color_manual(values=popPalette)+
  labs(title="Relationship between sequencing depth and number of polymorphic insertions identified",
  x = "Sequencing depth", y = "Number of insertions per sample")+
  theme_bw(base_size = 28)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_seq_depth_by_id_pol.pdf", width=18, height=18)
seq_depth_by_id
dev.off()

seq_depth_by_id_het <- ggplot(data = meta, aes(x = seq_depth, y = sum_het_sites_all, color=pop_id))+
  geom_point(size=5)+
  theme_bw()+
  scale_color_manual(values=popPalette)+
  labs(title="Relationship between coverage and number of polymorphic sites identified",
  x = "Sequencing depth", y = "Number of heterzoygous sites per sample")+
  theme_bw(base_size = 28)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-2_seq_depth_by_id_het_sites.pdf", width=18, height=18)
seq_depth_by_id_het
dev.off()

#### add mean sequencing depth
meta <- meta[, mean_seq_depth := mean(seq_depth), by=pop_id]
meta <- meta[, seq_depth_diff := seq_depth - mean_seq_depth]

meta <- meta[, mean_sum_ins := mean(sum_ins_per_sample), by=pop_id]
meta <- meta[, sum_ins_diff := sum_ins_per_sample - mean_sum_ins]
meta <- meta[, mean_sum_all := mean(sum_pol_per_sample), by=pop_id]
meta <- meta[, sum_all_diff := sum_pol_per_sample - mean_sum_all]
meta <- meta[, mean_sum_del := mean(sum_del_per_sample), by=pop_id]
meta <- meta[, sum_del_diff := sum_del_per_sample - mean_sum_del]

sum_ins_diff_by_id <- ggplot(data = meta, aes(x = seq_depth, y = sum_ins_diff, color=pop_id))+
  geom_point(size=5)+
  theme_bw()+
  labs(title="Influence of coverage on TE identification",
  x = "sequencing depth", y = "TE het insertions (sample) - TE het insertions (mean for population)")+
  theme_bw(base_size = 26)+
  scale_color_manual(values=popPalette)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black")+
  geom_smooth(data=subset(meta, seq_depth < 25), method='lm',formula=y~x,se=F, color="blue")
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_sum_ins_diff_by_id.pdf", width=18, height=16)
sum_ins_diff_by_id
dev.off()

sum_del_diff_by_id <- ggplot(data = meta, aes(x = seq_depth, y = sum_del_diff, color=pop_id))+
  geom_point(size=5)+
  theme_bw()+
  labs(title="Influence of coverage on TE identification",
  x = "sequencing depth", y = "TE het deletions (sample) - TE het deletions (mean for population)")+
  theme_bw(base_size = 26)+
  scale_color_manual(values=popPalette)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black")+
  geom_smooth(data=subset(meta, seq_depth < 25), method='lm',formula=y~x,se=F, color="blue")
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_sum_del_diff_by_id.pdf", width=18, height=16)
sum_del_diff_by_id
dev.off()

sum_all_diff_by_id <- ggplot(data = meta, aes(x = seq_depth, y = sum_all_diff, color=pop_id))+
  geom_point(size=5)+
  theme_bw()+
  labs(title="Influence of coverage on TE identification",
  x = "sequencing depth", y = "Nr het sites (sample) - nr het sites (population mean)")+
  theme_bw(base_size = 26)+
  scale_color_manual(values=popPalette)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black")+
  geom_smooth(data=subset(meta, seq_depth < 25), method='lm',formula=y~x,se=F, color="blue")
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_sum_all_diff_by_id.pdf", width=18, height=16)
sum_all_diff_by_id
dev.off()

##### is there a relationship for the number of heterozygous sites
meta <- meta[, mean_sum_hets_all := mean(sum_het_sites_all), by=pop_id]
meta <- meta[, sum_het_diff_all := sum_het_sites_all - mean_sum_hets_all]
sum_het_diff_by_id_reg_all <- ggplot(data = meta, aes(x = seq_depth, y = sum_het_diff_all, color=pop_id))+
  geom_point(size=5)+
  theme_bw()+
  labs(title="Influence of coverage on TE identification",
  x = "sequencing depth", y = "Nr het polymorphisms (sample) - nr het polymorphisms (population mean)")+
  theme_bw(base_size = 24)+
  scale_color_manual(values=popPalette)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black")+
  geom_smooth(data=subset(meta, seq_depth < 25), method='lm',formula=y~x,se=F, color="blue")+
  #stat_cor(aes(label = after_stat(rr.label)), color = "black", geom = "label")
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-11_sum_het_diff_by_id_all.pdf", width=18, height=16)
sum_het_diff_by_id_reg_all
dev.off()

meta <- meta[, mean_sum_hets_ins := mean(sum_het_sites_ins), by=pop_id]
meta <- meta[, sum_het_diff_ins := sum_het_sites_ins - mean_sum_hets_ins]
sum_het_diff_by_id_reg_ins <- ggplot(data = meta, aes(x = seq_depth, y = sum_het_diff_ins, color=pop_id))+
  geom_point(size=5)+
  theme_bw()+
  labs(title="Influence of coverage on TE identification",
  x = "sequencing depth", y = "Nr het insertions (sample) - nr het insertions (population mean)")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method='lm', formula= y~x, se=FALSE, color="black")+
  geom_smooth(data=subset(meta, seq_depth < 25), method='lm',formula=y~x,se=F, color="blue")+
  #stat_cor(aes(label = after_stat(rr.label)), color = "black", geom = "label")
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-2_sum_het_diff_by_id_ins.pdf", width=18, height=16)
sum_het_diff_by_id_reg_ins
dev.off()

################################################################################
################################################################################
################################################################################
################################################################################
##### sum up for different orders
sum_hets_by_order <- vcf[, lapply(.SD, function(x) length(which(x=="1"))), by=order, .SDcols=meta_ids]
sum_hets_by_superfamily <- vcf[, lapply(.SD, function(x) length(which(x=="1"))), by=superfamily, .SDcols=meta_ids]
sum_hets_by_family <- vcf[, lapply(.SD,function(x) length(which(x=="1"))), by=family, .SDcols=meta_ids]

t_sum_by_order <- as.data.frame(t(sum_hets_by_order[,-1]))
colnames(t_sum_by_order) <- sum_hets_by_order$order
t_sum_by_superfamily <- as.data.frame(t(sum_hets_by_superfamily[,-1]))
colnames(t_sum_by_superfamily) <- sum_hets_by_superfamily$superfamily
t_sum_by_family <- as.data.frame(t(sum_hets_by_family[,-1]))
colnames(t_sum_by_family) <- sum_hets_by_family$family

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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_MEI_sum_ins_per_sample.pdf", width=22, height=16)
box_TEs_per_sample
dev.off()

meta_sums_family <- cbind(meta_sub, t_sum_by_family)
meta_sums_family_melt <- melt(meta_sums_family)
meta_family <- meta_sums_family_melt[, mean_pop := mean(value), by=pop_id]
meta_family <- meta_sums_family_melt[, mean_fam := mean(value), by=variable]
meta_family <- meta_sums_family_melt[, mean_fam_pop := mean(value), by=.(variable, pop_id)]
meta_family_mean <- unique(meta_family[,c("clade", "pop_id", "variable", "mean_fam", "mean_pop", "mean_fam_pop"), with=FALSE])
mean_diff_fam <- meta_family_mean[, diff := mean_fam_pop - mean_fam]
mean_diff_fam <- meta_family_mean[, diff_percent := abs(mean_fam_pop - mean_fam) / mean_fam * 100]

mean_diff_fam[, c("superfamily", "number") := tstrsplit(variable, "-", fixed=FALSE, keep=c(1,2))]
mean_diff_fam <- mean_diff_fam[,superfamily := as.factor(superfamily)]
mean_piggy <- mean_diff_fam[superfamily == "PiggyBac"]

box_TEs_per_family <- ggplot(data = mean_piggy, aes(x = pop_id, y = diff_percent))+
  geom_point(aes(color=variable, size=5))+
  theme_bw()+
  facet_wrap(~variable)+
  labs(title="Distribution of PiggyBack families",x = "Population", y = "% number of mean insertions in the family relative to average number for the family")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-6_MEI_percent_diff_family_piggyback.pdf", width=22, height=16)
box_TEs_per_family
dev.off()

#####################################################################
#####################################################################
### Okay, how many het sites on average does each population have?
####################################################################
#####################################################################
#####################################################################

###### Not a very smart way to do all this
######## split vcfs by orders 
sum_ap <- as.data.table(colSums(vcf[,..meta_ap]))
sum_ask <- as.data.table(colSums(vcf[,..meta_as_kingiri]))
sum_asm <- as.data.table(colSums(vcf[,..meta_as_masoko]))
sum_crhoa <- as.data.table(colSums(vcf[,..meta_crhoa]))
sum_cch <- as.data.table(colSums(vcf[,..meta_cch]))
sum_cvm <- as.data.table(colSums(vcf[,..meta_cv_malombe]))
sum_cvs <- as.data.table(colSums(vcf[,..meta_cv_swarm]))
sum_czeb <- as.data.table(colSums(vcf[,..meta_czeb]))
sum_diplo <- as.data.table(colSums(vcf[,..meta_dip]))
sum_fr <- as.data.table(colSums(vcf[,..meta_fr]))
sum_fuel <- as.data.table(colSums(vcf[,..meta_fuel]))
sum_trew <- as.data.table(colSums(vcf[,..meta_trew]))
sum_mzeb <- as.data.table(colSums(vcf[,..meta_mzeb]))
sum_oto <- as.data.table(colSums(vcf[,..meta_oto]))
sum_rham <- as.data.table(colSums(vcf[,..meta_rham]))

sum_ap <- sum_ap[,population := "Alticorpus peterdaviesi.Cape_Maclear"]
sum_ask <- sum_ask[,population := "Astatotilapia calliptera.Lake_Kingiri"]
sum_asm <- sum_asm[,population := "Astatotilapia calliptera.Lake_Masoko"]
sum_crhoa <- sum_crhoa[,population := "Chilotilapia rhoadesii"]
sum_cch <- sum_cch[,population := "Copadichromis chrysonotus.Lake_Malombe"]
sum_cvm <- sum_cvm[,population := "Copadichromis virginalis.Lake_Malombe"]
sum_cvs <- sum_cvs[,population := "Copadichromis virginalis.Southwest_arm"]
sum_czeb <- sum_czeb[,population := "Cynotilapia zebroides.Cape_Maclear"]
sum_diplo <- sum_diplo[,population := "Diplotaxodon limnothrissa.Southwest_arm"]
sum_fr <- sum_fr[,population := "Fossorochromis rostratus.Lake_Malombe"]
sum_fuel <- sum_fuel[,population := "Labeotropheus fuelleborni.Chilumba"]
sum_trew <- sum_trew[,population := "Labeotropheus trewavasae.Chilumba"]
sum_mzeb <- sum_mzeb[,population := "Maylandia zebra.Cape_Maclear"]
sum_oto <- sum_oto[,population := "Otopharynx argyrosoma.Southeast_arm"]
sum_rham <- sum_rham[,population := "Rhamphochromis longiceps"]

sums_per_species <- rbind(sum_ap, sum_ask, sum_asm, sum_crhoa, sum_cch, sum_cvm, sum_cvs, sum_czeb, 
sum_diplo, sum_fr, sum_fuel, sum_trew, sum_mzeb, sum_oto, sum_rham)
sums_per_species <- sums_per_species[,population := as.factor(population)]

box_TEs_per_sample_pop <- ggplot(data = sums_per_species, aes(x = population, y = V1, color=population))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  labs(title="Number of polymorphic insertions per sample",x = "Population", y = "Number of insertions per sample")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-28_MEI_sum_ins_per_sample_by_pop_insdel.pdf", width=18, height=12)
box_TEs_per_sample_pop
dev.off()

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

