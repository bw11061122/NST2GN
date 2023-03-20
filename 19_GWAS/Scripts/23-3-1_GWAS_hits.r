#### 2023-03-01

# trait    insID	dist2Gene	VEP	geneID	geneSymbol
# colour    chr7:41166580-41166583_O	0	intron	ENSACLG00000014805	ntrk3a
# colour	chr3:8249800-8249807_O	0	intron	ENSACLG00000018516	nlgn2a
# Caudal peduncle depth	chr22:13852670-13852671_O	0	intron	ENSACLG00000013131	trps1

#### Can I identify whether there are any ME insertions which are next to these genes?

### load libraries 
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

### meta
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
meta <- meta[, pop_id := name_loc][(genus %in% c("Rhamphochromis", "Chilotilapia")) | is.na(name_loc), pop_id := names][]
meta <- meta[, pop_id:=as.factor(pop_id)] 
upd.cols = sapply(meta, is.factor)
meta <- meta[, names(meta)[upd.cols] := lapply(.SD, factor), .SDcols = upd.cols] ## this gets rid of levels which are absent in the dt 
meta_ids <- unlist(meta[,primary_id]) ### this has all the IDs that you are interested in

### vcf
vcf <- fread("23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-2-27_vcf_ins_del_all_final.txt", 
    sep = "\t", header=T)
vcf <- vcf[, order:=as.factor(order)]
vcf <- vcf[, family:=as.factor(family)]
vcf <- vcf[, superfamily:=as.factor(superfamily)]
upd.cols.vcf = sapply(vcf, is.factor)
vcf <- vcf[, names(vcf)[upd.cols.vcf] := lapply(.SD, factor), .SDcols = upd.cols.vcf] ## this gets rid of levels which are absent in the dt 
setnames(vcf, "#CHROM", "chromosome")

###### okay so what is the average distance between insertions on my chromosomes
vcf_ordered <- vcf[order(chromosome, POS),]
vcf_ordered <- vcf_ordered[ , difference_between_TE := (POS - shift(POS)), by = chromosome]

#### okay maybe I want to have the sum of 0_0, 0_1, 1_1, and missing ones etc for each insertion
vcf_num <- vcf[,c(meta_ids), with=FALSE]
diff_bn_TE <- unlist(vcf_ordered[, difference_between_TE])
max(diff_bn_TE, na.rm=TRUE)
min(diff_bn_TE, na.rm=TRUE)
mean(diff_bn_TE, na.rm=TRUE)
median(diff_bn_TE, na.rm=TRUE)

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-1_distance_between_insertions.pdf")
hist(vcf_ordered[, difference_between_TE], breaks=100)
dev.off()

##### how many insertions do I have per chromosome
sum_ins_per_chrom <- vcf[, .N, by=chromosome]
sum_ins_per_chrom <- sum_ins_per_chrom[chromosome %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
"chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
"chr20", "chr21", "chr22")]
chromosome_levels <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
chr_sizes <- fread("23-3-1_chrom.sizes.txt", header=F)
setnames(chr_sizes, "V1", "chromosome")
setnames(chr_sizes, "V2", "chromosome_size")
setkey(chr_sizes, chromosome)
setkey(sum_ins_per_chrom, chromosome)
sum_ins_chrom <- chr_sizes[sum_ins_per_chrom]
sum_ins_chrom <- sum_ins_chrom[, TE_density_per_chrom := N / chromosome_size]
sum_ins_chrom <- sum_ins_chrom[, TE_density_per_kb := TE_density_per_chrom * 1000]

col <- "chromosome"
sum_ins_chrom[, chromosome := factor(get(col), levels = chromosome_levels)]
TE_per_chrom <- ggplot(data = sum_ins_chrom, aes(x = chromosome, y = N, fill=chromosome))+
  geom_bar(stat="identity")+
  theme_bw()+
  labs(title="Number of polymorphisms per chromosome",x = "chromosome", y = "Nr insertions")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-1_nr_insertions_per_chrom.pdf")
TE_per_chrom
dev.off()

#### average density per chromosome
TE_per_kb <- ggplot(data = sum_ins_chrom, aes(x = chromosome, y = TE_density_per_kb, fill=chromosome))+
  geom_bar(stat="identity")+
  theme_bw()+
  labs(title="Number of polymorphisms per chromosome",x = "chromosome", y = "Nr insertions / kb")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-1_nr_insertions_per_kb.pdf", width=16, height=10)
TE_per_kb
dev.off()

#### average distance between TE per chromosome
mean_dist_TE <- vcf_ordered[, mean_dist_per_chr := mean(difference_between_TE, na.rm=TRUE), by=chromosome]
mean_dist_TE <- mean_dist_TE[,c("chromosome", "difference_between_TE", "mean_dist_per_chr"), with=FALSE]
mean_dist_TE <- mean_dist_TE[chromosome %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
"chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
"chr20", "chr21", "chr22")]
chromosome_levels <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
col <- "chromosome"
mean_dist_TE[, chromosome := factor(get(col), levels = chromosome_levels)]
mean_dist_TE[,mean_dist_plot := round(mean_dist_per_chr,1)]
TE_av_dist <- ggplot(data = mean_dist_TE, aes(x = chromosome, y = difference_between_TE, color=chromosome))+
  geom_boxplot()+
  ylim(c(0,2000))+
  labs(title="Average distance between neighbouring polymorphisms",x = "chromosome", y = "Distance between neighbouring TE")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-2_mean_dist_TE.pdf", width=16, height=10)
TE_av_dist
dev.off()

##### okay the jitter is the most uninformative thing I could possibly get 
TE_av_dist_jitter <- ggplot(data = mean_dist_TE, aes(x = chromosome, y = difference_between_TE, color=chromosome))+
  geom_boxplot()+
  geom_jitter()+
  ylim(c(0,2000))+
  labs(title="Average distance between neighbouring polymorphisms",x = "chromosome", y = "Distance between neighbouring TE")+
  theme_bw(base_size = 24)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-2_mean_dist_TE_jitter.pdf", width=16, height=10)
TE_av_dist_jitter
dev.off()

# trait    insID	dist2Gene	VEP	geneID	geneSymbol
# colour    chr7:41166580-41166583_O	0	intron	ENSACLG00000014805	ntrk3a
# colour	chr3:8249800-8249807_O	0	intron	ENSACLG00000018516	nlgn2a
# Caudal peduncle depth	chr22:13852670-13852671_O	0	intron	ENSACLG00000013131	trps1

#### what are the coordinates of these genes? (literally gene start and end)
#### I am taking this from the Ensembl genome browser, the genome assembly is called easter happy
#### ntrk3a 7: 41,141,843-41,234,380
#### nlgn2a 3: 8,167,129-8,404,685
#### trps1 22: 13,818,172-13,914,023 
vcf_chr3 <- vcf[chromosome == "chr3",]
vcf_chr7 <- vcf[chromosome == "chr7",]
vcf_chr22 <- vcf[chromosome == "chr22",]

ntrk_region <- c((41141843-5000): (41234380+5000))
nlgn_region <- c((8167129-5000): (8404685+5000))
trps_region <- c((13818172-5000): (13914023+5000))
vcf_ntrk <- vcf[chromosome == "chr7" & POS %in% c(ntrk_region),]
vcf_nlgn <- vcf[chromosome == "chr3" & POS %in% c(nlgn_region),]
vcf_trps <- vcf[chromosome == "chr22" & POS %in% c(trps_region),]

#### Are there any in exons? Which populations have them?
#### classify them in introns, exons, UTRs 
#### what is the allele frequency of each?
#### got this file from Pio

#### I can first get sums for each populations
meta_ap <- unlist(meta[name_loc=="Alticorpus peterdaviesi.Cape_Maclear",primary_id])
meta_as_kingiri <- unlist(meta[name_loc=="Astatotilapia calliptera.Lake_Kingiri",primary_id])
meta_as_masoko <- unlist(meta[name_loc=="Astatotilapia calliptera.Lake_Masoko",primary_id])
meta_crhoa <- unlist(meta[names=="Chilotilapia rhoadesii",primary_id])
meta_cch <- unlist(meta[name_loc=="Copadichromis chrysonotus.Lake_Malombe",primary_id])
meta_cv_malombe <- unlist(meta[name_loc=="Copadichromis virginalis.Lake_Malombe",primary_id])
meta_cv_swarm <- unlist(meta[name_loc=="Copadichromis virginalis.Southwest_arm",primary_id])
meta_czeb <- unlist(meta[name_loc=="Cynotilapia zebroides.Cape_Maclear",primary_id])
meta_dip <- unlist(meta[name_loc=="Diplotaxodon limnothrissa.Southwest_arm",primary_id])
meta_fr <- unlist(meta[name_loc=="Fossorochromis rostratus.Lake_Malombe",primary_id])
meta_fuel <- unlist(meta[name_loc=="Labeotropheus fuelleborni.Chilumba",primary_id])
meta_trew <- unlist(meta[name_loc=="Labeotropheus trewavasae.Chilumba",primary_id])
meta_mzeb <- unlist(meta[name_loc=="Maylandia zebra.Cape_Maclear",primary_id])
meta_oto <- unlist(meta[name_loc=="Otopharynx argyrosoma.Southeast_arm",primary_id])
meta_rham <- unlist(meta[names=="Rhamphochromis longiceps",primary_id])

sum_ap_ntrk <- apply(vcf_ntrk[,c(meta_ap),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_ask_ntrk <- apply(vcf_ntrk[,c(meta_as_kingiri),with=FALSE], 1, function(x) length(which(x=="1")))
sum_asm_ntrk <- apply(vcf_ntrk[,c(meta_as_masoko),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_crh_ntrk <- apply(vcf_ntrk[,c(meta_crhoa),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_cch_ntrk <- apply(vcf_ntrk[,c(meta_cch),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_cvm_ntrk <- apply(vcf_ntrk[,c(meta_cv_malombe),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_cvs_ntrk <- apply(vcf_ntrk[,c(meta_cv_swarm),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_czeb_ntrk <- apply(vcf_ntrk[,c(meta_czeb),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_dip_ntrk <- apply(vcf_ntrk[,c(meta_dip),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_fr_ntrk <- apply(vcf_ntrk[,c(meta_fr),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_fuel_ntrk <- apply(vcf_ntrk[,c(meta_fuel),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_trew_ntrk <- apply(vcf_ntrk[,c(meta_trew),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_mzeb_ntrk <- apply(vcf_ntrk[,c(meta_mzeb),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_oto_ntrk <- apply(vcf_ntrk[,c(meta_oto),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_rham_ntrk <- apply(vcf_ntrk[,c(meta_rham),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_total_ntrk <- apply(vcf_ntrk[,c(meta_ids),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
vcf_ntrk <- cbind(vcf_ntrk, sum_ap_ntrk, sum_ask_ntrk, sum_asm_ntrk, sum_crh_ntrk, sum_cch_ntrk, sum_cvm_ntrk, 
sum_cvs_ntrk, sum_czeb_ntrk, sum_dip_ntrk, sum_fr_ntrk, sum_fuel_ntrk, sum_trew_ntrk, sum_mzeb_ntrk, sum_oto_ntrk, 
sum_rham_ntrk, sum_total_ntrk)
vcf_ntrk_filtered <- vcf_ntrk[sum_total_ntrk >=24,]
vcf_ntrk_filtered <- vcf_ntrk_filtered[,gene := "ntrk"]

sum_ap_nlgn <- apply(vcf_nlgn[,c(meta_ap),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_ask_nlgn <- apply(vcf_nlgn[,c(meta_as_kingiri),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_asm_nlgn <- apply(vcf_nlgn[,c(meta_as_masoko),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_crh_nlgn <- apply(vcf_nlgn[,c(meta_crhoa),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_cch_nlgn <- apply(vcf_nlgn[,c(meta_cch),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_cvm_nlgn <- apply(vcf_nlgn[,c(meta_cv_malombe),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_cvs_nlgn <- apply(vcf_nlgn[,c(meta_cv_swarm),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_czeb_nlgn <- apply(vcf_nlgn[,c(meta_czeb),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_dip_nlgn <- apply(vcf_nlgn[,c(meta_dip),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_fr_nlgn <- apply(vcf_nlgn[,c(meta_fr),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_fuel_nlgn <- apply(vcf_nlgn[,c(meta_fuel),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_trew_nlgn <- apply(vcf_nlgn[,c(meta_trew),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_mzeb_nlgn <- apply(vcf_nlgn[,c(meta_mzeb),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_oto_nlgn <- apply(vcf_nlgn[,c(meta_oto),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_rham_nlgn <- apply(vcf_nlgn[,c(meta_rham),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_total_nlgn <- apply(vcf_nlgn[,c(meta_ids),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
vcf_nlgn <- cbind(vcf_nlgn, sum_ap_nlgn, sum_ask_nlgn, sum_asm_nlgn, sum_crh_nlgn, sum_cch_nlgn, sum_cvm_nlgn, 
sum_cvs_nlgn, sum_czeb_nlgn, sum_dip_nlgn, sum_fr_nlgn, sum_fuel_nlgn, sum_trew_nlgn, sum_mzeb_nlgn, sum_oto_nlgn, 
sum_rham_nlgn, sum_total_nlgn)
vcf_nlgn_filtered <- vcf_nlgn[sum_total_nlgn >=24,]
vcf_nlgn_filtered <- vcf_nlgn_filtered[,gene := "nlgn"]

sum_ap_trps <- apply(vcf_trps[,c(meta_ap),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_ask_trps <- apply(vcf_trps[,c(meta_as_kingiri),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_asm_trps <- apply(vcf_trps[,c(meta_as_masoko),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_crh_trps <- apply(vcf_trps[,c(meta_crhoa),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_cch_trps <- apply(vcf_trps[,c(meta_cch),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_cvm_trps <- apply(vcf_trps[,c(meta_cv_malombe),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_cvs_trps <- apply(vcf_trps[,c(meta_cv_swarm),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_czeb_trps <- apply(vcf_trps[,c(meta_czeb),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_dip_trps <- apply(vcf_trps[,c(meta_dip),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_fr_trps <- apply(vcf_trps[,c(meta_fr),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_fuel_trps <- apply(vcf_trps[,c(meta_fuel),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_trew_trps <- apply(vcf_trps[,c(meta_trew),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_mzeb_trps <- apply(vcf_trps[,c(meta_mzeb),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_oto_trps <- apply(vcf_trps[,c(meta_oto),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_rham_trps <- apply(vcf_trps[,c(meta_rham),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
sum_total_trps <- apply(vcf_trps[,c(meta_ids),with=FALSE], 1, function(x) length(which(x=="1" | x=="2")))
vcf_trps <- cbind(vcf_trps, sum_ap_trps, sum_ask_trps, sum_asm_trps, sum_crh_trps, sum_cch_trps, sum_cvm_trps, 
sum_cvs_trps, sum_czeb_trps, sum_dip_trps, sum_fr_trps, sum_fuel_trps, sum_trew_trps, sum_mzeb_trps, sum_oto_trps, 
sum_rham_trps, sum_total_trps)
vcf_trps_filtered <- vcf_trps[sum_total_trps >=24,]
vcf_trps_filtered <- vcf_trps_filtered[,gene := "trps"]

setnames(vcf_trps_filtered, c("sum_ap_trps", "sum_ask_trps", "sum_asm_trps", "sum_crh_trps", "sum_cch_trps", "sum_cvm_trps", 
"sum_cvs_trps", "sum_czeb_trps", "sum_dip_trps", "sum_fr_trps", "sum_fuel_trps", "sum_trew_trps", "sum_mzeb_trps", "sum_oto_trps", 
"sum_rham_trps", "sum_total_trps"), c("sum_ap", "sum_ask", "sum_asm", "sum_chr", "sum_cch", "sum_cvm", "sum_cvs", "sum_czeb",
"sum_dip", "sum_fr", "sum_fuel", "sum_trew", "sum_mzeb", "sum_oto", "sum_rham", "sum_total"))

setnames(vcf_nlgn_filtered, c("sum_ap_nlgn", "sum_ask_nlgn", "sum_asm_nlgn", "sum_crh_nlgn", "sum_cch_nlgn", "sum_cvm_nlgn", 
"sum_cvs_nlgn", "sum_czeb_nlgn", "sum_dip_nlgn", "sum_fr_nlgn", "sum_fuel_nlgn", "sum_trew_nlgn", "sum_mzeb_nlgn", "sum_oto_nlgn", 
"sum_rham_nlgn", "sum_total_nlgn"), c("sum_ap", "sum_ask", "sum_asm", "sum_chr", "sum_cch", "sum_cvm", "sum_cvs", "sum_czeb",
"sum_dip", "sum_fr", "sum_fuel", "sum_trew", "sum_mzeb", "sum_oto", "sum_rham", "sum_total"))

setnames(vcf_ntrk_filtered, c("sum_ap_ntrk", "sum_ask_ntrk", "sum_asm_ntrk", "sum_crh_ntrk", "sum_cch_ntrk", "sum_cvm_ntrk", 
"sum_cvs_ntrk", "sum_czeb_ntrk", "sum_dip_ntrk", "sum_fr_ntrk", "sum_fuel_ntrk", "sum_trew_ntrk", "sum_mzeb_ntrk", "sum_oto_ntrk", 
"sum_rham_ntrk", "sum_total_ntrk"), c("sum_ap", "sum_ask", "sum_asm", "sum_chr", "sum_cch", "sum_cvm", "sum_cvs", "sum_czeb",
"sum_dip", "sum_fr", "sum_fuel", "sum_trew", "sum_mzeb", "sum_oto", "sum_rham", "sum_total"))

vcf_filtered_gwas <- rbind(vcf_trps_filtered, vcf_nlgn_filtered, vcf_ntrk_filtered)
write.table(vcf_filtered_gwas, "23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-2_vcf_filtered_gwas.txt", 
row.names=F, quote=F, sep="\t")

### How many sites in each sample are heterozygous or homozygous?
sum_0_0_per_sample <- apply(vcf_num, 2, function(x) length(which(x=="0")))
sum_0_1_per_sample <- apply(vcf_num, 2, function(x) length(which(x=="1")))
sum_1_1_per_sample <- apply(vcf_num, 2, function(x) length(which(x=="2")))
sum_0_missing_per_sample <- apply(vcf_num, 2, function(x) length(which(x=="0.5")))
sum_1_missing_per_sample <- apply(vcf_num, 2, function(x) length(which(x=="1.5")))

### How many times is a given TE present in a heterozygous vs homozygous state?
sum_0_0_per_te <- apply(vcf_num, 1, function(x) length(which(x=="0")))
sum_0_1_per_te <- apply(vcf_num, 1, function(x) length(which(x=="1")))
sum_1_1_per_te <- apply(vcf_num, 1, function(x) length(which(x=="2")))
sum_0_missing_per_te <- apply(vcf_num, 1, function(x) length(which(x=="0.5")))
sum_1_missing_per_te <- apply(vcf_num, 1, function(x) length(which(x=="1.5")))



