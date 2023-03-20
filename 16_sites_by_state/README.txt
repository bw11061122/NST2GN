##### README
##### 2023-03-16 
##### Barbara Walkowiak 
##### NST2GN Part II Project 

##### PACKAGE VERSIONS / SOFTWARE
Rversion 4.2.2
ggpmisc_0.5.2      ggpp_0.5.1         RColorBrewer_1.1-3 factoextra_1.0.7  
FactoMineR_2.7     ggrepel_0.9.3      tibble_3.1.8       dplyr_1.1.0       
plyr_1.8.8         ggplot2_3.4.1      gridExtra_2.3      data.table_1.14.6 

bcftools 1.15.1    

megane_v1.0.1.beta.sif 

###### DATA
## FOLDER: Data
### Metadata: 
    cichlid_callset_metadata.xlsx ## full metadata
    23-2-16_cichlid_meta_standardized_dataset.csv ## final metadata used 
    23-3-10_180cichlid_geo_distirbution_map ## screenshot of map I used on the presentation to show geogrpahical distribution
### Vcf
    cichlid_all_MEI_biallelic.vcf.gz ## raw output from MEGANE (insertions)
    cichlid_all_MEA_biallelic.vcf.gz ## raw output from MEGANE (deletions)
    cichlid_all_MEI_biallelic_no_meta_encoded.vcf.gz ## output from MEGANE (removed metadata, encoded genotypes)
    cichlid_all_MEA_biallelic_no_meta_encoded.vcf.gz ## output from MEGANE (removed metadata, encoded genotypes)
### TE annotation library
    lib2_renamed_nodups_fixedunkowns_plusgraph_tab.rep

## File location on the cluster:
## MEGANE 
## Genome (AstCal1.2)
~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa) 
## paths to all cichlid files to generate the vcf (fofn from Bettina Fischer)
~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/data
## Generated vcf files 
~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/23-2-1_MEGANE_all_cichilds/vcf_for_phasing
## TE annotation libarary (MWCichlidTE-3): 
~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/Lib3

## CORE ANALYSIS
## FOLDER: Scripts
## in each folder, I included scripts which are required for a given step and output figure (if applicable)
## if the folder described the generation of a figure which required data files other than those in the data folder, I included these in the folder 

##### MEGANE pipeline (vcf calling)
##### FOLDER: 1_MEGANE
MEGANE step 0 megane_step0_a_calliptera.sh
Modify library 23-1-24_get_correct_lib2.sh
MEGANE step 1 prep 23-1-23_megane_launcher_step1_prep_a_calliptera.sh
MEGANE step 1 23-1-24_test_megane_launch.sh
MEGANE step 2 23-1-26_megane_step2.sh
MEGANE step 3 23-1-27_megane_step3.sh

##### Choice of metadata
##### FOLDER: 2_metadata_180_cichlids
23-2-16_standardized_meta_final.R 
## in this file, I chose the final set of 180 cichlids (12 individuals x 15 populations) together with associated metadata

##### Generation of the subsetted vcf file
##### FOLDER: 3_vcf_180_cichlids
23-2-16_1_generate_subsetted_vcf.r ## subset the main vcf file for insertions in the 180 cichlids I am looking at
23-2-16_2_generate_subsetted_vcf_info.r ## split INFO columns in the subsetted vcf to obtain data about TE classification 
23-2-22_3_generate_subsetted_vcf_correct_ins.r ## ensure each insertion is present in at least one sample (genotypes "0/1", "./1", "1/1")
23-2-24_4_merge_vcf_insertions_deletions.r ## merge insertion data together with deletion data

##### Basic stats metadata
##### FOLDER: 4_metadata_180_cichlids_analysis
23-3-9_metadata_analysis.r
## basic analysis of subsetted metadata (# locations, populations, species)

##### Anlysis of the distribution of TE insertion SITES (Fig. 1)
##### FOLDER: 5_site_analysis
23-3-9_vcf_analysis.r

##### Analysis of the distribution of TE INSERTIONS (fig. 2)
##### FOLDER: 6_insertion_analysis
23-2-27_sums.r

##### PCA (fig. 3)
##### FOLDER: 7_PCA
23-2-22_PCA_on_filtered_final.r

##### Heatmaps (fig. 4)
##### FOLDER: 8_Heatmaps 
## Note: scripts run locally (not on the cluster)
23-3-10_pheatmap.R
23-3-11_heatmap_insertions.R
23-3-11_heatmap_deletions.R

##### Minor allele frequency (fig. 5, fig. 6)
##### FOLDER: 9_MAF  
23-3-9_minor_allele_freq.r

##### Minor allele frequency (MAF) heatmaps
##### FOLDER: 10_MAF_heatmaps
## Note: scripts run locally (not on the cluster)
23-3-14_heatmap_maf_05_by_fam_by_pop.R
23-3-13_heatmap_maf.R
23-3-13_maf_analysis_for_heatmaps.R ## this script should have complete code to generate both tables used for heatmaps 

##### Families which contribute most insertions across populations (table 1)
##### FOLDER: 11_families_top5_present_in_all_populations
23-3-6_differences_across_families.r
23-3-11_differences_bn_fam_ranking.r

##### Families which are not present in across populations (table 2)
##### FOLDER: 12_families_not_present_in_all_populations
23-3-12_top5_families.R

##### Kruskall-Wallist test (fig. 8)
##### FOLDER: 14_Kruskal-Wallis_test
23-3-8_differences_between_families_v2.r

##### Relationship between number of insertions and % shared in all families (fig. 8)
23-3-9_percent_ins_contributing.r

##### Comparison to SNPs (PCA - figure 3, MAF - figure 6)
23-3-12_snps.sh ## how to get subsets of biallelic SNPs 
(using https://github.com/MoritzBlumer/inversion_scripts/tree/main/windowed_pca)

## SUPPLEMENTARY FIGURES
##### Count of different genotypes (fig. S1)
23-3-11_counts_sites_by_state.r

##### Sequencing depth vs number of insertions identified (fig. S2)


##### PCA for 4 main TE orders (fig. S3)
23-3-8_PCA_orders_superfamilies.r

##### Additional clustering analysis (UMAP, tSNE) (fig. S4)
23-2-17_umap_tSNE.r

##### Kruskal-Wallist test: number of insertions vs significance (fig. S5)
23-3-9_percent_ins_contributing.r

##### GWAS hits 
23-3-6_gwas_update.r

##### GWAS hits
23-3-12_filtering_count.sh (filtering)

