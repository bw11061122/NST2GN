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

In this folder, I generate SNP files for comparison with the TE data.
I generated 3 SNP files which can be found in the /Data subfolder:
23-3-12_vcf_snps_set180_encoded.vcf.gz 
is the first file with 517436 polymorphic SNPs generated from chr4
23-3-15_snps_subset_for_comparison_v2.txt
is the second file (103487 SNPs from chromosomes: 4, 6, 10, 12, 14 each)
23-3-16_snps_subset_for_comparison_v3.txt
is the second file (103487 SNPs from chromosomes: 2, 5, 9, 16, 19 each) 

The subsetted files for chromosomes can be found on the cluster
~/rds/rds-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/snps

Additional input: list of sample names is located in the NST2GN_FINAL/Data/Metadata folder

Scripts used to analyse these files have names with dates corresponding to the date in the file name:

Note: I would like to thank Moritz Blumer for the help with generating the subsetted SNP vcf files 
(I used https://github.com/MoritzBlumer/inversion_scripts/tree/main/windowed_pca to get the initial files)
