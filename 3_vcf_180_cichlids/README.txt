##### README
##### 2023-03-16 
##### Barbara Walkowiak 
##### NST2GN Part II Project 

##### PACKAGE VERSIONS / SOFTWARE
Rversion 4.2.2
ggpmisc_0.5.2      ggpp_0.5.1         RColorBrewer_1.1-3 factoextra_1.0.7  
FactoMineR_2.7     ggrepel_0.9.3      tibble_3.1.8       dplyr_1.1.0       
plyr_1.8.8         ggplot2_3.4.1      gridExtra_2.3      data.table_1.14.6 


In this folder, I include scripts which I used the vcf file on the entire callset to generate a subsetted vcf containing polymorphic insertions (shared with and absent form the refernece). Input and output files can be found in NST2GN/Data/vcf folder

##### Scripts
23-2-16_1_generate_subsetted_vcf.r 
	## subset the main vcf file for insertions in the 180 cichlids I am looking at
23-2-16_2_generate_subsetted_vcf_info.r 
	## split INFO columns in the subsetted vcf to obtain data about TE classification 
23-2-22_3_generate_subsetted_vcf_correct_ins.r 
	## ensure each insertion is present in at least one sample (genotypes "0/1", "./1", "1/1")
23-2-24_4_merge_vcf_insertions_deletions.r 
	## merge insertion data together with deletion data (include insertions shared with the reference) 
