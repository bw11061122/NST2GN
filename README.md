# NST2GN Part II Genetics Project
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

This folder contains 20 folders with scripts and associated results and a data folder which together should contain all data that is relevant to my Part II project. In each folder, there is a README file which briefly explains the purpose of the scripts contained in that folder and any additional relevant information. 
In each of the 1-20 folders, I included separate folders for scripts and results and (where relevant) additional data used to generate results. Otherwise, all results should be reproducible with the data included in the main Data folder.  

Brief summary of each folder:
1. Running the MEGANE pipeline to identify polymorphic insertions
2. Choice of 180 samples used for downstream analysis 
3. Subsetting of the general vcf file to extract polymorphic insertions in the chosen 180 samples
4. Analysis of the general features of metadata
5. Analysis of the number of polymorphic TE sites contributed by each TE order across populations
6. Analysis of total numbers of polymorphic insertions across populations
7. Analysis of population structure (PCA)
8. Analysis of the differential patterns of activity of TE superfamilies and families using heatmaps
9. Analysis of minor allele frequency for the entire dataset and for each order 
10. Further analysis of MAF (raw MAF values and % MAF sites > 0.05 in each TE family and superfamily) 
11. Analysis of the distribution of 16 families not detected in all populations
12. Analysis of families contributing most insertions in each population
13. Comparison to equal-sized sets of SNPs 
14. Correlation between coverage of a given sample number of sites identified (quality control) 
15. Kruskal-Wallis test (number of insertions identified in a given family (normalised to a given population) ~ population, by TE family) 

Additional data (not in the report, potentially relevant for me / for the work of the group / supplementary)

16. Analysis of the distribution of different genotypes identified 
17. Population clustering (UMAP, tSNE)
18. Analysis of the relationship between number of insertions and % sites shared
19. Analysis of GWAS hits (from a GWAS performed in the group) 
 
Data: metadata, final vcf files (not available on github due to size), TE annotation library (MWCichlidTE-3)

