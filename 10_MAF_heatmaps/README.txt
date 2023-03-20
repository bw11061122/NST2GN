##### README
##### 2023-03-16 
##### Barbara Walkowiak 
##### NST2GN Part II Project 

##### PACKAGE VERSIONS / SOFTWARE
Rversion 4.2.2
ggpmisc_0.5.2      ggpp_0.5.1         RColorBrewer_1.1-3 factoextra_1.0.7  
FactoMineR_2.7     ggrepel_0.9.3      tibble_3.1.8       dplyr_1.1.0       
plyr_1.8.8         ggplot2_3.4.1      gridExtra_2.3      data.table_1.14.6 

In this folder, I generate hetmaps of minor allele frequency value and the fraction of sites from each superfamily which is associated with MAF > 0.05 (ie common) across each population. The aggregate MAF data used to generate the heat maps can be found in the Data subfolder. 

Note: files with indication "percent_sites_05" refer to files used for analysis of the number of insertion sites with MAF > 0.05. Other files refer to the analysis of raw MAF values. 

Note: scripts were run locally and the paths correspond to the paths on my PC at the time of doing the project. The input files (metadata and vcf file) were the same as those present on the cluster and the results can be reproduced using these 