##### README
##### 2023-03-16 
##### Barbara Walkowiak 
##### NST2GN Part II Project 

This folder includes the key data used in the project which was used for downstream analysis and to obtain results generated in other folders: 1) Metadata, 2) vcf file of polymorphic TE insertions, 3) TE annotation library (MWCichlidTE-3). 
Note that data on SNPs (used for comparison during the project) is present in the separate folder 13_comparison to SNPs 

## FOLDER: Data
### Metadata: 
    cichlid_callset_metadata.xlsx ## full metadata
    23-2-16_cichlid_meta_standardized_dataset.csv ## final metadata used (180 samples)
    23-3-1_meta_samplesids.txt ## list of sample names (180 samples)
    23-3-10_180cichlid_geo_distirbution_map 
    ## screenshot of map I used on the presentation to show geographical distribution
    23-3-15_meta_table_report.txt (txt format), 23-3-15_report_meta_table1.xlsx (xlsx format)
    ## table included in my final report to show sequencing depth by species by location 
    23-3-1_chrom_sizes.txt ## sizes of chromosomes of the AstCal1.2 genome 
### vcf (these are available in the folder on the cluster, but not on github)
    cichlid_all_MEI_biallelic.vcf.gz ## raw output from MEGANE (insertions)
    cichlid_all_MEA_biallelic.vcf.gz ## raw output from MEGANE (deletions)
    cichlid_all_MEI_biallelic_no_meta_encoded.vcf.gz ## output from MEGANE (removed metadata, encoded genotypes)
    cichlid_all_MEA_biallelic_no_meta_encoded.vcf.gz ## output from MEGANE (removed metadata, encoded genotypes)
    23-2-22_vcf_filtered_evidence_for_ins.txt ## final vcd file with only insertions agent from the reference 
    23-2-27_vcf_ins_del_all_final.txt ## final vcf file that I used for downstream analysis
### TE annotation library
    lib2_renamed_nodups_fixedunkowns_plusgraph_tab.rep

##### File location on the cluster:
## Genome (AstCal1.2)
~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames.fa) 
## paths to all cichlid files to generate the vcf (fofn from Bettina Fischer)
~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/data
## Generated vcf files 
~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/23-2-1_MEGANE_all_cichilds/vcf_for_phasing
## TE annotation libarary (MWCichlidTE-3): 
~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/barbara/Lib3
