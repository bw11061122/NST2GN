#### [2023-2-27]
##### Trying to maybe do this smarter
##### Repeat sums but now do this no the merged file

# trait	insID	dist2Gene	VEP	geneID	geneSymbol
# colour	chr7:41166580-41166583_O	0	intron	ENSACLG00000014805	ntrk3a #### did not find this - may be more complicated because it has an inversion 
# colour	chr3:8249800-8249807_O	0	intron	ENSACLG00000018516	nlgn2a #### is in the dataset ### hAT_AC-17
# Caudal peduncle depth	chr22:13852670-13852671_O	0	intron	ENSACLG00000013131	trps1 #### unknown 19


### identify the insertions which correspond to the sites you got from Bettina
awk '$2 == "chr7" ' 23-2-1_MEGANE_all_cichlids/vcf_for_phasing/cichlid_all_MEI_biallelic_no_meta_encoded.vcf \
> 23-2-1_MEGANE_all_cichlids/vcf_for_phasing/23-3-6_vcf_for_phasing_chr7.vcf ## get just for chromosome 7
grep '411665' 23-2-1_MEGANE_all_cichlids/vcf_for_phasing/23-3-6_vcf_for_phasing_chr7.vcf > \
23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-6_ntrk3a_hits.txt

awk '$2 == "chr22" ' 23-2-1_MEGANE_all_cichlids/vcf_for_phasing/cichlid_all_MEI_biallelic_no_meta_encoded.vcf \
> 23-2-1_MEGANE_all_cichlids/vcf_for_phasing/23-3-6_vcf_for_phasing_chr22.vcf ## get just for chromosome 7
grep '138526' 23-2-1_MEGANE_all_cichlids/vcf_for_phasing/23-3-6_vcf_for_phasing_chr22.vcf > \
23-2-1_MEGANE_all_cichlids/23-2-16_subset/23-3-6_ntrk3a_hits.txt