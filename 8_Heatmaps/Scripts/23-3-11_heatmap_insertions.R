sum_ins_ap <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ap]
sum_ins_ask <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ask]
sum_ins_asm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_asm]
sum_ins_cch <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cch]
sum_ins_crhoa <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_crhoa]
sum_ins_cvm <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cvs]
sum_ins_cvs <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cvm]
sum_ins_czeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_czeb]
sum_ins_diplo <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_diplo]
sum_ins_fr <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_fr]
sum_ins_fuel <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_fuel]
sum_ins_trew <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_trew]
sum_ins_mzeb <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_mzeb]
sum_ins_oto <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_oto]
sum_ins_rham <- vcf[type=="MEI", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_rham]

vcf_sub_ins <- vcf[type=="MEI",c("ID", "family", "order", "superfamily", "type")]
vcf_dt_sums_ins <- as.data.table(cbind(vcf_sub_ins, sum_ins_ap, sum_ins_ask, sum_ins_asm, 
                                   sum_ins_crhoa, sum_ins_cch, sum_ins_cvm, sum_ins_cvs, sum_ins_czeb, 
                                   sum_ins_diplo, sum_ins_fr, sum_ins_fuel,  sum_ins_trew, sum_ins_mzeb, 
                                   sum_ins_oto, sum_ins_rham))

### this is still too large: i will plot sum of insertions for each family 
vcf_df_fam_ins <- as.data.frame(cbind(vcf[type=="MEI","family", with=FALSE], sum_ins_ap, sum_ins_ask, sum_ins_asm, 
                                  sum_ins_crhoa, sum_ins_cch, sum_ins_cvm, sum_ins_cvs, sum_ins_czeb, 
                                  sum_ins_diplo, sum_ins_fr, sum_ins_fuel,  sum_ins_trew, sum_ins_mzeb, 
                                  sum_ins_oto, sum_ins_rham))
sums <- c("sum_ins_ap", "sum_ins_ask", "sum_ins_asm", 
          "sum_ins_crhoa", "sum_ins_cch", "sum_ins_cvm", "sum_ins_cvs", "sum_ins_czeb", 
          "sum_ins_diplo", "sum_ins_fr", "sum_ins_fuel",  "sum_ins_trew", "sum_ins_mzeb", 
          "sum_ins_oto", "sum_ins_rham")
vcf_df_agg_ins <- aggregate(.~family, vcf_df_fam_ins, sum) ## sum by family 

### generate a dataset where you have summed up all families 
### with superfamily and order annotation
fam_names <- vcf_df_agg_ins$family
df_ins_ann <- vcf %>% filter(type=="MEI") %>% select(c(order, superfamily, family)) %>% filter(family %in% fam_names) 
df_ins_ann <- unique(df_ins_ann) 
vcf_df_ins_merge <- merge(df_ins_ann, as.data.table(vcf_df_agg_ins), by="family")
vcf_df_ins_merge <-vcf_df_ins_merge[order(order, superfamily)]
vcf_df_num_ins <- vcf_df_ins_merge %>% select(sums) %>% as.matrix()
colnames(vcf_df_num_ins) <- c("Alticorpus peterdaviesi.Cape_Maclear",
                              "Astatotilapia calliptera.Lake_Kingiri", "Astatotilapia calliptera.Lake_Masoko", "Chilotilapia rhoadesii",
                              "Copadichromis chrysonotus.Lake_Malombe", "Copadichromis virginalis.Lake_Malombe", "Copadichromis virginalis.Southwest_arm", 
                              "Cynotilapia zebroides.Cape_Maclear", "Diplotaxodon limnothrissa.Southwest_arm", "Fossorochromis rostratus.Lake_Malombe", 
                              "Labeotropheus fuelleborni.Chilumba", "Labeotropheus trewavasae.Chilumba", "Maylandia zebra.Cape_Maclear", 
                              "Otopharynx argyrosoma.Southeast_arm", "Rhamphochromis longiceps")                        
rownames(vcf_df_num_ins) <- vcf_df_ins_merge$family
pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-3-10_heatmap_try_sum.pdf", width=30, height=40)
heatmap(vcf_df_num_ins, Rowv=NA)
dev.off()

## column annotation 
col_ann <- meta %>% select(clade, pop_id)
col_ann <- unique(col_ann)
col_ann <- col_ann %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "pop_id")

## row annotation 
row_ann <- vcf_df_ins_merge %>% select(order, superfamily, family) %>% column_to_rownames(var="family")
row_ann$ann <- paste(row_ann$order, row_ann$superfamily, sep="_")
row_ann1 <- row_ann %>% select(ann)
colnames(row_ann1) <- "order_superfamily"
row_ann_order <- row_ann %>% select(order)
row_ann_order <- droplevels(row_ann_order)
annotation_colors = list(
  clade = c(AstCal="#AED581", Mbuna="#5E35B1", Utaka="#1B5E20", Benthic="#B71C1C", 
            Deep="#1A237E", Rhampho="#EF6C00", Diplo="#827717"),
  order = c(DNA="#EF5350", LINE="#2196F3", LTR="#43A047", SINE="#9575CD",
            RC="#D1C4E9"))

## modify the stupid column names function in pheatmap to get rotated labels
## courtesy of https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps 
library(grid) 
draw_colnames_90 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                  ns=asNamespace("pheatmap"))
## heatmap 
pheatmap(vcf_df_num_ins, scale="row",annotation_row = row_ann_order,
         annotation_col = col_ann,
         annotation_colors = annotation_colors,
         main="Insertions absent from the reference", cluster_rows=F, 
         cellwidth = 30, 
         cellheight = 0.8,
         cluster_cols=T, show_rownames = F,
         show_annotation_name = F,
         cexCol=8, fontsize = 20)


