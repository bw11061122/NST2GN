sum_del_ap <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ap]
sum_del_ask <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_ask]
sum_del_asm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_asm]
sum_del_cch <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cch]
sum_del_crhoa <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_crhoa]
sum_del_cvm <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cvs]
sum_del_cvs <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_cvm]
sum_del_czeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_czeb]
sum_del_diplo <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_diplo]
sum_del_fr <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_fr]
sum_del_fuel <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_fuel]
sum_del_trew <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_trew]
sum_del_mzeb <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_mzeb]
sum_del_oto <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_oto]
sum_del_rham <- vcf[type=="MEA", apply(.SD, 1, function(x) length(which(x=="1"))+ 2*length(which(x=="2"))), .SDcols=meta_rham]

vcf_sub_del <- vcf[type=="MEA",c("ID", "family", "order", "superfamily", "type")]
vcf_dt_sums_del <- as.data.table(cbind(vcf_sub_del, sum_del_ap, sum_del_ask, sum_del_asm, 
                                       sum_del_crhoa, sum_del_cch, sum_del_cvm, sum_del_cvs, sum_del_czeb, 
                                       sum_del_diplo, sum_del_fr, sum_del_fuel,  sum_del_trew, sum_del_mzeb, 
                                       sum_del_oto, sum_del_rham))

### this is still too large: i will plot sum of delertions for each family 
vcf_df_fam_del <- as.data.frame(cbind(vcf[type=="MEA","family", with=FALSE], sum_del_ap, sum_del_ask, sum_del_asm, 
                                      sum_del_crhoa, sum_del_cch, sum_del_cvm, sum_del_cvs, sum_del_czeb, 
                                      sum_del_diplo, sum_del_fr, sum_del_fuel,  sum_del_trew, sum_del_mzeb, 
                                      sum_del_oto, sum_del_rham))
sums <- c("sum_del_ap", "sum_del_ask", "sum_del_asm", 
          "sum_del_crhoa", "sum_del_cch", "sum_del_cvm", "sum_del_cvs", "sum_del_czeb", 
          "sum_del_diplo", "sum_del_fr", "sum_del_fuel",  "sum_del_trew", "sum_del_mzeb", 
          "sum_del_oto", "sum_del_rham")
vcf_df_agg_del <- aggregate(.~family, vcf_df_fam_del, sum) ## sum by family 

### generate a dataset where you have summed up all families 
### with superfamily and order annotation
fam_names_del <- vcf_df_agg_del$family
df_del_ann <- vcf %>% filter(type=="MEA") %>% select(c(order, superfamily, family)) %>% filter(family %in% fam_names_del) 
df_del_ann <- unique(df_del_ann) 
vcf_df_del_merge <- merge(df_del_ann, as.data.table(vcf_df_agg_del), by="family")
vcf_df_del_merge$combo <- paste(vcf_df_del_merge$order, vcf_df_del_merge$family, sep="_")
### will need to remove the incorrect matchings
c_rem <- c("DNA_Rex-Babar-6", "DNA_Rex-Babar-13", 
           "DNA_Rex-Babar-16", "DNA_Rex-Babar-12",
           "DNA_Dong-R4-1", "DNA_Dong-R4-2",
           "DNA_L2-13", "DNA_L2-17", "DNA_L2-6",
           "DNA_Ngaro-1", "DNA_Gypsy-35",
           "DNA_Pao-36", "DNA_L1-3", "DNA_L1-18",
           "DNA_Unknown-12", "DNA_Unknown-63",
           "DNA_Unknown-24", "DNA_L2-32", "DNA_Rex-Babar-17")
vcf_df_del_merge <- vcf_df_del_merge %>% filter(!(combo %in% c_rem))
vcf_df_del_merge <-vcf_df_del_merge[order(order, superfamily)]
vcf_df_num_del <- vcf_df_del_merge %>% select(sums) %>% as.matrix()
colnames(vcf_df_num_del) <- c("A. peterdaviesi",
                                   "A. calliptera (LK)", "A. calliptera (LM)", "C. rhoadesii",
                                   "C. chrysonotus", "C. virginalis (LMB)", "C. virginalis (SWA)", 
                                   "C. zebroides", "D. limnothrissa", "F. rostratus", 
                                   "L. fuelleborni", "L. trewavasae", "M. zebra", 
                                   "O. argyrosoma", "R. longiceps") ## pheatmap package
rownames(vcf_df_num_del) <- vcf_df_del_merge$family

## column annotation 
col_ann <- meta %>% select(clade, pop_id)
col_ann <- unique(col_ann)
col_ann <- col_ann %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "pop_id")

## row annotation 
row_ann_del <- vcf_df_del_merge %>% select(order, family) 
row_ann_order_del <- row_ann_del %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "family")
row_ann_order_del <- droplevels(row_ann_order_del)
annotation_colors = list(
  clade = c(AstCal="#AED581", Mbuna="#5E35B1", Utaka="#1B5E20", Benthic="#B71C1C", 
            Deep="#1A237E", Rhampho="#EF6C00", Diplo="#827717"),
  order = c(DNA="#EF5350", LINE="#2196F3", LTR="#43A047", SINE="#9575CD",
            RC="#D1C4E9"))

## modify the column names function in pheatmap to get rotated labels
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
pheatmap(vcf_df_num_del, scale="row",annotation_row = row_ann_order_del,
         annotation_col = col_ann,
         annotation_colors = annotation_colors,
         main="Comparison of TE families across populations", cluster_rows=F, 
         cellwidth = 30, 
         cellheight = 0.8,
         cluster_cols=T, show_rownames = F,
         show_annotation_name = T,
         cexCol=8, fontsize = 20)

vcf_dt_sums_del <- as.data.table(cbind(vcf_sub_del, sum_del_ap, sum_del_ask, sum_del_asm, 
                                       sum_del_crhoa, sum_del_cch, sum_del_cvm, sum_del_cvs, sum_del_czeb, 
                                       sum_del_diplo, sum_del_fr, sum_del_fuel,  sum_del_trew, sum_del_mzeb, 
                                       sum_del_oto, sum_del_rham))

### this is still too large: i will plot sum of delertions for each family 
vcf_df_sfam_del <- as.data.frame(cbind(vcf[type=="MEA","superfamily", with=FALSE], sum_del_ap, sum_del_ask, sum_del_asm, 
                                      sum_del_crhoa, sum_del_cch, sum_del_cvm, sum_del_cvs, sum_del_czeb, 
                                      sum_del_diplo, sum_del_fr, sum_del_fuel,  sum_del_trew, sum_del_mzeb, 
                                      sum_del_oto, sum_del_rham))
sums <- c("sum_del_ap", "sum_del_ask", "sum_del_asm", 
          "sum_del_crhoa", "sum_del_cch", "sum_del_cvm", "sum_del_cvs", "sum_del_czeb", 
          "sum_del_diplo", "sum_del_fr", "sum_del_fuel",  "sum_del_trew", "sum_del_mzeb", 
          "sum_del_oto", "sum_del_rham")
vcf_df_agg_sfam_del <- aggregate(.~superfamily, vcf_df_sfam_del, sum) ## sum by family 
sfam_names <- vcf_df_agg_sfam_del$superfamily
df_del_ann_sfam <- vcf %>% filter(type=="MEI") %>% select(c(order, superfamily)) %>% filter(superfamily %in% sfam_names) 
df_del_ann_sfam  <- unique(df_del_ann_sfam) 
vcf_df_del_sfam_merge <- merge(df_del_ann_sfam, as.data.table(vcf_df_agg_sfam_del), by="superfamily")
vcf_df_del_sfam_merge <-vcf_df_del_sfam_merge[order(order, superfamily)]
vcf_df_num_sfam_del <- vcf_df_del_sfam_merge %>% select(sums) 
vcf_df_num_sfam_del <- vcf_df_num_sfam_del%>% as.matrix()
rownames(vcf_df_num_sfam_del) <- vcf_df_del_sfam_merge$superfamily
colnames(vcf_df_num_sfam_del) <- c("A. peterdaviesi",
                                   "A. calliptera (LK)", "A. calliptera (LM)", "C. rhoadesii",
                                   "C. chrysonotus", "C. virginalis (LMB)", "C. virginalis (SWA)", 
                                   "C. zebroides", "D. limnothrissa", "F. rostratus", 
                                   "L. fuelleborni", "L. trewavasae", "M. zebra", 
                                   "O. argyrosoma", "R. longiceps") ## pheatmap package
col_ann <- meta %>% select(clade, pop_id)
col_ann <- unique(col_ann)
col_ann <- col_ann %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "pop_id")

## row annotation 
row_ann_sfam <- vcf_df_del_sfam_merge %>% select(order, superfamily) %>% column_to_rownames(var="superfamily")

## set colors 
annotation_colors = list(
  clade = c(AstCal="#AED581", Mbuna="#5E35B1", Utaka="#1B5E20", Benthic="#B71C1C", 
            Deep="#1A237E", Rhampho="#EF6C00", Diplo="#827717"),
  order = c(DNA="#EF5350", LINE="#2196F3", LTR="#43A047", SINE="#9575CD",
            RC="#D1C4E9"))

## modify the column names function in pheatmap to get rotated labels
## courtesy of https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps 
library(grid) 
draw_colnames_90 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                  ns=asNamespace("pheatmap"))
## heat map 
pheatmap(vcf_df_num_sfam_del, scale="row", annotation_row = row_ann_sfam,
         annotation_col = col_ann,
         cellwidth=25,
         cellheight=20,
         annotation_colors = annotation_colors,
         main="Comparison of TE superfamilies across populations", cluster_rows=F, 
         cluster_cols=T, show_rownames = T,
         fontsize=22, cexCol=8)
