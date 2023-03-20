## 23-2-16
# install.packages("umap")
# install.packages("tsne")
library(umap)
library(data.table)
library(tsne)
library(datasets)
library(plotly)

## UMAP
### tutorial here: https://cran.r-project.org/web/packages/umap/vignettes/umap.html 
### define the UMAP function (from the tutorial)
plot.iris <- function(x, labels,
                      main="A UMAP visualization",
                      colors=c("pink", "red", "grey", "green", "orange", "violet", "blue", "darkblue", 
                               "palegreen1", "royalblue1", "palevioletred", "paleturquoise", "royalblue4",
                               "red4", "rosybrown", "darkorange3", "goldenrod", "forestgreen", "violetred", 
                               "tomato4"),
                      pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
                      cex.main=1, cex.legend=0.85) {
  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  } 
  xylim <- range(layout)
  xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u <- unique(labels)
  legend.pos <- "topleft"
  legend.text <- as.character(labels.u)
  if (add) {
    legend.pos <- "bottomleft"
    legend.text <- paste(as.character(labels.u), legend.suffix)
  }
  
  legend(legend.pos, legend=legend.text, inset=0.03,          
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

#### read the data (rows: cichlids, columns: insertions)

## read the metadata
meta <- fread("~/Desktop/Genetics_project/23-1-23/Data/meta/23-2-16_cichlid_meta_standardized_dataset.csv",
                              sep=",", header = T) ## this is the full metadat
meta <- meta[, clade:=as.factor(clade)] 
meta <- meta[, sex:=as.factor(sex)] 
meta <- meta[, location:=as.factor(location)] 
meta <- meta[, sublocation:=as.factor(sublocation)] 
meta <- meta[, species:=as.factor(species)] 
meta <- meta[, genus:=as.factor(genus)] 
meta <- meta[, names:=paste(genus, species)] 
meta <- meta[, name_loc:=paste0(names,".",location)] 
meta <- meta[, names:=as.factor(names)] 
meta <- meta[, name_loc:=as.factor(name_loc)] 
meta <- meta[, population_id := population][(genus %in% c("Rhamphochromis", "Chilotilapia")) | is.na(population), id := names][]

#### subset meta for specific clades
meta_ap <- meta[name_loc=="Alticorpus peterdaviesi.Cape_Maclear",]
meta_as_kingiri <- meta[name_loc=="Astatotilapia calliptera.Lake_Kingiri",]
meta_as_masoko <- meta[name_loc=="Astatotilapia calliptera.Lake_Masoko",]
meta_crhoa <- meta[names=="Chilotilapia rhoadesii",]
meta_cch <- meta[name_loc=="Copadichromis chrysonotus.Lake_Malombe",]
meta_cv_malombe <- meta[name_loc=="Copadichromis virginalis.Lake_Malombe",]
meta_cv_swarm <- meta[name_loc=="Copadichromis virginalis.Southwest_arm",]
meta_czeb <- meta[name_loc=="Cynotilapia zebroides.Cape_Maclear",]
meta_dip <- meta[name_loc=="Diplotaxodon limnothrissa.Southwest_arm",]
meta_fr <- meta[name_loc=="Fossorochromis rostratus.Lake_Malombe",]
meta_fuel <- meta[name_loc=="Labeotropheus fuelleborni.Chilumba",]
meta_trew <- meta[name_loc=="Labeotropheus trewavasae.Chilumba",]
meta_mzeb <- meta[name_loc=="Maylandia zebra.Cape_Maclear",]
meta_oto <- meta[name_loc=="Otopharynx argyrosoma.Southeast_arm",]
meta_rham <- meta[names=="Rhamphochromis longiceps",]

#### get the ids of cichlids in each clade that you will analyse 
meta_ap_names <- unlist(meta_ap[,primary_id]) 
meta_as_kingiri_names <- unlist(meta_as_kingiri[,primary_id]) 
meta_as_masoko_names <- unlist(meta_as_masoko[,primary_id]) 
meta_crhoa_names <- unlist(meta_crhoa[,primary_id]) 
meta_cch_names <- unlist(meta_cch[,primary_id]) 
meta_cv_malombe_names <- unlist(meta_cv_malombe[,primary_id]) 
meta_cv_swarm_names <- unlist(meta_cv_swarm[,primary_id]) 
meta_czeb_names <- unlist(meta_czeb[,primary_id]) 
meta_dip_names <- unlist(meta_dip[,primary_id]) 
meta_fr_names <- unlist(meta_fr[,primary_id])
meta_fuel_names <- unlist(meta_fuel[,primary_id]) 
meta_trew_names <- unlist(meta_trew[,primary_id]) 
meta_mzeb_names <- unlist(meta_mzeb[,primary_id]) 
meta_oto_names <- unlist(meta_oto[,primary_id]) 
meta_rham_names <- unlist(meta_rham[,primary_id]) 

### umap 
## don't do on all columns now 
vcf_sub_data <- vcf_sub[, c(1:644951)]
vcf_sub_data <- as.data.table(lapply(vcf_sub_data, as.numeric))
vcf.umap100 <- umap(vcf_sub_data, n_neighbors=100) # I am wondering if 12 is correct given that I have 12 individuals / population
vcf.umap50 <- umap(vcf_sub_data, n_neighbors=50) # I am wondering if 12 is correct given that I have 12 individuals / population
vcf.umap10 <- umap(vcf_sub_data, n_neighbors=10) # I am wondering if 12 is correct given that I have 12 individuals / population
vcf.umap5 <- umap(vcf_sub_data, n_neighbors=5) 
ids_of_interest <- unlist(vcf_sub[,primary_id])
meta_of_int <- meta_all[primary_id %in% ids_of_interest, ]
vcf.labels <- meta_of_int[,clade]
plot.iris(vcf.umap100, vcf.labels)
plot.iris(vcf.umap50, vcf.labels)
plot.iris(vcf.umap10, vcf.labels)
plot.iris(vcf.umap5, vcf.labels)

iris.data <- iris[, grep("Sepal|Petal", colnames(iris))]
iris.labels <- iris[, "Species"]
iris.umap <- umap(iris.data)
plot.iris(iris.umap, iris.labels)

vcf_sub_data <- vcf_sub[, c(1:644951)]
vcf_sub_data <- as.data.table(lapply(vcf_sub_data, as.numeric))
vcf.umap100 <- umap(vcf_sub_data, n_neighbors=100) # I am wondering if 12 is correct given that I have 12 individuals / population
vcf.umap50 <- umap(vcf_sub_data, n_neighbors=50) # I am wondering if 12 is correct given that I have 12 individuals / population
vcf.umap10 <- umap(vcf_sub_data, n_neighbors=10) # I am wondering if 12 is correct given that I have 12 individuals / population
vcf.umap5 <- umap(vcf_sub_data, n_neighbors=5) 
ids_of_interest <- unlist(vcf_sub[,primary_id])
vcf.labels <- meta_all[primary_id %in% ids_of_interest, name_loc]
plot.iris(vcf.umap100, vcf.labels)
plot.iris(vcf.umap50, vcf.labels)
plot.iris(vcf.umap10, vcf.labels)
plot.iris(vcf.umap5, vcf.labels)

vcf_mock_data <- vcf_sub[,c(1:100)]
vcf_mock_data <- as.data.table(lapply(vcf_sub_data, as.numeric))
vcf.umap100_mock <- umap(vcf_mock_data, n_neighbors=100) # I am wondering if 12 is correct given that I have 12 individuals / population
vcf.umap50_mock <- umap(vcf_mock_data, n_neighbors=50) # I am wondering if 12 is correct given that I have 12 individuals / population
vcf.umap10_mock <- umap(vcf_mock_data, n_neighbors=10) # I am wondering if 12 is correct given that I have 12 individuals / population
vcf.umap5_mock <- umap(vcf_mock_data, n_neighbors=5) 
vcf.umap_mock <- umap(vcf_mock_data) 

plot.iris(vcf.umap100_mock, vcf.labels)
plot.iris(vcf.umap50_mock, vcf.labels)
plot.iris(vcf.umap10_mock, vcf.labels)
plot.iris(vcf.umap5_mock, vcf.labels)

plot.iris(vcf.umap_mock, vcf.labels)

vcf.umap.mock <- umap(vcf_mock_data, n_neighbors=20)
plot.iris(vcf.umap.mock, vcf.labels)

########### UMAP new tutorial here: https://datavizpyr.com/how-to-make-umap-plot-in-r/
library(tidyverse)
library(palmerpenguins)
#install.packages("umap")
#install.packages("palmerpenguins")
library(umap)
theme_set(theme_bw(18))

#### New UMAP
vcf <- fread("~/Desktop/Genetics_project/23-2-20/23-2-22_vcf_filtered_evidence_for_ins.txt",
            header=T, sep="\t")
### your meta is in meta_of_int
vcf_num <- vcf[,7:186]
set.seed(1)
vcf_num_sample <- vcf_num[sample(.N, 10000)]
vcf_trans <- as.data.frame(t(vcf_num_sample))
vcf_trans$id <- c(1:180)
meta_of_int[,id := 1:180]
umap <- vcf_trans %>%
  select(where(is.numeric)) %>%
  scale() %>% 
  umap()
umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(id=row_number())%>%
  inner_join(meta_of_int, by="id")
umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = name_loc)) +
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_10k_seed1.pdf", height=16, width=30)

set.seed(5)

vcf_num_sample <- vcf_num[sample(.N, 40000)]
vcf_trans <- as.data.frame(t(vcf_num_sample))
vcf_trans$id <- c(1:180)
meta[,id := 1:180]
######################################################################
### trying different neighbour numbers for example
library(RColorBrewer)

nn100_umap <- vcf_trans %>%
  select(where(is.numeric)) %>%
  umap(n_neighbors=100)
umap_df100 <- nn100_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(id=row_number())%>%
  inner_join(meta, by="id")
umap_df100 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = population_id)) +
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_40k_seed5_nn100.pdf", height=10, width=14)

##### check number of neighbours 
set.seed(5)
nn12_umap <- vcf_trans %>%
  select(where(is.numeric)) %>%
  umap(n_neighbors=12)
umap_df12 <- nn12_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(id=row_number())%>%
  inner_join(meta, by="id")
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colourCount = length(unique(meta[,population_id]))
umap_df12 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color=population_id))+
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  scale_color_manual(values = getPalette(colourCount))+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_40k_seed5_nn12.pdf", height=10, width=14)

nn12_umap <- vcf_trans %>%
  select(where(is.numeric)) %>%
  umap(n_neighbors=12)
umap_df12 <- nn12_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(id=row_number())%>%
  inner_join(meta, by="id")
umap_df12 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = population_id)) +
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_40k_seed5_nn12_ns.pdf", height=10, width=14)

nn5_umap <- vcf_trans %>%
  select(where(is.numeric)) %>%
  umap(n_neighbors=5)
umap_df5 <- nn5_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(id=row_number())%>%
  inner_join(meta, by="id")
umap_df5 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = population_id)) +
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_40k_seed5_nn5.pdf", height=10, width=14)

umap_df5 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = clade)) +
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_10k_seed5_nn5_clade.pdf", height=10, width=14)

umap_df12 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = clade)) +
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_40k_seed5_nn12_clade.pdf", height=10, width=14)

umap_df100 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = clade)) +
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_10k_seed5_nn100_clade.pdf", height=10, width=14)

##### maybe n-neighbours 20?
nn30_umap <- vcf_trans %>%
  select(where(is.numeric)) %>%
  scale() %>% 
  umap(n_neighbors=30)
umap_df30 <- nn30_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(id=row_number())%>%
  inner_join(meta_of_int, by="id")
umap_df30 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = name_loc)) +
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_10k_seed5_nn30.pdf", height=10, width=14)
umap_df30 %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = clade)) +
  geom_point(size=3, alpha=0.5)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme_bw(base_size=20)+
  theme(legend.position="right")
ggsave("~/Desktop/UMAP_test_10k_seed5_nn30_clade.pdf", height=10, width=14)

#############################################################
set.seed(1)
vcf_sub <- vcf_num[sample(.N, 40000)] # I think 100,000 sounds reasonable
vcf_trans <- transpose(vcf_sub)

vcf.tsne <- tsne(vcf_trans)
tsne_plot_clade <- vcf.tsne %>% 
  as.data.frame() %>% 
  mutate(type=as.factor(cichlid_clades$clade)) %>% 
  ggplot(aes(x=V1,y=V2,color=type))+
  geom_point()+
  theme_bw()+
  labs(title="tSNE insertions",x="tSNE dimension 1", y = "tSNE dimenstion 2")

pdf("~/Desktop/23-2-22_tSNE_by_clade.pdf")
tsne_plot_clade
dev.off()

tsne_plot_pop <- vcf.tsne2 %>% 
  as.data.frame() %>% 
  mutate(type=as.factor(cichlid_clades$population)) %>% 
  ggplot(aes(x=V1,y=V2,color=type))+
  geom_point()+
  theme_bw()+
  labs(title="tSNE insertions",x="tSNE dimension 1", y = "tSNE dimenstion 2")

pdf("23-2-1_MEGANE_all_cichlids/23-2-16_subset/results_PCA/23-2-22_tSNE_by_pop.pdf")
tsne_plot_pop
dev.off()
