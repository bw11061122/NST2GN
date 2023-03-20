#### 2023-03-11
#### matrix
library(ggplot2)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(reshape)
library(tidyr)
library(tidyverse)

df <- read.table("~/Desktop/23-3-11_families_not_in_all_final.txt", sep="\t", header=T)
df <- df %>% select(!(sum_ins_all))
df <- df %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "family")
df.matrix <- as.matrix(df)
heatmap(df.matrix, Colv = NA, Rowv = NA, xlab="population", ylab="family")

