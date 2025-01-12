#' @description: WKS1KO+Inh19 vs Inh19
#######################
###
# 1.Initialize----
###
#######################
## 1.1 load packages:
#######################
library(Seurat)
library(ggpubr)
library(patchwork)
library(harmony)
#######################
## 1.2 setting dic:
#######################
setwd("...")
#######################
###
# 2. RunHarmony
###
#######################
DefaultAssay(ko.combined) <- "integrated"
ko.combined <- RunPCA(ko.combined, npcs = 50)
ko.combined@reductions

ko.combined <- RunHarmony(ko.combined, group.by.vars = "orig.ident", assay.use = "integrated",plot_convergence = TRUE)

# 提取信息
hm.em <- Embeddings(hg.combined, 'harmony')
hm.em[1:5, 1:5]

#######################
###
# 3. Pipeline
###
#######################
rsl <- seq(0.2, 1.2, by = 0.1)
PCs <- 40
ko.combined <- RunUMAP(ko.combined, reduction = "harmony", dims = 1:PCs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:PCs) %>%
  FindClusters(resolution = rsl, verbose = FALSE) 


