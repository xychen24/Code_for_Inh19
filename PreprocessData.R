#' @description WKS1KO+Inh19 vs Inh19

########################
###
# 1.Initialize----
###
########################
## 1.1 load packages:
########################
library(Seurat)
library(ggpubr)
library(DoubletFinder)
#######################
## 1.2 setting dic:
#######################
setwd("...")

#######################
###
# 2.Handle data seperately----
###
#######################
## 2.1 input data & prefilter
#######################
s07.data <- Read10X(data.dir = "../04.Matrix", gene.column=1)
colnames(s07.data) <- paste0("s07_", colnames(s07.data))
s07 <- CreateSeuratObject(counts = s07.data, project = "s07", min.cells = 3, min.features = 200)
dim(s07)


s08.data <- Read10X(data.dir = "../04.Matrix", gene.column=1)
colnames(s08.data) <- paste0("s08_", colnames(s08.data))
s08 <- CreateSeuratObject(counts = s08.data, project = "s08", min.cells = 3, min.features = 200)
dim(s08)

#######################
## 2.2 filter
#######################
# Mito gene
s07[["percent.mt"]] <- PercentageFeatureSet(s07, pattern = "^MT-")
s08[["percent.mt"]] <- PercentageFeatureSet(s08, pattern = "^MT-")


# filter
s07.seu <- subset(s07, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 9000) 
s08.seu <- subset(s08, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 9000) 
dim(s07.seu)
dim(s08.seu)


#######################
## 2.3 SCTransform
#######################
rsl <- seq(0.5, 1.5, by = 0.1)

# s07
s07.seu <- SCTransform(s07.seu, vst.flavor = "v2", vars.to.regress = c("percent.mt")) %>% 
  RunPCA(npcs = 50) %>%
  FindNeighbors(dims = 1:40) %>%
  FindClusters(resolution = rsl, verbose = FALSE) %>%
  RunUMAP(dims = 1:40)

# s08
s08.seu <- SCTransform(s08.seu, vst.flavor = "v2", vars.to.regress = c("percent.mt")) %>% 
  RunPCA(npcs = 50) %>%
  FindNeighbors(dims = 1:40) %>%
  FindClusters(resolution = rsl, verbose = FALSE) %>%
  RunUMAP(dims = 1:40)


#######################
## 2.4 remove doublet
#######################
# s07
bcmvn_s07 <- paramSweep_v3(s07.seu, PCs = 1:40, sct = TRUE) %>%
  summarizeSweep(GT = FALSE) %>%
  find.pK()
pK_bcmvn <- bcmvn_s07$pK[which.max(bcmvn_s07$BCmetric)] %>% as.character() %>% as.numeric()
annotations <- s07.seu@meta.data$SCT_snn_res.0.5
homotypic.prop <- modelHomotypic(annotations)

doublet.rate <- 0.0375
nExp_poi <- round(doublet.rate*nrow(s07.seu@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
s07.seu1 <- doubletFinder_v3(s07.seu, PCs = 1:40, pN = 0.25,
                             pK = pK_bcmvn,
                             nExp = nExp_poi.adj,
                             reuse.pANN = FALSE,
                             sct = TRUE)
colnames(s07.seu1[[]])
table(s07.seu1$DF.classifications_0.25_0.29_174)
s07.seu2 <- subset(s07.seu1, subset = DF.classifications_0.25_0.29_174 == "Singlet")

# s08
bcmvn_s08 <- paramSweep_v3(s08.seu, PCs = 1:40, sct = TRUE) %>%
  summarizeSweep(GT = FALSE) %>%
  find.pK()
pK_bcmvn <- bcmvn_s08$pK[which.max(bcmvn_s08$BCmetric)] %>% as.character() %>% as.numeric()
annotations <- s08.seu@meta.data$SCT_snn_res.0.5
homotypic.prop <- modelHomotypic(annotations)
doublet.rate <- 0.0237
nExp_poi <- round(doublet.rate*nrow(s08.seu@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
s08.seu1 <- doubletFinder_v3(s08.seu, PCs = 1:40, pN = 0.25,
                             pK = pK_bcmvn,
                             nExp = nExp_poi.adj,
                             reuse.pANN = FALSE,
                             sct = TRUE)
colnames(s08.seu1[[]])
table(s08.seu1$DF.classifications_0.25_0.01_68)
s08.seu2 <- subset(s08.seu1, subset = DF.classifications_0.25_0.01_68 == "Singlet")

#######################
## 2.5 Remove resolution
#######################
rsl <- seq(0.5, 1.5, by = 0.1)

index <- match(paste0("SCT_snn_res.", rsl), colnames(s07.seu2@meta.data))
s07.seu2@meta.data <- s07.seu2@meta.data[,-index]
colnames(s07.seu2[[]])

index <- match(paste0("SCT_snn_res.", rsl), colnames(s08.seu2@meta.data))
s08.seu2@meta.data <- s08.seu2@meta.data[,-index]
colnames(s08.seu2[[]])
