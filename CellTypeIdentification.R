#' @description WKS1KO+Inh19 vs Inh19

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
library(dplyr)
library(export)
#######################
## 1.2 setting dic:
#######################
setwd("...")
#######################
###
# 2.Visualize
###
#######################
# integrated_snn_res.1
DefaultAssay(ko.combined)
Idents(ko.combined) <- ko.combined$integrated_snn_res.1
DimPlot(ko.combined, reduction = "umap", 
        group.by = "integrated_snn_res.1",
        label = TRUE) + NoLegend()

#######################
###
# 3. CellTypeIdentification
###
#######################
## 3.1 FindAllMarkers in 16 clusters
#######################
DefaultAssay(ko.combined) <- "SCT"
Idents(ko.combined) <- ko.combined$integrated_snn_res.1
ko.combined <- PrepSCTFindMarkers(ko.combined)
clusters.all.marker <- FindAllMarkers(ko.combined, assay = "SCT", slot = "data",
                                      test.use = "MAST")
write.csv(clusters.all.marker, file = "clusters.all.marker_MAST.csv")
save(clusters.all.marker, file = "clusters.all.marker.RData")

clusters.sig.marker <- clusters.all.marker[which(clusters.all.marker$p_val_adj<0.05),]
clusters.sig.marker.top <- clusters.sig.marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(ko.combined, features = unique(clusters.sig.marker.top$gene), size = 2) + NoLegend()
write.csv(table(ko.combined$integrated_snn_res.1), file = "cellnumber_16clusters.csv")

#######################
## 3.2 RenameIdents
#######################
Idents(ko.combined) <- ko.combined$integrated_snn_res.1
ko.combined <- RenameIdents(ko.combined,
                            "0" = "β cells",
                            "1" = "β cells",
                            "2" = "α cells",
                            "3" = "α cells",
                            "4" = "β cells",
                            "5" = "Pancreatic progenitors",
                            "6" = "β cells",
                            "7" = "α cells",
                            "8" = "δ cells",
                            "9" = "EC cells",
                            "10" = "Proliferation cells",
                            "11" = "EC cells",
                            "12" = "α cells",
                            "13" = "Proliferation cells",
                            "14" = "ε cells",
                            "15" = "Pancreatic progenitors")

#save(ko.combined, file = "ko.combined_after_rename.RData")
#load(file = "ko.combined_after_rename.RData")

ko.combined$cell.type <- Idents(ko.combined)
color_change <- c("β cells"='#70AD47',
                  "α cells"='#ED7D31',
                  "ε cells"="#DC143C",
                  "δ cells"="#264478",
                  "Pancreatic progenitors"='#9E480E',
                  #"Polyhormonal cells"='#5B9BD5',
                  "Proliferation cells" = '#A5A5A5',
                  "EC cells"='#4472C4'
                  # "α-like cells"='#FFC000',
)
ko.combined$cell.type <- factor(ko.combined$cell.type, 
                                levels = c("β cells","α cells", "δ cells", "ε cells", "EC cells",
                                           "Pancreatic progenitors","Proliferation cells"))
prop.table(table(ko.combined$cell.type,ko.combined$orig.ident),margin = 2)
VlnPlot(ko.combined, features = c("INS","GCG"))


ab.combined <- ko.combined[,ko.combined$cell.type %in% c("α cells","β cells")]
prop.table(table(ab.combined$integrated_snn_res.1,ab.combined$orig.ident),margin = 2)
VlnPlot(ab.combined, features = c("INS","GCG"),group.by = "integrated_snn_res.1")


DimPlot(ab.combined, 
        reduction = "umap",
        group.by = "cell.type",
        label = FALSE,
        cols = color_change,
        pt.size = 0.9)
#######################
## 3.3 FindAllMarkers in celltypes
#######################
Idents(ko.combined) <- ko.combined$cell.type
cell.type.marker <- FindAllMarkers(ko.combined, assay = "SCT", slot = "data",
                                      test.use = "MAST")
head(cell.type.marker)
cell.type.marker.sig <- cell.type.marker[which(cell.type.marker$p_val_adj < 0.05),]
top.marker.sig <- cell.type.marker.sig %>% group_by(cell.type.marker.sig$cluster) %>% top_n(10, wt = avg_log2FC)


