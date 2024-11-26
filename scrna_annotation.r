library(Seurat)
library(ggplot2)
library(data.table)
library(patchwork)
DATA_DIR <- "/Users/chris/Desktop/DSC_291/Final_Project_Data"
setwd(DATA_DIR)

combinedSeurat <- readRDS("scRNA_after_umap.rds")
combinedSeurat <- FindNeighbors(combinedSeurat, dims = 20)
combinedSeurat <- FindClusters(combinedSeurat, resolution = 0.5)
combinedSeurat <- JoinLayers(combinedSeurat)
all.markers <- FindAllMarkers(combinedSeurat,
  only.pos = TRUE, min.pct = 0.3,
  logfc.threshold = 0.25
)
write.table(all.markers,
  file = "scrna_all_markers.txt", quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE
)
