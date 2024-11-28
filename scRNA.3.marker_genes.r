library(Seurat)
library(ggplot2)
library(data.table)
library(patchwork)
DATA_DIR <- "/Users/chris/Desktop/DSC_291/Final_Project_Data"
setwd(DATA_DIR)
combinedSeurat <- readRDS("scRNA_after_umap.rds")
combinedSeurat <- FindNeighbors(combinedSeurat,
  dims = 1:20,
  reduction = "harmony"
)
# We need to find optimal resolution for the optimal number of clusters.
# Tried 0.5, but really noisy clustering. Not working.
# Will use Calinski-Harabasz
test_resolutions <- seq(0.1, 1.9, by = 0.1)
n_clusters <- numeric(length(test_resolutions))
ch_scores <- numeric(length(test_resolutions))
reduction_matrix <- Embeddings(combinedSeurat, reduction = "harmony")[, 1:20]
for (i in seq_along(test_resolutions)) {
  seurat_obj <- FindClusters(combinedSeurat,
    resolution = test_resolutions[i],
    verbose = FALSE
  )
  n_clusters[i] <- length(unique(Idents(seurat_obj)))
  # Calinski-Harabasz Index
  if (n_clusters[i] > 1) {
    clusters <- Idents(seurat_obj)
    overall_mean <- colMeans(reduction_matrix)
    bc_ss <- 0
    wc_ss <- 0
    for (cluster in unique(clusters)) {
      cluster_points <- reduction_matrix[clusters == cluster, , drop = FALSE]
      n_points <- nrow(cluster_points)
      cluster_mean <- colMeans(cluster_points)
      bc_ss <- bc_ss + n_points * sum((cluster_mean - overall_mean)^2)
      wc_ss <- wc_ss + sum(sweep(cluster_points, 2, cluster_mean)^2)
    }
    ch_scores[i] <- (bc_ss / (n_clusters[i] - 1)) / (wc_ss / (nrow(reduction_matrix) - n_clusters[i]))
  } else {
    ch_scores[i] <- NA
  }
}
# Run using the best resolution.
optimal_idx <- which.max(ch_scores)
optimal_resolution <- test_resolutions[optimal_idx]
print(paste("Total Clusters:", n_clusters[optimal_idx]))
combinedSeurat <- FindClusters(combinedSeurat,
  resolution = optimal_resolution,
  verbose = TRUE
)
combinedSeurat <- JoinLayers(combinedSeurat)
# Check the Dimplot
DimPlot(combinedSeurat,
  reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, label.size = 3, repel = TRUE
)
# Marker Genes
all.markers <- FindAllMarkers(combinedSeurat,
  only.pos = TRUE, min.pct = 0.3,
  logfc.threshold = 0.25
)
write.table(all.markers,
  file = paste0("scrna_all_markers_optRes", optimal_resolution, ".txt"),
  quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE
)
optimal_resolution
saveRDS(combinedSeurat, "scRNA_after_clustering.rds")
