library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)
DATA_DIR <- "/Users/chris/Desktop/DSC_291/Final_Project_Data"
FIGURES_DIR <- "/Users/chris/Desktop/DSC_291/final_project/dsc291final/Figures_scATAC"
setwd(DATA_DIR)

combinedSeurat <- readRDS("scRNA_after_annotations.RDS")
combinedSignac <- readRDS("scATAC_after_lsi.RDS")

# Transfer Annotations with combinedSeurat as reference.
variable_genes <- VariableFeatures(combinedSeurat) # we will only consider the variable features found in scRNA-seq for annotations.
gene_activities <- GeneActivity(
  object = combinedSignac,
  features = variable_genes, # Will use all genes in the annotation
  extend.upstream = 2000, # 2kb upstream
  extend.downstream = 0 # Gene body only
)
saveRDS(gene_activities, "scATAC_gene_activity.RDS")
combinedSignac[["RNA"]] <- CreateAssayObject(counts = gene_activities)
DefaultAssay(combinedSignac) <- "RNA"
combinedSignac <- NormalizeData(combinedSignac)
transfer_anchors <- FindTransferAnchors(
  reference = combinedSeurat,
  query = combinedSignac,
  normalization.method = "LogNormalize",
  reference.reduction = "harmony",
  recompute.residuals = FALSE,
  dims = 1:20
)
predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = combinedSeurat$cell_type,
  weight.reduction = combinedSignac[["lsi"]],
  dims = 1:20
)
# Add main prediction and the maximum prediction score
combinedSignac$predicted.celltype <- predictions$predicted.id
combinedSignac$prediction.score <- predictions$prediction.score.max

# Prediction Quality/Summaries/Verifications.
print("Number of cells with predictions:")
print(length(combinedSignac$predicted.celltype))
print("Prediction score summary:")
print(summary(combinedSignac$prediction.score))
print("Distribution of predicted cell types:")
print(table(combinedSignac$predicted.celltype))

# We will remove low quality predictions.
Idents(combinedSignac) <- "predicted.celltype"
combinedSignac <- combinedSignac[, combinedSignac$prediction.score > 0.5]

# UMAP Clustering with Annotations
combinedSignac <- RunUMAP(
  object = combinedSignac,
  dims = 1:20, # Using LSI dimensions as input
  n.components = 20, # Specifying 20 UMAP dimensions
  reduction = "lsi", # Using LSI as input reduction
  reduction.name = "umap"
)
combinedSignac <- FindMultiModalNeighbors(
  object = combinedSignac,
  reduction.list = list("umap", "lsi"),
  dims.list = list(1:20, 2:20),
  modality.weight.name = c("RNA.weight", "peaks.weight"),
  verbose = TRUE
)

# Dimplot
clusters <- DimPlot(subset(combinedSignac, group == "Tumor"), label = TRUE, repel = TRUE, reduction = "umap", dims = c(2, 3))
ggsave(paste0(FIGURES_DIR, "/4.Annotated_with_anchors_clusters.png"),
  plot = clusters,
  width = 200, height = 200, units = "mm"
)

# Linking Peaks
DefaultAssay(combinedSignac) <- "peaks"

library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
peaks <- granges(combinedSignac[["peaks"]])
seqlevels(peaks) <- paste0("chr", c(1:22, "X", "Y", "M"))
seqinfo(peaks) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(peaks)]

combinedSignac <- RegionStats(
  object = combinedSignac,
  regions = peaks,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
peak_links <- LinkPeaks(
  object = combinedSignac,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = NULL, # Set to specific genes if you want to subset
  min.cells = 10, # Minimum number of cells that a peak should be detected in
  distance = 500000 # Maximum distance between peak and gene (500kb)
)
saveRDS(peak_links, "peak_links.RDS")

combinedSignac <- SortIdents(combinedSignac)
