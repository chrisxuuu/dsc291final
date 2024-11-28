#### Libraries
library(Seurat)
library(ggplot2)
library(harmony)
library(data.table)
library(patchwork)
DATA_DIR <- "/Users/chris/Desktop/DSC_291/Final_Project_Data"
FIGURES_DIR <- "/Users/chris/Desktop/DSC_291/final_project/dsc291final/Figures_scRNA"
PERCENT_MT <- 5
nFEATURE_RNA <- 500
nCOUNT_RNA.min <- 1000
nCOUNT_RNA.max <- 20000
setwd(DATA_DIR)

#### Violin Plot
combinedSeurat <- readRDS("scRNA_filter_data.rds")
table(combinedSeurat$orig.ident)
violin <- VlnPlot(combinedSeurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3, pt.size = 0, group.by = "orig.ident"
)

#### Feature Scatter
fscatter <- FeatureScatter(combinedSeurat,
  group.by = "orig.ident", feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)
violin / fscatter
ggsave(paste0(FIGURES_DIR, "/1.QC_plots.png"), width = 230, units = "mm")

#### Clustering and Harmony Batch Correction
combinedSeurat <- NormalizeData(combinedSeurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
combinedSeurat <- FindVariableFeatures(combinedSeurat,
  selection.method = "vst",
  nfeatures = 2000
)
combinedSeurat <- ScaleData(combinedSeurat,
  features = VariableFeatures(combinedSeurat)
)
# Prepare PCA for Harmony
combinedSeurat <- RunPCA(combinedSeurat,
  features = VariableFeatures(object = combinedSeurat),
  verbose = FALSE, npcs = 20
)
combinedSeurat <- FindVariableFeatures(combinedSeurat,
  selection.method = "vst", nfeatures = 2000
)
var_feat_plot <- VariableFeaturePlot(combinedSeurat) +
  theme_bw()
top10 <- head(x = VariableFeatures(combinedSeurat), 10)
var_feat_plot <- LabelPoints(
  plot = var_feat_plot, points = top10,
  repel = TRUE, xnudge = 0, ynudge = 0
)
ggsave(paste0(FIGURES_DIR, "/2.Variable_Features_plot.png"), width = 200, height = 150, units = "mm")

# Harmony
combinedSeurat <- RunHarmony(combinedSeurat, "orig.ident", theta = 3)

# For some reason Quick-TRANSfer did not converge.
# We check PCA plot before and after harmony batch correction.
before_pca <- DimPlot(combinedSeurat,
  reduction = "pca", group.by = "orig.ident",
  pt.size = 0.1
) +
  ggtitle("Before Harmony") +
  theme_minimal()
after_pca <- DimPlot(combinedSeurat,
  reduction = "harmony", group.by = "orig.ident",
  pt.size = 0.1
) +
  ggtitle("After Harmony") +
  theme_minimal()
pca_plots <- before_pca + after_pca
ggsave(paste0(FIGURES_DIR, "/3.PCA_plots_harmony.png"), width = 400, height = 150, units = "mm")
# There might not be batch effect. Tested theta = 3 for stronger correction.
# Elbow Plot
elbow_plot <- ElbowPlot(combinedSeurat)
ggsave(paste0(FIGURES_DIR, "/4.Elbow_plot.png"), width = 200, height = 150, units = "mm")
# UMAP
combinedSeurat <- RunUMAP(combinedSeurat, dims = 1:20, reduction = "harmony")
DimPlot(combinedSeurat, reduction = "umap", split.by = "orig.ident")
ggsave(paste0(FIGURES_DIR, "/5.UMAP_plot.png"), width = 400, height = 150, units = "mm")
saveRDS(combinedSeurat, file = "scRNA_after_umap.rds")
