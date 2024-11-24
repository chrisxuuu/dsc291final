#### Libraries
library(Seurat)
library(ggplot2)
library(harmony)
library(data.table)
library(patchwork)
DATA_DIR <- "/Users/chris/Desktop/DSC_291/Final_Project_Data"
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
ggsave("1.QC_plots.png", width = 230, units = "mm")

#### Clustering
