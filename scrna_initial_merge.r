#### Quality Check
#### Libraries
library(Seurat)
library(ggplot2)
library(harmony)
library(data.table)
DATA_DIR <- "/Users/chris/Desktop/DSC_291/Final_Project_Data"
PERCENT_MT <- 5
nFEATURE_RNA <- 500
nCOUNT_RNA.min <- 1000
nCOUNT_RNA.max <- 20000
setwd(DATA_DIR)

#### Load scRNA Data
group_data <- lapply(c("HCC", "HCCF"), function(group) {
  sc_data <- Read10X(data.dir = paste0("./scRNA_data/", group))
  sc_data <- CreateSeuratObject(counts = sc_data, project = group)
  sc_data <- subset(
    x = sc_data,
    subset = nFeature_RNA > nFEATURE_RNA &
      nCount_RNA > nCOUNT_RNA.min &
      nCount_RNA < nCOUNT_RNA.max
  )
  sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")
  sc_data <- subset(sc_data, subset = percent.mt < PERCENT_MT)
  sc_data
})
combinedSeurat <- merge(
  x = group_data[[1]], y = group_data[[2]],
  add.cell.ids = c("HCC", "HCCF")
)
saveRDS(combinedSeurat, "scRNA_filter_data.rds")
