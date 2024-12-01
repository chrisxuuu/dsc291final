# Load libraries and set directory
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(data.table)

DATA_DIR <- "/Users/chris/Desktop/DSC_291/Final_Project_Data"
setwd(DATA_DIR)

# Read data
counts <- Matrix::readMM("./scATAC_data/GSE227265_BinarizedPeaks.mtx.gz")
meta_data <- fread("./scATAC_data/GSE227265_TFMotifActivity.csv.gz", header = TRUE)

# Get cell barcodes from metadata
meta_barcodes <- colnames(meta_data)
meta_barcodes <- meta_barcodes[2:length(meta_barcodes)]
colnames(counts) <- meta_barcodes
writeLines(meta_barcodes, "valid_barcodes.txt")

rownames(counts) <- meta_barcodes
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  fragments = "./scATAC_data/GSE227265_fragments_AllSamples.tsv.gz",
  min.cells = 10,
  min.features = 200
)

# Create Seurat object
atac_seurat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)
