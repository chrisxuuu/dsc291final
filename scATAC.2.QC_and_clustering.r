library(Signac)
library(Seurat)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(dtplyr)
library(patchwork)
library(AnnotationHub)
DATA_DIR <- "/Users/chris/Desktop/DSC_291/Final_Project_Data"
FIGURES_DIR <- "/Users/chris/Desktop/DSC_291/final_project/dsc291final/Figures_scATAC"
setwd(DATA_DIR)

# Count Matrix
counts <- readRDS("computed_peaks_count_mtx.RDS")

# Meta Data (Idents)
tumor_barcodes <- fread(
  file = "./scATAC_data/GSE227265_TFMotifActivityMalignantCells.csv.gz",
  header = TRUE,
  nrows = 1
) %>%
  colnames() %>%
  gsub("\\.", "-", .)
tumor_barcodes <- tumor_barcodes[-1]
adjacent_barcodes <- fread(
  file = "./scATAC_data/GSE227265_TFMotifActivity.csv.gz",
  header = TRUE,
  nrows = 1
) %>%
  colnames()
adjacent_barcodes <- adjacent_barcodes[-1]
metadata <- data.table(
  barcode = adjacent_barcodes,
  group = ifelse(adjacent_barcodes %in% tumor_barcodes,
    "Tumor", "Adjacent"
  )
)
table(metadata$group)

# Make the Signac Object
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "./scATAC_data/GSE227265_fragments_AllSamples.tsv.gz",
  min.cells = 10,
  min.features = 200
)
# Ignore Warnings, we extended each region by +/- 150 bp after MACS2.
combinedSignac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
# Check a summary of the object.
combinedSignac[["peaks"]] # We have around 810502 features.
peaks.keep <- seqnames(granges(combinedSignac)) %in% standardChromosomes(granges(combinedSignac))
combinedSignac <- combinedSignac[as.vector(peaks.keep), ]

# Add Gene Annotations
# We will use EnsDb annotation for human GRCh38 since data is on the
# Illumina NovaSeq 6000 platform (GPL24676) in 2023.
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("EnsDb", "Homo Sapiens", "GRCh38"))
annotations <- ahDb[[length(ahDb)]]
annotations <- genes(annotations) # We will keep both coding and non-coding for annotations.
annotations <- annotations[seqnames(annotations) != "LRG_432"] # Remove LRG_432 since it is causing errors.
seqlevelsStyle(annotations) <- "UCSC"
standard_chroms <- paste0("chr", c(1:22, "X", "Y"))
annotations <- keepSeqlevels(annotations, standard_chroms, pruning.mode = "coarse")
Annotation(combinedSignac) <- annotations

# QC
# Nucleosome
combinedSignac <- NucleosomeSignal(object = combinedSignac)
# TSS
combinedSignac <- TSSEnrichment(object = combinedSignac)
# FRIP
fragment_path <- "./scATAC_data/GSE227265_fragments_AllSamples.tsv.gz"
stats <- CountFragments(fragment_path)
counts_dt <- as.data.table(stats)
counts_dt <- counts_dt[CB %in% combinedSignac$barcode]
combinedSignac$pct_reads_in_peaks <-
  counts_dt[match(combinedSignac$barcode, CB)]$frequency_count /
    counts_dt[match(combinedSignac$barcode, CB)]$reads_count * 100
saveRDS(combinedSignac, "scATAC_QC_Stats.RDS")

# Scatter Plot
combinedSignac <- readRDS("scATAC_QC_Stats.RDS")
scatter <- DensityScatter(combinedSignac, x = "nCount_peaks", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)
# From the scatterplot, we set TSS > 3.5 filter.
combinedSignac$nucleosome_group <- ifelse(combinedSignac$nucleosome_signal > 4, "NS > 4", "NS < 4")
nucleo_hist <- FragmentHistogram(object = combinedSignac, group.by = "nucleosome_group")
ggsave(paste0(FIGURES_DIR, "/1.QC_plots.png"),
  plot = scatter + nucleo_hist,
  height = 130, width = 130 * 2, units = "mm"
)
violin <- VlnPlot(
  object = combinedSignac,
  features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks"),
  group.by = "group",
  pt.size = 0.1,
  ncol = 4
)
ggsave(paste0(FIGURES_DIR, "/2.QC_violin.png"),
  plot = violin,
  height = 130, width = 130 * 2, units = "mm"
)
# FILTERS: nucleosome_signal < 4
#          FRIP > 0.15
#          TSS enrichment score > 3.5
combinedSignac <- subset(
  x = combinedSignac,
  subset = nCount_peaks > 1000 & nCount_peaks < 20000 &
    nucleosome_signal < 4 & nucleosome_signal > 0.5 &
    TSS.enrichment > 3.5
)

# Clustering
combinedSignac <- RunTFIDF(combinedSignac)
combinedSignac <- FindTopFeatures(combinedSignac, min.cutoff = "q0")
combinedSignac <- RunSVD(combinedSignac)
combinedSignac <- RunUMAP(object = combinedSignac, reduction = "lsi", dims = 2:30) # We skip the first dim
combinedSignac <- FindNeighbors(object = combinedSignac, reduction = "lsi", dims = 2:30)
combinedSignac <- FindClusters(object = combinedSignac, verbose = FALSE, algorithm = 3)
dimplot <- DimPlot(object = combinedSignac, label = TRUE) + NoLegend() +
  ggtitle("LSI Clusters")
ggsave(paste0(FIGURES_DIR, "/3.Clustering_dimplot.png"), plot = dimplot, width = 200, height = 200, units = "mm")
saveRDS(combinedSignac, "scATAC_after_lsi.RDS")
