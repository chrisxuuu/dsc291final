# Load libraries and set directory
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(data.table)
library(Matrix)
DATA_DIR <- "/Users/chris/Desktop/DSC_291/Final_Project_Data"
setwd(DATA_DIR)

# Peaks Processing (paper uses MACS2 to call)
mem.maxVSize(vsize = Inf) # Note: Requires a bit over 32GB ram.
frags <- fread("./scATAC_data/GSE227265_fragments_AllSamples.tsv.gz",
  col.names = c("chr", "start", "end", "barcode", "count")
)
frag_ranges <- GRanges(
  seqnames = frags$chr,
  ranges = IRanges(start = frags$start, end = frags$end),
  barcode = frags$barcode
)
system("macs2 callpeak -f BED -g hs --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -q 0.05 -t ./scATAC_data/GSE227265_fragments_AllSamples.tsv.gz -n initial_peaks")
peaks <- fread("initial_peaks_summits.bed")
peaks[, c("V2", "V3") := list(V2 - 150, V3 + 150)] # From the paper, each peak was extended by 150 for each side.
peaks_gr <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))

# Remove ENCODE Blacklist Peaks Problematic Regions.
# From https://www.encodeproject.org/files/ENCFF356LFX/
blacklist <- fread("./scATAC_data/ENCFF356LFX.bed.gz")
blacklist <- GRanges(
  seqnames = blacklist$V1,
  ranges = IRanges(start = blacklist$V2, end = blacklist$V3)
)
peaks_gr <- peaks_gr[!overlapsAny(peaks_gr, blacklist)]
peaks_gr$score <- countOverlaps(peaks_gr, frag_ranges)

# We will insert the peaks into the count matrix (based on coverage)
counts <- Matrix::readMM("./scATAC_data/GSE227265_BinarizedPeaks.mtx.gz")
n_peaks <- nrow(counts)
peaks_gr <- peaks_gr[order(-peaks_gr$score)][1:n_peaks]
peak_names <- paste0(seqnames(peaks_gr), ":", start(peaks_gr), "-", end(peaks_gr))
writeLines(peak_names, "computed_peak_names.txt")

# We will rebuild the count matrix
computed_peaks <- Matrix(0,
  nrow = length(peaks_gr),
  ncol = ncol(counts),
  sparse = TRUE
)
rownames(computed_peaks) <- peak_names
# For column names, we extract from metadata for orginal barcodes.
# Get cell barcodes from metadata
meta_data <- fread("/Users/chris/Desktop/DSC_291/Final_Project_Data/scATAC_data/GSE227265_TFMotifActivity.csv.gz", header = TRUE)
meta_barcodes <- colnames(meta_data)
meta_barcodes <- meta_barcodes[2:length(meta_barcodes)]
colnames(counts) <- meta_barcodes
# Find Counts while matching barcodes to original matrix.
overlaps <- findOverlaps(peaks_gr, frag_ranges)
peak_indices <- queryHits(overlaps)
cell_indices <- match(
  frag_ranges$barcode[subjectHits(overlaps)],
  colnames(counts)
)
valid_matches <- !is.na(cell_indices)
peak_indices <- peak_indices[valid_matches]
cell_indices <- cell_indices[valid_matches]
unique_coords <- unique(data.table(i = peak_indices, j = cell_indices))
computed_peaks <- Matrix::sparseMatrix(
  i = unique_coords$i,
  j = unique_coords$j,
  x = 1,
  dims = dim(computed_peaks)
)
colnames(computed_peaks) <- meta_barcodes
rownames(computed_peaks) <- peak_names
saveRDS(computed_peaks, "computed_peaks_count_mtx.RDS") # Note: Rows are unsorted.

# Too Slow, use the vectorized solution above.
if (FALSE) {
  for (i in seq_along(overlaps)) {
    peak_idx <- queryHits(overlaps[i])
    cell_idx <- match(
      frag_ranges$barcode[subjectHits(overlaps[i])],
      colnames(counts)
    )
    if (!is.na(cell_idx)) {
      computed_peaks[peak_idx, cell_idx] <- 1
    }
  }
}
