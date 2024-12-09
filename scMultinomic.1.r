library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
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
gene_activities <- readRDS("scATAC_gene_activity.RDS")
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
library(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(combinedSignac) <- "peaks"
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
  genes.use = NULL,
  min.cells = 10,
  distance = 500000
)
saveRDS(peak_links, "peak_links.RDS")
saveRDS(combinedSignac, "scATAC_after_link_peaks.RDS")

# =========================================================================#
# Coverage Plot
deg_list <- fread("scrna_all_markers_optRes0.2.txt", header = TRUE)
deg_genes <- unique(deg_list$gene)
# Get annotations to map gene names to genomic coordinates
annotations <- Annotation(peak_links)
deg_coords <- annotations[annotations$gene_name %in% deg_genes, ]
peak_links <- readRDS("peak_links.RDS")
combinedSignac <- readRDS("scATAC_after_link_peaks.RDS")
DefaultAssay(peak_links) <- "peaks"


find_nearby_peaks <- function(gene_info, peaks) {
  chr <- as.character(seqnames(gene_info))
  start_pos <- as.numeric(unlist(strsplit(GRangesToString(gene_info), "-"))[2])
  end_pos <- as.numeric(unlist(strsplit(GRangesToString(gene_info), "-"))[3])
  # Find peaks on the same chromosome (modified for new format with colon)
  chr_peaks <- peaks[grep(paste0("^", chr, ":"), peaks)]
  # Extract start positions from peak names (modified for new format)
  peak_starts <- as.numeric(gsub(".*?:(\\d+)-\\d+", "\\1", chr_peaks))
  # Find peaks within range of the gene
  nearby_indices <- which(peak_starts >= (start_pos - 8000) & peak_starts <= (end_pos + 5000))
  return(chr_peaks[nearby_indices])
}

# Modify the syntax of peaks
assay <- peak_links[["peaks"]]
new_names <- gsub("^(chr[^-]+)-", "\\1:", rownames(assay@counts))
rownames(assay@counts) <- new_names
rownames(assay@data) <- new_names
rownames(assay@meta.features) <- new_names
if (length(assay@var.features) > 0) {
  assay@var.features <- gsub("^(chr[^-]+)-", "\\1:", assay@var.features)
}
peak_links[["peaks"]] <- assay

create_coverage_plot <- function(peak_links, region, nearby_peaks) {
  links_data <- Links(peak_links)
  gene_name <- region$gene_name
  links_df <- as.data.frame(mcols(links_data))
  gene_links <- links_df[links_df$gene == gene_name, ]
  if (nrow(gene_links) == 0) {
    return(NULL)
  }
  coverage_matrix <- GetAssayData(peak_links, layer = "counts", assay = "peaks")
  peak_coverage <- rowMeans(coverage_matrix[nearby_peaks, , drop = FALSE])
  plot_data <- data.frame()
  for (peak in names(peak_coverage)) {
    # Parse peak coordinates
    peak_coords <- StringToGRanges(gsub(":", "-", peak))
    start_pos <- start(peak_coords)
    end_pos <- end(peak_coords)
    # Add row to plot data
    plot_data <- rbind(plot_data, data.frame(
      start = start_pos,
      end = end_pos,
      coverage = as.numeric(peak_coverage[peak])
    ))
  }
  links_plot_data <- data.frame()
  for (i in seq_len(nrow(gene_links))) {
    peak_coords <- StringToGRanges(gsub(":", "-", gene_links$peak[i]))
    peak_center <- start(peak_coords) + (end(peak_coords) - start(peak_coords)) / 2
    gene_center <- start(region) + (end(region) - start(region)) / 2
    links_plot_data <- rbind(links_plot_data, data.frame(
      x = peak_center,
      xend = gene_center,
      score = gene_links$score[i]
    ))
  }

  gene_model_data <- data.frame(
    start = start(region),
    end = end(region),
    y = -0.5,
    yend = -0.5
  )
  max_coverage <- max(plot_data$coverage, na.rm = TRUE)
  y_max <- max_coverage * 1.2
  p <- ggplot() +
    geom_rect(
      data = plot_data,
      aes(
        xmin = start,
        xmax = end,
        ymin = 0,
        ymax = coverage
      ),
      fill = "blue",
      alpha = 0.5
    ) +
    geom_curve(
      data = links_plot_data,
      aes(
        x = x,
        xend = xend,
        y = max_coverage * 0.5,
        yend = max_coverage * 0.5,
        alpha = score
      ),
      curvature = 0.2,
      color = "red",
      linewidth = 0.5
    ) +
    geom_segment(
      data = gene_model_data,
      aes(
        x = start,
        xend = end,
        y = y,
        yend = yend
      ),
      color = "black",
      linewidth = 2
    ) +
    # Add baseline
    geom_hline(yintercept = 0, color = "black") +
    # Customize theme and labels
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.position = "none"
    ) +
    labs(
      x = paste0(as.character(seqnames(region)), " position"),
      y = "Average Coverage",
      title = paste0("Coverage Plot for ", region$gene_name)
    ) +
    scale_x_continuous(
      limits = c(
        start(region) - 8000,
        end(region) + 5000
      )
    ) +
    scale_y_continuous(
      limits = c(-1, y_max)
    )
  return(p)
}

for (i in 1:length(deg_coords)) {
  gene <- deg_coords$gene_name[i]
  gene_coords <- deg_coords[i]
  # Find nearby peaks
  nearby_peaks <- find_nearby_peaks(gene_coords, rownames(peak_links))
  if (length(nearby_peaks) > 0) {
    p <- create_coverage_plot(
      peak_links = peak_links,
      region = gene_coords,
      nearby_peaks = nearby_peaks
    )
    if (!is.null(p)) {
      ggsave(
        filename = paste0(FIGURES_DIR, "/coverage_plot_", gene, ".png"),
        plot = p,
        width = 12,
        height = 8,
        bg = "white"
      )
    }
  }
}

# Initialize counters
total_genes <- length(deg_coords)
genes_with_links <- 0
for (i in 1:length(deg_coords)) {
  gene <- deg_coords$gene_name[i]
  gene_coords <- deg_coords[i]

  # Get links data
  links_data <- Links(peak_links)
  links_df <- as.data.frame(mcols(links_data))
  gene_links <- links_df[links_df$gene == gene, ]

  # Check if gene has links
  if (nrow(gene_links) > 0) {
    genes_with_links <- genes_with_links + 1
  }
}
# Calculate proportion
proportion <- genes_with_links / total_genes
cat(sprintf("Total genes: %d\n", total_genes))
cat(sprintf("Genes with links: %d\n", genes_with_links))
cat(sprintf("Proportion: %.2f%%\n", proportion * 100))
