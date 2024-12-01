# DSC 291 Final Project
Hanchang Cai, Chris Xu

Edit: Nov. 30, 2024

## Proposal
Integrated analysis of scRNA-seq and scATAC-seq data for HCC patients. The overlaps of Differentially Accessible Regions (DARs) and Differentially Expressed Genes (DEGs) will further extend the interpretations from traditional scRNA-seq DEGs by adding 
chromatin accessibility information. This can reveal cell-type specific enhancer-gene relationships, which can lead to novel therapeutic targets. 

Inspired by https://pmc.ncbi.nlm.nih.gov/articles/PMC11539654/.

#### Hypothesis: 
We hypothesize that the DARs from scATAC-seq and DEGs from scRNA-seq will have overlaps.

#### Data:
scRNA-seq: GSE156625 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156625)

scATAC-seq: ~~GSE227265 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227265)~~ 

GSE188289 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188289)

#### Methods:
* scRNA-Seq Processing:
  * Seurat Package.
  * QC cell and expressions using feature counts and cell counts, Harmony batch correction.
  * Cell-type clustering UMAP, t-SNE, to classify cells.
  * Auto and Manual cell-type annotations clusters.
  * ~~Select a few cell types to focus on (depending on the number and cell-types from the scATAC-seq data)~~
  * Find DEG's between fetal and primary tumor.
* scATAC-seq Processing:
  * Signac Package.
  * QC using strength of the nucleosome signal per cell, TSS enrichment.
  * Cell-type clustering UMAP, annotations.
  * ~~Select the cell-types matching the selection from scRNA-seq Processing.~~
  * We will find anchor points between scRNA-seq and scATAC-seq using transcriptional activity.
    * Find transcriptional activity of each gene in the scATAC-seq dataset using GeneActivity()
    * This generates a score for each gene's activity.
    * Anchor points between the two datasets are based on dimentionality reduction results using FindTransferAnchors().
      1. Perform dimentionality reduction on both datasets.
      2. Find pairs of cells that are within a certain neighbourhood of each other. These become anchors.
      3. Filter the anchors, remove low confidence ones (based on the number of neighbourhoods away (k.filter))
        * Note: The "size" of each neighbourhood is controlled by max.features parameter.
      4. The rest of the anchors are assigned a score. Network map using the anchor and the cells.
  * Transfer cell-type annotations from scRNA-seq to scATAC-seq using anchors.
  * Call the peaks, then link peaks to gene.
* ~~Correlation analysis of expression matricies for DARs and DEGs from scRNA-Seq and scATAC-seq to identify enhancers.~~

Repo Link: https://github.com/chrisxuuu/dsc291final
