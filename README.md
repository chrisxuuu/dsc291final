# DSC 291 Final Project

Chris Xu, Hanchang Cai

Edit: Dec. 8, 2024

## Introduction:

Integrated analysis of scRNA-seq and scATAC-seq data for HCC patients. The overlaps of Differentially Accessible Regions (DARs) and Differentially Expressed Genes (DEGs) will further extend the interpretations from traditional scRNA-seq DEGs by adding 
chromatin accessibility information. This can reveal cell-type specific enhancer-gene relationships, which can lead to novel therapeutic targets. 

Inspired by: Bai, Y., Deng, X., Chen, D., Han, S., Lin, Z., Li, Z., Tong, W., Li, J., Wang, T., Liu, X., Liu, Z., Cui, Z., & Zhang, Y. (2024). Integrative Analysis Based on ATAC-Seq and RNA-seq reveals a novel oncogene PRPF3 in hepatocellular carcinoma. Clinical Epigenetics, 16(1). https://doi.org/10.1186/s13148-024-01769-w.

#### Software: 

* R version 4.4.1 (2024-06-14)
  * Signac v1.14.0
  * Seurat v5.1.0
  * GenomicRanges v1.56.2
  * BSgenome v1.72.0 (BSgenome.Hsapiens.UCSC.hg38 v1.4.5)
* MACS2 v2.2.9.1 (Gaspar, J. M. (2018). Improved peak-calling with MACS2. https://doi.org/10.1101/496521)

#### Installation Notes:

* All R packages can be found on CRAN or BioConductor.
* MACS2 was installed using pip 24.3.1 on Python 3.10.

#### Datasets Used:

* scRNA-seq: GSE156625 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156625)
* scATAC-seq: GSE227265 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227265)
* Panglao DB Marker Gene Annotation (https://www.panglaodb.se/markers.html)

#### SOP:

1. Download all Supplementary Files from the scRNA-seq and scATAC-seq Accessions and Panglao DB Marker Gene Annotations.
2. Start with scRNA-seq. Rename 10x files to standard 10x names. Seperate into two folders for 
   HCC and HCCF.
3. Run scRNA scripts in order of the script numbers.
4. Run scATAC scripts in order of the script numbers.
5. Run scMultinomic script.

#### Repo Link: 

https://github.com/chrisxuuu/dsc291fina
