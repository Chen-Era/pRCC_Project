# pRCC_Project

A pipeline for Integrated Single-Cell and Spatial Transcriptomics Reveals Tertiary Lymphoid Structures in High-Grade Papillary Renal Cell Carcinoma with Inferior Vena Cava Tumor Thrombus data analysis, from raw data processing to downstream analysis.

## Structure

- `exec.R`: Main execution script
- `import.R`: Data import and preprocessing
- `qc.R`: Quality control
- `normalization.R`: Data normalization and dimensionality reduction
- `clustering.R`: Cell clustering and marker gene identification
- `cell_comm.R`: Cell-cell communication analysis
- `tf_network.R`: Transcription factor regulatory network analysis
- `diff_expr.R`: Differential gene expression analysis
- `pathway.R`: Gene set and pathway enrichment analysis

## Dependencies

- Seurat (v3.1.4+)
- dplyr
- ggplot2
- CellPhoneDB
- pySCENIC
- Monocle3
- QuSAGE
- clusterProfiler
