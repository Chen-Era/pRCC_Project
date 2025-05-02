# Main execution script for scRNA-seq data analysis

library(Seurat)
library(dplyr)
library(ggplot2)

# Execute data import and preprocessing
source("import.R")

# Execute quality control
source("qc.R")

# Execute data normalization and dimensionality reduction
source("normalization.R")

# Execute cell clustering and marker gene identification
source("clustering.R")

# Execute cell-cell communication analysis
source("cell_comm.R")

# Execute transcription factor regulatory network analysis
source("tf_network.R")

# Execute differential gene expression analysis
source("diff_expr.R")

# Execute gene set and pathway enrichment analysis
source("pathway.R") 