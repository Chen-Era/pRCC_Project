# pRCC_Project

## System Requirements

- Operating system: macOS or Linux recommended.
- R: 4.x (tested with recent versions).
- Python: 3.8+ for pySCENIC.
- Memory/CPU: SCENIC and fastMNN can be compute-intensive; use sufficient RAM and threads.

## Dependencies (per script)

General recommendations: install Bioconductor packages via BiocManager; Python environment for pySCENIC.

- Seruat3Normalization.R: Seurat, dplyr, ggplot2 (Seurat ecosystem).
- Seurat3Cluster.R: SingleCellExperiment, scater, plyr, reshape2, Seurat, mclust, dplyr.
- MNN.R: batchelor, SingleCellExperiment, scater, Seurat (plus plotting dependencies).
- ClusterProfile-GO.R: clusterProfiler, org.Hs.eg.db, dplyr (Bioconductor).
- ClusterProfile-Pathway.R: clusterProfiler, org.Hs.eg.db, dplyr (Bioconductor).
- GSVA.ssGSEA.R: GSVA, limma, stringr.
- Monocle2.Pseudotime.R: monocle, Seurat, reshape2, dplyr, ggplot2.
- inferCNV.R: ggplot2.
- QuSAGE.R: qusage, Seurat.
- SCENIC.sh: pySCENIC (arboreto, ctx, etc.); requires organism-specific motif/cistarget databases.

R installation examples (selected):

- install.packages(c("Seurat", "dplyr", "ggplot2", "plyr", "reshape2", "mclust", "limma", "stringr"))
- if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
- BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "GSVA", "qusage", "monocle", "SingleCellExperiment", "scater", "batchelor"))

Python/pySCENIC:

- pip install pyscenic

