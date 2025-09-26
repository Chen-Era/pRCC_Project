# pRCC_Project

An integrative script collection for single-cell and spatial transcriptomics analysis in papillary renal cell carcinoma (pRCC). The repository covers end-to-end steps from raw matrix import, quality control, normalization/dimensionality reduction, clustering and marker discovery, batch integration, functional enrichment, pathway analysis, CNV metrics, pseudotime inference, transcription factor regulatory network (SCENIC), and gene set activity analysis (GSVA/QuSAGE).

This README provides dependencies, CLI usage, minimal examples, I/O conventions, and recommended practices to facilitate reproducible research and manuscript supplementation.

## Scope and Modules

- Seruat3Normalization.R: Import, filtering, cell cycle scoring, normalization, variable feature selection, scaling, PCA/UMAP/t-SNE; saves a processed Seurat object (RDS).
- Seurat3Cluster.R: Seurat-based graph clustering, cluster statistics, marker gene discovery; outputs summary tables and RDS.
- MNN.R: Batch integration with batchelor::fastMNN, UMAP/t-SNE, clustering, and results saving.
- ClusterProfile-GO.R: Gene Ontology enrichment (BP/MF/CC) via clusterProfiler; outputs CSV.
- ClusterProfile-Pathway.R: KEGG pathway enrichment via clusterProfiler; outputs CSV.
- GSVA.ssGSEA.R: Single-sample GSEA (ssGSEA) using GSVA; outputs the activity matrix.
- Monocle2.Pseudotime.R: Monocle2-based pseudotime trajectory from a Seurat RDS; generates trajectory plots.
- inferCNV.R: CNV-related metrics and visualization from provided matrices and group labels; predicts tumor cells and outputs summaries and plots.
- QuSAGE.R: QuSAGE gene set activity per cluster from a Seurat RDS; outputs density curves, confidence interval plots and tables.
- SCENIC.sh: pySCENIC workflow wrapper (GRN inference, motif enrichment, AUCell scoring).

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

## CLI Usage and Minimal Examples

Unless otherwise stated, scripts create the output directory if it does not exist. Results are written to the specified output_dir or to ./Result within the working directory.

1) Normalization and DR (Seruat3Normalization.R)
- Usage: Rscript Seruat3Normalization.R <MatrixDir> <sample.txt> <output_dir> <prefix> <genecolumn> <cellMinGene> <cellMaxGene> <geneExpMinCell> <MaxMTPercent> <regressOut>
- Notes:
  - MatrixDir: path to 10x matrix or a sparse matrix directory
  - sample.txt: sample metadata file
  - prefix: output prefix
  - genecolumn: gene column name
  - other parameters: filtering and regression options
- Output: processed Seurat RDS with PCA/UMAP/t-SNE.

2) Clustering and Markers (Seurat3Cluster.R)
- Usage: Rscript Seurat3Cluster.R <input_rds> <prefix> <output_dir> <resolution> <marker_method> <logFC_threshold> <min.pct> <only.pos>
- Output:
  - <prefix>_GraphClust.Statistics.txt
  - <prefix>_GraphClust.Summary_Cell.txt
  - <prefix>_GraphClust.AllMarkerGenes.txt
  - <prefix>_GraphClust.seuset.rds

3) Batch Integration (MNN.R)
- Usage: Rscript MNN.R <input_rds> <output_dir> <minCell> <k> <resolution>
- Output: integrated object, UMAP/t-SNE embeddings, and clustering results (filenames follow script outputs).

4) GO Enrichment (ClusterProfile-GO.R)
- Usage: Rscript ClusterProfile-GO.R <gene_list.txt> [output_dir]
- Output: GO_enrichment.csv (includes BP/MF/CC).

5) KEGG Pathway (ClusterProfile-Pathway.R)
- Usage: Rscript ClusterProfile-Pathway.R <gene_list.txt> [output_dir] [organism:hsa]
- Output: KEGG_enrichment.csv.

6) ssGSEA (GSVA.ssGSEA.R)
- Usage: Rscript GSVA.ssGSEA.R <gene_set.gmt|txt> <expression_matrix.txt> [output_dir]
- Output: ssgsea.txt (rows: samples/cells; columns: gene sets).

7) Pseudotime (Monocle2.Pseudotime.R)
- Usage: Rscript Monocle2.Pseudotime.R <seuset.rds> <output_dir> [gene_file]
- Output: trajectory_state.png, trajectory_cluster.png, trajectory_pseudotime.png.

8) CNV Metrics (inferCNV.R)
- Usage: Rscript inferCNV.R <observations.txt> <references.txt> <group.txt> <output_dir>
- Output:
  - CNV_cells.txt, CNV_clusters.txt
  - correlation_against_meanSquare.pdf, CNVPlot.png

9) QuSAGE Gene Set Activity (QuSAGE.R)
- Usage: Rscript QuSAGE.R <input_rds> <geneset_files_comma_separated> <output_dir>
- Notes: gene set files may be GMT or two-column txt (type, gene); multiple files comma-separated.
- Output: under <output_dir>/QuSAGE/<geneset_name>/ per cluster: table_*.txt, plotDC_*.png/pdf, plotCI_*.png.

10) SCENIC (SCENIC.sh, pySCENIC based)
- Typical steps: GRN inference (grn), motif enrichment (ctx), AUCell scoring (aucell).
- Recommendation: run `bash SCENIC.sh -h` for full parameters and database configuration; prepare organism-specific feather rankings and motif annotations.
- Output: SCENIC intermediate and final results in the specified directory (candidate networks, enrichment outputs, AUCell scores).

## I/O Conventions

- Gene lists: one-column SYMBOL text; GO/KEGG scripts internally map to ENTREZID (human: org.Hs.eg.db).
- Expression matrices: tab-delimited text; orientation depends on script (GSVA/ssGSEA expects columns as samples/cells).
- Seurat objects: RDS format; use this repository’s normalization script to ensure slot compatibility.
- CNV matrices: space-delimited text; observations: tumor candidates; references: reference cells; group.txt: two columns (cellID, group/cluster).

## Recommended Workflow (Example)

1) Normalization/DR → 2) Clustering/Markers → 3) Batch Integration (if needed) → 4) Functional/Pathway enrichment → 5) GSVA/QuSAGE → 6) Trajectory/SCENIC → 7) CNV metrics integration and visualization.

Modules are composable; adapt the order as needed for your study design.

## Reproducibility and Best Practices

- Dependency installation: some scripts only check for missing packages and do not auto-install; install per the instructions above if errors occur.
- Version compatibility: Seurat v3/v4 slot names differ; QuSAGE/Monocle2 are handled with basic compatibility but prefer objects produced by these scripts.
- Performance: SCENIC and MNN are compute-intensive; validate on small subsets first. Increase RAM and threads for full runs.
- Paths and encoding: quote paths with spaces; ensure UTF-8 text without BOM.
- Reproducibility: set seeds externally, record session info, and maintain output directories by stage (./Result/ for outputs; ./Test/ for small-scale validation before full runs).

## Notes for Manuscript Supplement

- Report software and package versions used (R, Python, Seurat, clusterProfiler, GSVA, qusage, monocle, batchelor, pySCENIC).
- Provide exact CLI commands and parameter values used in the study (see examples above).
- Document input data sources and preprocessing (e.g., QC thresholds: cellMinGene/cellMaxGene/MaxMTPercent; marker discovery method and thresholds; integration k and resolution settings).
- Include a schematic workflow if desired (Normalization → Clustering → Integration → Enrichment → Activity → Trajectory → CNV → SCENIC).
- File naming convention: keep the default filenames produced by scripts to ease cross-referencing in supplementary tables and figures.

## Additional Remarks

- The script name "Seruat3Normalization.R" reflects the file present in the repository (typo in name preserved for compatibility). If you prefer, rename it to "Seurat3Normalization.R" and update references accordingly.
- For environment locking, consider R renv or Conda (environment.yml) to pin versions and ensure reproducibility.

## Repository

GitHub: https://github.com/Chen-Era/pRCC_Project
