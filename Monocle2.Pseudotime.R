library(monocle)
library(Seurat)
library(reshape2)
library(dplyr)

# CLI: <seuset.rds> <output_dir> [gene_file]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript Monocle2.Pseudotime.R <seuset.rds> <output_dir> [gene_file]\n")
  quit(save = "no", status = 1)
}

rdsFile <- args[1]
outpath <- args[2]
geneFile <- if (length(args) >= 3) args[3] else NULL
if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

seuset = readRDS(rdsFile)

seuset <- AddMetaData(object = seuset, metadata = seuset@active.ident, col.name = "Cluster")
seuCDSAll = importCDS(seuset, import_all = TRUE)
if(!is.null(geneFile) && file.exists(geneFile)){
    genelist = unique(read.table(geneFile, sep = "\t", header = TRUE)[,1])
}else{
    genelist = VariableFeatures(seuset)
}

seuCDSAll <- setOrderingFilter(seuCDSAll, genelist)
seuCDSAll <- estimateSizeFactors(seuCDSAll)
seuCDS <- reduceDimension(seuCDSAll, max_components = 2, method = "DDRTree")
seuCDS <- orderCells(seuCDS)

p1 <- plot_cell_trajectory(seuCDS, color_by = "State")
p2 <- plot_cell_trajectory(seuCDS, color_by = "Cluster")
p3 <- plot_cell_trajectory(seuCDS, color_by = "Pseudotime")

ggplot2::ggsave(file.path(outpath, "trajectory_state.png"), p1, width = 8, height = 6)
ggplot2::ggsave(file.path(outpath, "trajectory_cluster.png"), p2, width = 8, height = 6)
ggplot2::ggsave(file.path(outpath, "trajectory_pseudotime.png"), p3, width = 8, height = 6)
cat(paste0("Monocle2 pseudotime plots saved to ", outpath, "\n"))

