options(stringAsFactors=FALSE)

# CLI: <inputfile.rds> <prefix> <output_dir> [resolution=0.8] [marker_method=wilcox] [logfc=0.25] [minpct=0.1]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript Seurat3Cluster.R <inputfile.rds> <prefix> <output_dir> [resolution=0.8] [marker_method=wilcox] [logfc=0.25] [minpct=0.1] [onlypos=TRUE]\n")
  quit(save = "no", status = 1)
}

rdsFile <- args[1]
prefix <- args[2]
outpath <- if (length(args) >= 3) args[3] else "./Result"
resolution <- if (length(args) >= 4) as.numeric(args[4]) else 0.8
MarkerGeneMethod <- if (length(args) >= 5) args[5] else "wilcox"
threshold <- if (length(args) >= 6) as.numeric(args[6]) else 0.25
minpct <- if (length(args) >= 7) as.numeric(args[7]) else 0.1
OnlyPos <- if (length(args) >= 8) as.logical(args[8]) else TRUE
if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

library(SingleCellExperiment)
library(scater)
library(plyr)
library(reshape2)
library(Seurat)
library(mclust)
library(dplyr)

print("Start")
seuset = readRDS(rdsFile)
assay = DefaultAssay(seuset)

dir.create(file.path(outpath, "GraphClust"), showWarnings = FALSE, recursive = TRUE)

if(!is.null(seuset@commands$RunUMAP)){
    dims = seuset@commands$RunUMAP@params$dims
}else if(!is.null(seuset@commands$RunTSNE)){
    dims = seuset@commands$RunTSNE@params$dims
}
seuset <- FindNeighbors(seuset, reduction = "pca", dims = dims)
seuset <- FindClusters(seuset, resolution = resolution)

table = table(Idents(seuset), seuset@meta.data$orig.ident)
print(table)
summary_table = dcast(data.frame(table), Var1 ~ Var2)
colnames(summary_table)[1] = "Cluster"
write.table(summary_table, file = file.path(outpath, paste0(prefix, "_GraphClust.Statistics.txt")), sep = "\t", row.names = FALSE, quote = FALSE)

data = data.frame(Cell = colnames(seuset), Cluster = Idents(seuset))
write.table(data, file = file.path(outpath, paste0(prefix, "_GraphClust.Summary_Cell.txt")), sep = "\t", row.names = FALSE, quote = FALSE)

allmarkers <- FindAllMarkers(object = seuset, only.pos = OnlyPos, min.pct = minpct, logfc.threshold = threshold, test.use = MarkerGeneMethod)
write.table(data.frame(gene = allmarkers$gene, allmarkers[,1:6]), file = file.path(outpath, paste0(prefix, "_GraphClust.AllMarkerGenes.txt")), sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(seuset, file = file.path(outpath, paste0(prefix, "_GraphClust.seuset.rds")))
cat(paste0("Graph-based clustering completed. Outputs saved to ", outpath, "\n"))

