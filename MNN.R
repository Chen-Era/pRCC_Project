# Dependencies check
req_pkgs <- c("Seurat","reshape2","dplyr","ggthemes","stringr","SingleCellExperiment","scater","scran","batchelor")
missing <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop(paste("Missing packages:", paste(missing, collapse = ", "), "- please install them before running."))
}
library(Seurat)
library(reshape2)
library(dplyr)
library(ggthemes)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(scran)
library(batchelor)

# CLI: <input_rds> <output_dir> [minCell=10] [k=5] [resolution=0.8]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript MNN.R <input_rds> <output_dir> [minCell=10] [k=5] [resolution=0.8]\n")
  quit(save = "no", status = 1)
}

inrds <- args[1]
output_dir <- args[2]
minCell <- if (length(args) >= 3) as.integer(args[3]) else 10
kValue <- if (length(args) >= 4) as.integer(args[4]) else 5
resolution <- if (length(args) >= 5) as.numeric(args[5]) else 0.8
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

set.seed(2023)

pcadim = 1:10
computePCs = 50
nneighbors = 30
varGeneNum = 2000

resolution = 0.8

seuset= readRDS(inrds)
regressvar = seuset@commands$ScaleData.RNA$vars.to.regress

seuset.list <- SplitObject(seuset, split.by = "orig.ident")

cellNum = lapply(seuset.list, ncol)
print(cellNum)
print(lapply(seuset.list,nrow))
seuset.list = seuset.list[cellNum>minCell]

if (length(x = seuset.list) < 2) {
    stop("RDS file must contain multiple Seurat objects for MNN", call. = FALSE)
}
seuset = merge(
    x = seuset.list[[1]],
    y = seuset.list[2:length(x = seuset.list)]
)
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = varGeneNum)
scalegene = VariableFeatures(seuset)
seuset <- ScaleData(object = seuset, features = scalegene, vars.to.regress = regressvar)

objects.sce <- lapply(
    X = seuset.list,
    FUN = function(x, f) { 
      return(as.SingleCellExperiment(x = subset(x = x, features = VariableFeatures(seuset))))
    },
    f = VariableFeatures(seuset)
)

print("RunMNN")
out <- do.call(fastMNN, c(objects.sce, list(k = kValue, d = computePCs)));

outcorrected = reducedDim(out)
rownames(outcorrected)=colnames(seuset)
colnames(outcorrected)=paste0("MNN_",1:ncol(outcorrected))
featureloading = as.matrix(rowData(out))
colnames(featureloading)=paste0("MNN_",1:ncol(featureloading))
seuset[["mnn"]] <- CreateDimReducObject(
    embeddings = outcorrected,
    loadings = featureloading,
    assay = DefaultAssay(object = seuset),
    key = "mnn_"
)

print("RunUMAP")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors, reduction.name = "umap3d", n.components = 3, reduction.key = "umap3d_",reduction="mnn")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors,reduction="mnn")

print("RunTSNE")
seuset <- RunTSNE(object = seuset,dims = pcadim,reduction.name = "tsne3d",dim.embed=3, reduction.key = "tSNE3d_",check_duplicates = FALSE,reduction="mnn")
seuset <- RunTSNE(object = seuset,dims = pcadim,check_duplicates = FALSE,reduction="mnn")

seuset <- FindNeighbors(seuset, reduction = "mnn", dims = pcadim)
seuset <- FindClusters(seuset, resolution = resolution)

# Outputs
write.table(data.frame(Cell = colnames(seuset), Cluster = Idents(seuset)),
            file = file.path(output_dir, "MNN_Clusters.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(seuset, file = file.path(output_dir, "seuset.mnn.rds"))
cat(paste0("MNN integration completed. Outputs saved to ", output_dir, "\n"))

