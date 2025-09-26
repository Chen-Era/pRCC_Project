# CLI: <MatrixDir> <sample.txt> <output_dir> [prefix=SampleName] [genecolumn=2] [cellMinGene=200] [cellMaxGene=10000] [geneExpMinCell=3] [MaxMTPercent=0.5] [regressOut=nCount_RNA,percent.mito]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript Seruat3Normalization.R <MatrixDir> <sample.txt> <output_dir> [prefix=SampleName] [genecolumn=2] [cellMinGene=200] [cellMaxGene=10000] [geneExpMinCell=3] [MaxMTPercent=0.5] [regressOut=nCount_RNA,percent.mito]\n")
  quit(save = "no", status = 1)
}

countsFile <- args[1]
samplelist <- args[2]
outpath <- args[3]
prefix <- if (length(args) >= 4) args[4] else "SampleName"
genecolumn <- if (length(args) >= 5) as.integer(args[5]) else 2
cellMinGene <- if (length(args) >= 6) as.integer(args[6]) else 200
cellMaxGene <- if (length(args) >= 7) as.integer(args[7]) else 10000
geneExpMinCell <- if (length(args) >= 8) as.integer(args[8]) else 3
MaxMTPercent <- if (length(args) >= 9) as.numeric(args[9]) else 0.5
regressOut <- if (length(args) >= 10) args[10] else "nCount_RNA,percent.mito"

NormType <- "scale"
ScaleModel <- "linear"
runFast <- TRUE
if (NormType=="scale"){
  pcadim = 1:10; computePCs = 50; nneighbors = 30; varGeneNum = 2000
} else if (NormType=="sctransform"){
  pcadim = 1:30; computePCs = 50; nneighbors = 30; varGeneNum = 3000
}

library(reshape2)
library(Seurat)
library(dplyr)
library(stringr)

print("Start")
if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)
setwd(outpath)

allcounts = Read10X(data.dir = countsFile, gene.column = genecolumn)

dir.create("1.Normalization", showWarnings = FALSE)
dir.create("2.PCAAnalysis", showWarnings = FALSE)

allcounts = allcounts[,colSums(allcounts)>0, drop = FALSE]

# Mito genes file expected at ./MT.txt (can be customized)
mito_file <- file.path(outpath, "MT.txt")
if (!file.exists(mito_file)) stop("Mito gene list file MT.txt not found in output directory.")
mito.genes = read.table(mito_file, sep = "\t", header = FALSE)[,1]
mito.genes = mito.genes[is.finite(match(mito.genes, rownames(allcounts)))]
percent.mito <- Matrix::colSums(allcounts[mito.genes, ]) / Matrix::colSums(allcounts)

seuset <- CreateSeuratObject(counts = allcounts, project = prefix, min.cells = geneExpMinCell, min.features = 1)

print("Remove mitochondrial genes")
seuset = subset(seuset, features = setdiff(row.names(seuset), mito.genes))
seuset <- AddMetaData(object = seuset, metadata = percent.mito, col.name = "percent.mito")

sampleInfo = read.table(samplelist, sep = "\t", header = FALSE, row.names = 1)
samples = sapply(strsplit(colnames(seuset), "-"), function(v) return(v[2]))
sampleInfo = sampleInfo[samples,]

names(sampleInfo) = colnames(seuset)
seuset <- AddMetaData(object = seuset, metadata = sampleInfo, col.name = "orig.ident")
write.table(data.frame(CellName = rownames(seuset@meta.data), seuset@meta.data), file = "1.Normalization/CellInfosRaw.txt", sep = "\t", quote = FALSE, row.names = FALSE)

print(paste0("nGene filter range:", cellMinGene, "-", cellMaxGene))
seuset <- subset(x = seuset, subset = nFeature_RNA > cellMinGene & nFeature_RNA < cellMaxGene)
seuset <- subset(x = seuset, subset = percent.mito < MaxMTPercent)

seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = varGeneNum)

ccRDS = file.path(outpath, "CCGenes.rds")
if (!file.exists(ccRDS)) stop("Cell cycle gene file CCGenes.rds not found in output directory.")
ccgene = readRDS(ccRDS)
seuset = CellCycleScoring(seuset, s.features = ccgene$s.genes, g2m.features = ccgene$g2m.genes, set.ident = FALSE)
seuset$CC.Difference <- seuset$S.Score - seuset$G2M.Score

if(NormType=="scale"){
  top10 <- head(VariableFeatures(seuset), 10)
  write.table(data.frame(Gene = rownames(seuset@hvg.info), seuset@hvg.info), file = "1.Normalization/hvgInfo.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(VariableFeatures(seuset), file = "1.Normalization/varGene.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  scalegene = VariableFeatures(seuset)
  print("Use variable genes to scale data")
  seuset <- ScaleData(object = seuset, features = scalegene, vars.to.regress = regressOut, model.use = ScaleModel)
}
write.table(data.frame(CellName = rownames(seuset@meta.data), seuset@meta.data), file = "1.Normalization/CellInfosFilter.txt", sep = "\t", row.names = FALSE, quote = FALSE)

seuset <- RunPCA(object = seuset, features = VariableFeatures(seuset), npcs = computePCs, do.print = TRUE, ndims.print = 1:computePCs, nfeatures.print = 10, seed.use = 2021)

pcgene <- c()
for(i in 1:computePCs){
    topgenes = TopFeatures(object = seuset[["pca"]], dim = i, balanced = TRUE, nfeatures = 60)
    topgenes = data.frame(PC = paste0("PC", i), topgenes$positive, topgenes$negative)
    pcgene = rbind(pcgene, topgenes)
}
write.table(pcgene, file = "2.PCAAnalysis/pcGene.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.frame(CellName = rownames(seuset@reductions$pca@cell.embeddings), seuset@reductions$pca@cell.embeddings), file = "2.PCAAnalysis/pca.txt", sep = "\t", quote = FALSE, row.names = FALSE)
print("Run UMAP/TSNE")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors, reduction.name = "umap3d", n.components = 3, reduction.key = "umap3d_")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors)
seuset <- RunTSNE(object = seuset, dims = pcadim, reduction.name = "tsne3d", dim.embed = 3, reduction.key = "tSNE3d_", check_duplicates = FALSE)
seuset <- RunTSNE(object = seuset, dims = pcadim, check_duplicates = FALSE)

saveRDS(seuset, file = file.path(outpath, paste0(prefix, ".seuset.rds")))
cat(paste0("Normalization and PCA completed. Outputs saved to ", outpath, "\n"))
