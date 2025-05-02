# Data normalization and dimensionality reduction

library(Seurat)
library(dplyr)
library(ggplot2)

# Load filtered data if not in memory
if (!exists("combined_filtered")) {
  combined_filtered <- readRDS("data/combined_filtered.rds")
}

# Normalize data
seurat_obj <- NormalizeData(combined_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Plot variable features
p1 <- VariableFeaturePlot(seurat_obj)
top10 <- head(VariableFeatures(seurat_obj), 10)
p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
ggsave("figures/variable_features.png", p2, width = 10, height = 8)

# Scale data
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes, vars.to.regress = c("nCount_RNA", "percent.mt"))

# Run PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Visualize PCA results
p3 <- DimPlot(seurat_obj, reduction = "pca")
ggsave("figures/pca_plot.png", p3, width = 8, height = 6)

# Generate PCA heatmap
p4 <- DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
ggsave("figures/pca_heatmap.png", p4, width = 12, height = 10)

# Determine PC number
p5 <- ElbowPlot(seurat_obj, ndims = 40)
ggsave("figures/pca_elbow_plot.png", p5, width = 8, height = 6)

# Use first 10 PCs for subsequent analysis
n_pcs <- 10

# Run UMAP and t-SNE
seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:n_pcs)

# Visualize UMAP results
p6 <- DimPlot(seurat_obj, reduction = "umap")
ggsave("figures/umap_plot.png", p6, width = 8, height = 6)

# Visualize t-SNE results
p7 <- DimPlot(seurat_obj, reduction = "tsne")
ggsave("figures/tsne_plot.png", p7, width = 8, height = 6)

# Save processed object
saveRDS(seurat_obj, file = "data/seurat_processed.rds") 