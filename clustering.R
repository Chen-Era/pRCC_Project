# Cell clustering and marker gene identification

library(Seurat)
library(dplyr)
library(ggplot2)

# Load processed data if not in memory
if (!exists("seurat_obj")) {
  seurat_obj <- readRDS("data/seurat_processed.rds")
}

# Build nearest neighbor graph
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)

# Cluster cells using Louvain algorithm at different resolutions
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
}

# Select appropriate resolution
Idents(seurat_obj) <- "RNA_snn_res.0.6"

# Visualize clustering results
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("UMAP clustering results")
ggsave("figures/umap_clusters.png", p1, width = 10, height = 8)

p2 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE) + ggtitle("t-SNE clustering results")
ggsave("figures/tsne_clusters.png", p2, width = 10, height = 8)

# Identify marker genes for each cluster
all_markers <- FindAllMarkers(seurat_obj, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = 0.25, 
                             test.use = "wilcox")

# Save all marker genes
write.csv(all_markers, "results/all_cluster_markers.csv", row.names = FALSE)

# Select top 10 marker genes per cluster
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Generate heatmap of top marker genes
p3 <- DoHeatmap(seurat_obj, features = unique(top10_markers$gene)) + NoLegend()
ggsave("figures/top10_markers_heatmap.png", p3, width = 12, height = 10)

# Generate violin plots for selected marker genes
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

p4 <- VlnPlot(seurat_obj, features = unique(head(top5_markers$gene, 9)), ncol = 3)
ggsave("figures/markers_violin_plots.png", p4, width = 15, height = 12)

# Generate feature plots for selected marker genes
p5 <- FeaturePlot(seurat_obj, features = unique(head(top5_markers$gene, 9)), ncol = 3)
ggsave("figures/markers_feature_plots.png", p5, width = 15, height = 12)

# Create cell type annotation mapping
cell_type_annotations <- data.frame(
  cluster = 0:10,
  cell_type = c("T cells", "B cells", "Monocytes", "NK cells", "Macrophages", 
                "Dendritic cells", "Epithelial cells", "Endothelial cells", "Fibroblasts", 
                "Erythrocytes", "Stromal cells")
)

# Add annotations to Seurat object
seurat_obj$cell_type <- plyr::mapvalues(
  x = Idents(seurat_obj),
  from = as.character(cell_type_annotations$cluster),
  to = cell_type_annotations$cell_type
)

# Save annotation file
write.csv(cell_type_annotations, "results/cell_type_annotations.csv", row.names = FALSE)

# Plot with cell type annotations
Idents(seurat_obj) <- "cell_type"
p6 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("UMAP cell type annotations")
ggsave("figures/umap_cell_types.png", p6, width = 10, height = 8)

# Save annotated Seurat object
saveRDS(seurat_obj, file = "data/seurat_annotated.rds") 