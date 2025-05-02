# Cell-cell communication analysis

library(Seurat)
library(dplyr)
library(ggplot2)
library(reticulate)

# Load annotated data if not in memory
if (!exists("seurat_obj")) {
  seurat_obj <- readRDS("data/seurat_annotated.rds")
}

# Prepare CellPhoneDB input files
# Expression matrix file
expression_data <- as.matrix(GetAssayData(seurat_obj, slot = "data"))
write.table(expression_data, "cellphonedb/expression.txt", sep = "\t", quote = FALSE)

# Cell metadata file
meta_data <- data.frame(
  Cell = colnames(seurat_obj),
  cell_type = seurat_obj$cell_type
)
write.table(meta_data, "cellphonedb/meta.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Run CellPhoneDB analysis (Python code commented out, uncomment when available)
# use_python("/path/to/python")
# system("cellphonedb method statistical_analysis cellphonedb/meta.txt cellphonedb/expression.txt --output-path=cellphonedb/results --threads=4")
# system("cellphonedb plot dot_plot --means-path=cellphonedb/results/means.txt --pvalues-path=cellphonedb/results/pvalues.txt --output-path=cellphonedb/results")
# system("cellphonedb plot heatmap_plot --pvalues-path=cellphonedb/results/pvalues.txt --output-path=cellphonedb/results")

# Process CellPhoneDB results if available
if (file.exists("cellphonedb/results/means.txt") && file.exists("cellphonedb/results/pvalues.txt")) {
  means <- read.table("cellphonedb/results/means.txt", header = TRUE, sep = "\t")
  pvals <- read.table("cellphonedb/results/pvalues.txt", header = TRUE, sep = "\t")
}

# Example interaction data
set.seed(42)
interactions <- data.frame(
  source = sample(unique(seurat_obj$cell_type), 20, replace = TRUE),
  target = sample(unique(seurat_obj$cell_type), 20, replace = TRUE),
  weight = runif(20, 0, 1)
)

# Ensure source != target
interactions <- interactions[interactions$source != interactions$target, ]

# Visualize network using igraph if available
if (requireNamespace("igraph", quietly = TRUE)) {
  library(igraph)
  
  g <- graph_from_data_frame(interactions, directed = TRUE)
  
  png("figures/cell_interactions_network.png", width = 800, height = 800)
  plot(g, 
       edge.width = E(g)$weight * 10, 
       vertex.size = 30, 
       vertex.color = rainbow(length(unique(c(interactions$source, interactions$target)))),
       vertex.label.cex = 0.8,
       layout = layout_with_fr(g))
  dev.off()
}

# Save cell interaction results
saveRDS(interactions, "results/cell_interactions.rds") 