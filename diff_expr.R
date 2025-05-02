# Differential gene expression analysis

library(Seurat)
library(dplyr)
library(ggplot2)

# Load annotated data if not in memory
if (!exists("seurat_obj")) {
  seurat_obj <- readRDS("data/seurat_annotated.rds")
}

# Set identity based on cell type
Idents(seurat_obj) <- "cell_type"

# Define cell type pairs for comparison
cell_type_comparisons <- list(
  c("T cells", "B cells"),
  c("Macrophages", "Monocytes"),
  c("NK cells", "T cells")
)

# Perform differential expression analysis for each pair
for (comp in cell_type_comparisons) {
  # Check if both cell types exist
  if (all(comp %in% unique(seurat_obj$cell_type))) {
    # Perform differential expression analysis
    markers <- FindMarkers(
      seurat_obj,
      ident.1 = comp[1],
      ident.2 = comp[2],
      min.pct = 0.1,
      logfc.threshold = 0.25,
      test.use = "wilcox"
    )
    
    # Add gene name column
    markers$gene <- rownames(markers)
    
    # Save results
    write.csv(markers, 
             file = sprintf("results/diff_expr_%s_vs_%s.csv", 
                           gsub(" ", "_", comp[1]), 
                           gsub(" ", "_", comp[2])),
             row.names = FALSE)
    
    # Create volcano plot
    markers$diff_expressed <- "Not Significant"
    markers$diff_expressed[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.25] <- "Upregulated"
    markers$diff_expressed[markers$p_val_adj < 0.05 & markers$avg_log2FC < -0.25] <- "Downregulated"
    
    # Add labels for extreme genes
    markers$label <- NA
    top_genes <- markers %>%
      filter(p_val_adj < 0.05) %>%
      group_by(diff_expressed) %>%
      top_n(10, wt = abs(avg_log2FC)) %>%
      pull(gene)
    
    markers$label[markers$gene %in% top_genes] <- markers$gene[markers$gene %in% top_genes]
    
    # Volcano plot
    p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diff_expressed)) +
      geom_point(alpha = 0.7) +
      geom_text_repel(aes(label = label), na.rm = TRUE, size = 3, max.overlaps = 20) +
      scale_color_manual(values = c("Not Significant" = "grey", 
                                    "Upregulated" = "red", 
                                    "Downregulated" = "blue")) +
      theme_minimal() +
      labs(title = sprintf("Differential Expression: %s vs %s", comp[1], comp[2]),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value",
           color = "Expression Change")
    
    ggsave(sprintf("figures/volcano_%s_vs_%s.png", 
                  gsub(" ", "_", comp[1]), 
                  gsub(" ", "_", comp[2])),
           p, width = 10, height = 8)
    
    # Heatmap of top differentially expressed genes
    top_markers <- markers %>%
      filter(p_val_adj < 0.05) %>%
      top_n(20, wt = abs(avg_log2FC))
    
    if(nrow(top_markers) > 0) {
      cells_1 <- WhichCells(seurat_obj, idents = comp[1])
      cells_2 <- WhichCells(seurat_obj, idents = comp[2])
      cells <- c(cells_1, cells_2)
      
      p2 <- DoHeatmap(seurat_obj, 
                     features = top_markers$gene, 
                     cells = cells, 
                     label = TRUE) +
        ggtitle(sprintf("Heatmap: %s vs %s", comp[1], comp[2]))
      
      ggsave(sprintf("figures/heatmap_%s_vs_%s.png", 
                    gsub(" ", "_", comp[1]), 
                    gsub(" ", "_", comp[2])),
             p2, width = 12, height = 10)
    }
  }
}

# Compare specific cell type vs all others
target_cell_type <- "T cells"
if (target_cell_type %in% unique(seurat_obj$cell_type)) {
  markers_vs_all <- FindMarkers(
    seurat_obj,
    ident.1 = target_cell_type,
    ident.2 = NULL,  # Compare with all other cell types
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  
  markers_vs_all$gene <- rownames(markers_vs_all)
  
  # Save results
  write.csv(markers_vs_all, 
           file = sprintf("results/diff_expr_%s_vs_others.csv", 
                         gsub(" ", "_", target_cell_type)),
           row.names = FALSE)
  
  # Create dot plot
  top_markers_vs_all <- markers_vs_all %>%
    filter(p_val_adj < 0.05) %>%
    top_n(20, wt = avg_log2FC)
  
  if(nrow(top_markers_vs_all) > 0) {
    p3 <- DotPlot(seurat_obj, features = top_markers_vs_all$gene) + 
      coord_flip() +
      ggtitle(sprintf("%s Specific Genes", target_cell_type))
    
    ggsave(sprintf("figures/dotplot_%s_vs_others.png", 
                  gsub(" ", "_", target_cell_type)),
           p3, width = 12, height = 10)
  }
}

# Save summary of differential expression analysis
summary <- data.frame(
  comparison = paste0(unlist(lapply(cell_type_comparisons, function(x) paste0(x[1], " vs ", x[2])))),
  significant_genes = NA
)

for (i in 1:nrow(summary)) {
  comp <- strsplit(as.character(summary$comparison[i]), " vs ")[[1]]
  file_path <- sprintf("results/diff_expr_%s_vs_%s.csv", 
                      gsub(" ", "_", comp[1]), 
                      gsub(" ", "_", comp[2]))
  
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    summary$significant_genes[i] <- sum(data$p_val_adj < 0.05)
  }
}

write.csv(summary, "results/diff_expr_summary.csv", row.names = FALSE) 