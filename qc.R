# Quality control

library(Seurat)
library(dplyr)
library(ggplot2)

# Load data if not in memory
if (!exists("combined")) {
  combined <- readRDS("data/combined_raw.rds")
}

# Plot QC metrics distribution
p1 <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/qc_violin_plots.png", p1, width = 12, height = 6)

# Plot feature relationships
p2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p3 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
p4 <- FeatureScatter(combined, feature1 = "nFeature_RNA", feature2 = "percent.mt")

CombinePlots(plots = list(p2, p3, p4), ncol = 3)
ggsave("figures/qc_feature_relationships.png", width = 15, height = 5)

# Filter low-quality cells
# 保留表达200-5000个基因，UMI计数<30000，线粒体基因比例<20%的细胞
combined_filtered <- subset(combined, 
                           subset = nFeature_RNA > 200 & 
                                    nFeature_RNA < 5000 & 
                                    nCount_RNA < 30000 & 
                                    percent.mt < 20)

# Plot filtered QC metrics
p5 <- VlnPlot(combined_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/qc_filtered_violin_plots.png", p5, width = 12, height = 6)

# Save filtered data
saveRDS(combined_filtered, file = "data/combined_filtered.rds") 