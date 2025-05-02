# Data import and preprocessing

library(Seurat)
library(dplyr)

# Parameters
cellranger_output_dir <- "path/to/cellranger/output"
sample_ids <- c("sample1", "sample2", "sample3")

# Create list to store Seurat objects
seurat_list <- list()

# Import 10x Genomics Cell Ranger output data
for (sample_id in sample_ids) {
  data_path <- file.path(cellranger_output_dir, sample_id, "outs", "filtered_feature_bc_matrix")
  counts <- Read10X(data.dir = data_path)
  
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  seurat_obj$sample <- sample_id
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_list[[sample_id]] <- seurat_obj
}

# Merge samples if multiple
if (length(seurat_list) > 1) {
  combined <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.ids = sample_ids,
    project = "pRCC_Project"
  )
} else {
  combined <- seurat_list[[1]]
}

# Save merged object
saveRDS(combined, file = "data/combined_raw.rds") 