if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) 
    BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("enrichplot", quietly = TRUE))
    BiocManager::install("enrichplot")
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")

# Dependencies check (do not auto-install)
req_pkgs <- c("clusterProfiler","org.Hs.eg.db","dplyr")
missing <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop(paste("Missing packages:", paste(missing, collapse = ", "), "- please install them before running."))
}
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# CLI arguments: <gene_list.txt> [output_dir]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript ClusterProfile-GO.R <gene_list.txt> [output_dir]\n")
  cat("Example: Rscript ClusterProfile-GO.R genes.txt ./Result\n")
  quit(save = "no", status = 1)
}

gene_list_file <- args[1]
output_dir <- if (length(args) >= 2) args[2] else "./Result"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_file <- file.path(output_dir, "GO_enrichment.csv")


genes <- readLines(gene_list_file)
genes <- unique(genes)


gene_map <- bitr(genes, fromType = "SYMBOL",
                 toType = c("ENTREZID", "ENSEMBL"), 
                 OrgDb = org.Hs.eg.db)  


if (nrow(gene_map) == 0) {
  stop("No valid gene ID mapping found. Please check input gene symbols and the OrgDb.")
}

mapped_genes <- gene_map$ENTREZID

ego_BP <- enrichGO(gene          = mapped_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP", 
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_MF <- enrichGO(gene          = mapped_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_CC <- enrichGO(gene          = mapped_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

if (!is.null(ego_BP)) {
  ego_BP$ONTOLOGY <- "BP"
}
if (!is.null(ego_MF)) {
  ego_MF$ONTOLOGY <- "MF"
}
if (!is.null(ego_CC)) {
  ego_CC$ONTOLOGY <- "CC"
}

all_results <- rbind(as.data.frame(ego_BP), 
                    as.data.frame(ego_MF),
                    as.data.frame(ego_CC))

if (nrow(all_results) == 0) {
  all_results <- data.frame(
    ID = character(),
    Description = character(),
    ONTOLOGY = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    pAdjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer()
  )
}

write.csv(all_results, file = output_file, row.names = FALSE)
cat(paste0("GO enrichment results saved to: ", output_file, "\n"))