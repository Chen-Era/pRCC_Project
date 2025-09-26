# KEGG pathway enrichment script

# Dependencies check (do not auto-install)
req_pkgs <- c("clusterProfiler","org.Hs.eg.db","dplyr")
missing <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop(paste("Missing packages:", paste(missing, collapse = ", "), "- please install them before running."))
}
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# CLI arguments: <gene_list.txt> [output_dir] [organism]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript ClusterProfile-Pathway.R <gene_list.txt> [output_dir] [organism:hsa]\n")
  cat("Example: Rscript ClusterProfile-Pathway.R genes.txt ./Result hsa\n")
  quit(save = "no", status = 1)
}

gene_list_file <- args[1]
output_dir <- if (length(args) >= 2) args[2] else "./Result"
organism <- if (length(args) >= 3) args[3] else "hsa"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_file <- file.path(output_dir, "KEGG_enrichment.csv")

# 读取基因列表
genes <- readLines(gene_list_file)
genes <- unique(genes)  # 确保基因ID唯一

# 基因ID映射为ENTREZID
gene_map <- bitr(genes, fromType = "SYMBOL",  # 根据输入基因ID类型修改fromType
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)  # 根据物种修改OrgDb

# 检查是否有映射的基因
if (nrow(gene_map) == 0) {
  stop("No valid gene ID mapping found. Please check input gene symbols and the OrgDb.")
}

mapped_genes <- gene_map$ENTREZID

kk <- enrichKEGG(gene = mapped_genes,
                 organism = organism,
                 keyType = 'kegg',
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

if (!is.null(kk)) {
  kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
} else {
  kk <- data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    pAdjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer()
  )
}

write.csv(as.data.frame(kk), file = output_file, row.names = FALSE)
cat(paste0("KEGG enrichment results saved to: ", output_file, "\n"))