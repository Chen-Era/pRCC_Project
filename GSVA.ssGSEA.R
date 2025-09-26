# Dependencies
req_pkgs <- c("GSVA","limma","stringr")
missing <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop(paste("Missing packages:", paste(missing, collapse = ", "), "- please install them before running."))
}
library(GSVA)
library(limma)
library(stringr)

# CLI: <gene_set.gmt|txt> <expression_matrix.txt> [output_dir]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript GSVA.ssGSEA.R <gene_set.gmt|txt> <expression_matrix.txt> [output_dir]\n")
  cat("Example: Rscript GSVA.ssGSEA.R gene.gmt expData.txt ./Result\n")
  quit(save = "no", status = 1)
}

gmtFile <- args[1]
inputFile <- args[2]
output_dir <- if (length(args) >= 3) args[3] else "./Result"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_file <- file.path(output_dir, "ssgsea.txt")

gmtread = function (file)
{
     if (!grepl("\\.gmt$", file)[1]) {
        data = read.table(file,header = T,sep="\t")
	types = as.character(unique(data[,1]))
	geneSetDB = list()
	for(i in 1:length(types)){
		type = types[i]
		geneset = as.character(unique(data[data[,1]==type,2]))
		geneSetDB[[i]]=geneset
	}
	names(geneSetDB) = types 
     }else{
     	geneSetDB = readLines(file)
     	geneSetDB = strsplit(geneSetDB, "\t")
     	names(geneSetDB) = sapply(geneSetDB, "[", 1)
    	geneSetDB = lapply(geneSetDB, "[", -1:-2)
     	geneSetDB = lapply(geneSetDB, function(x) {
        	x[which(x != "")]
     	})
     }
     return(geneSetDB)
}

geneset = gmtread(gmtFile)

data <- read.table(inputFile, sep = "\t", header = TRUE, check.names = FALSE)

gsva_matrix <- gsva(data, geneset, method = "ssgsea", ssgsea.norm = FALSE)
gsva_matrix = t(gsva_matrix)

write.table(gsva_matrix, file = output_file, sep = "\t", quote = FALSE)
cat(paste0("ssGSEA results saved to: ", output_file, "\n"))

