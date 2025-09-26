
# Script: QuSAGE.R — QuSAGE analysis CLI for Seurat RDS and gene set files
suppressPackageStartupMessages({
  library(qusage)
  library(Seurat)
})

# Dependency check
required_pkgs <- c("qusage", "Seurat")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, FUN.VALUE = TRUE, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(paste0("Missing required packages: ", paste(missing_pkgs, collapse = ", "), 
              ". Please install them before running (e.g., install.packages())."))
}

# CLI: <input_rds> <geneset_files_comma_separated> <output_dir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript QuSAGE.R <input_rds> <geneset_files_comma_separated> <output_dir>\n")
  cat("Example: Rscript QuSAGE.R cells.rds H.all.v7.2.symbols.gmt,custom_sets.txt ./Result\n")
  quit(save = "no", status = 1)
}
file_rds <- args[1]
list_gmt <- strsplit(args[2], ",", fixed = TRUE)[[1]]
outdir <- args[3]
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Function: Get expression matrix and cluster list from rds file
# Input:    file_rds
# Output:   eset (expression matrix)
#           ident (list of clusters for each cell)
GetInfoFromRDS <- function (file_rds) {
  seuset <- readRDS(file_rds)
  date()
  # Extract expression matrix robustly across Seurat versions
  mat <- tryCatch(Seurat::GetAssayData(seuset, slot = "data"), error = function(e) NULL)
  if (is.null(mat)) mat <- tryCatch(seuset@data, error = function(e) NULL)
  if (is.null(mat)) stop("Could not extract expression matrix from Seurat object or RDS.")
  eset <<- as.matrix(mat)
  # Extract cluster identities robustly across Seurat versions
  idv <- tryCatch(as.character(Seurat::Idents(seuset)), error = function(e) NULL)
  if (is.null(idv)) idv <- tryCatch(as.character(seuset@ident), error = function(e) NULL)
  if (is.null(idv)) stop("Could not extract cluster identities from Seurat object.")
  ident <<- idv
  date()
  print("Loaded RDS successfully.")
}

# Function: Write gene set info into MSIG.genesets from gmt file or txt file with GeneSets info
# Input:    file_gmt
# Output:   MSIG.genesets
GetGeneSets <- function (file_gmt) {
  #print(paste('file_gmt:', file_gmt))
  if (grepl("txt$", file_gmt, perl = T, ignore.case = T)) {
    info <- read.table(file_gmt)
    colnames(info) <- c("type", "gene")
    MSIG.genesets <<- list()
    for (i in unique(info$type)) {
      set <- list(as.character(info[info$type==i,2]))
      MSIG.genesets <<- c(MSIG.genesets, set)
    }
    names(MSIG.genesets) <<- unique(info$type)
    #print(MSIG.genesets)
  }
  if ( grepl("gmt$", file_gmt, perl = T, ignore.case = T)) {
    MSIG.genesets <<- read.gmt(file_gmt)
    length(MSIG.genesets)
  }
  #MSIG.genesets
}


# Function: Create color schemes for plot
SetPlotArg <- function () {
  colors40 <<- c("darkblue", "orange", "green3", "red", "purple", "saddlebrown",
		"pink1", "darkgrey", "yellow", "deepskyblue", "darkgreen",
		"blue", "darkmagenta", "red4", "yellow4", "aquamarine2",
		"magenta1", "lemonchiffon1", "grey21", "tomato1", "hotpink4",
		"deeppink1", "deepskyblue4", "slategrey", "greenyellow",
		"mediumpurple", "darkorange1", "orchid3", "peachpuff1",
		"lightgreen", "cornflowerblue", "brown3", "tan", "darkorchid4",
		"hotpink", "blue3", "chartreuse4", "khaki1", "cyan", "salmon4")
  colors2  <<- c("black", rainbow(60))
  colors80 <<- c(colors40, rep("grey", 40))
  lwds50 <<- c(rep(5,49), 10)
}

# Instruction: This function is modified code based on plotDensityCurve from qusage library
plotDC <- function (QSarray, path.index = 1:numPathways(QSarray), zeroLine = TRUE,
		    addVIF = !is.null(QSarray$vif), col = NULL, plot = TRUE,
		    add = FALSE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
		    type = "l", ...) {
    if (is.character(path.index)) {
        path.index = match(path.index, names(QSarray$pathways))
    }
    if (all(is.na(QSarray$path.mean[path.index]))) {
        stop("no non-missing pathways in path.index")
    }
    scaleFactor = pdfScaleFactor(QSarray, addVIF = addVIF)
    if (is.null(xlim)) {
        xlim = range(calcBayesCI(QSarray, addVIF = addVIF)[, path.index], na.rm = T)
        xlim = xlim + (xlim - mean(xlim)) * 2
    }
    if (is.null(ylim)) {
        ylim = range(t(QSarray$path.PDF[, path.index])/scaleFactor[path.index], na.rm = T)
    }
    if (is.null(xlab))  xlab = "Gene Set Activity"
    if (is.null(ylab))  ylab = "Density"
    if (is.null(col))   col = par("col")
    if (!add & plot) {
        plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
    }
    if (length(col) != length(path.index)) {
        col = rep(col, length.out = length(path.index))
    }
    if (zeroLine & plot) abline(v = 0, lty = 2)
    retVal = list()
    for (i in 1:length(path.index)) {
        path = path.index[i]
        x = getXcoords(QSarray, path, addVIF = addVIF)
        y = QSarray$path.PDF[, path]/scaleFactor[path]
        if (plot) {
            lines(x, y, col = col[i], type = type, ...)
        }
        retVal[[i]] = data.frame(x, y)
    }
    names(retVal) = colnames(QSarray$path.PDF)[path.index]
    invisible(retVal)
}

# Instruction: This function is the original code from qusage library
pdfScaleFactor <- function(QSarray, addVIF=!is.null(QSarray$vif)){
  if(is.null(QSarray$vif) && addVIF){stop("vif is undefined for QSarray object. addVIF can not be set to true.")}
  sif = sapply(1:numPathways(QSarray),function(i){ifelse(addVIF,sqrt(QSarray$vif[i]),1)})
  sif[is.na(sif)] = 1
  pdfSum = colSums(QSarray$path.PDF)

  ##get the pdf range
  ranges = QSarray$ranges * 2

  ##the scale factor is essentially the distance between points in the x coordinates times the (current) sum of the pdf
  scaleFactor = (ranges*sif) / (QSarray$n.points-1) * pdfSum
  scaleFactor
}

# Function: Create QuSAGE results include DensityCurvesPlots, ConfidenceIntervalPlots and DataTable based on qs.results
# Input:    qs.results (output of qusage)
#           table_name (file name of DataTable)
#           plotDC_name(file name of DensityCurvesPlot)
#           plotCI_name(file name of ConfidenceIntervaPlots)
#           cluster    (cluster ID)
GetResultFiles <- function (
    qs.results,
    table_name,
    plotDC_name,
    plotCI_name,
    cluster
  ) {
  clusID <- paste("Cluster", cluster)
  p.vals <- pdf.pVal(qs.results)
  #class(p.vals)
  #print(p.vals)
  q.vals <- p.adjust(p.vals, method="fdr")
  #class(q.vals)
  topSets <- order(p.vals)
  #class(topSets)
  #print(topSets[1:25])

  #print(MSIG.genesets)
  table <- qsTable(qs.results, number = (length(MSIG.genesets)))
  #print(table[(1:5),(1:3)])
  write.table(table, file=table_name, sep="\t", row.names = F)
  table.order <- as.numeric(row.names(table))
  #print(table.order)


  ordered.legend = character()
  for (i in 1:length(table.order)) {
    ordered.legend[i] <- names(MSIG.genesets)[(table.order[i])]
    if (i > 5) next
    print(paste(table.order[i],names(MSIG.genesets)[table.order[i]],ordered.legend[i],sep=" "))
  }
  #print(names(MSIG.genesets))
  #print(ordered.legend)

  curve.num = length(table.order)
  top.num = curve.num

  png(file=plotDC_name, width = 1800, height = 1000)
  par(mgp = c(3,1,0), mar=c(6,6,3,2))
  plotDC(qs.results, path.index=table.order[top.num:1], lwd=3, col=colors80[top.num:1],
		    cex.axis=2, cex.xaxis=2, cex.lab=2, main=clusID, cex.main=2)
  plotDC(qs.results, path.index=table.order[1], lwd=9, col=colors80[1],
		    cex.axis=2, cex.xaxis=2, cex.lab=2,
		    add=T)
  legend("topleft", legend = ordered.legend[1:top.num], text.col = colors80, bty="n")#xpd=T)
  dev.off()

  pdf(file=sub("png", "pdf", plotDC_name), width = 27, height = 15)
  par(mgp = c(3,1,0), mar=c(6,6,3,3))
  plotDensityCurves(qs.results, path.index=table.order[top.num:1], lwd=5, col=colors80[top.num:1],
		    cex.axis=2, cex.xaxis=2, cex.lab=2,
		    xlab="Gene Set Activity", ylab="Density", main=clusID, cex.main=2)
  plotDensityCurves(qs.results, path.index=table.order[1],
		    lwd=10, col=colors80[1], cex.axis=2, cex.xaxis=2, cex.lab=2,
		    add=T)
  legend("topleft", legend = ordered.legend[1:top.num], text.col = colors80,
	 bty="o", bg="transparent") #xpd = T)
  dev.off()


  png(file=plotCI_name, width = 1800, height = 1000)
  par(mgp = c(6,1,0))
  plotCIs(qs.results, path.index=topSets[1:top.num], mar=c(30,15,3,1), lwd=3,
	  col=colors2[1:top.num], cex.axis=2, cex.xaxis=3, cex.lab=2,
	  ylab="Gene Set Activity", main=clusID, cex.main=2)
  dev.off()
}

# Function: Set output file names inclue DensityCurvesPlot, ConfidenceIntervalPlot and DataTable
# Input:    cluster (cluster name)
SetResultFileName <- function (cluster) {
  label_rds <<- "InputData"
  label_gmt <<- sub(".gmt$", "", file_gmt, perl = T, ignore.case = T)
  label_gmt <<- sub(".symbols$", "", label_gmt, perl = T, ignore.case = T)
  label_gmt <<- sub(".v6.2$", "", label_gmt, perl = T, ignore.case = T)
  label_gmt <<- sub("^.+/", "", label_gmt, perl = T, ignore.case = T)
  label <- paste(label_gmt, "Cluster", cluster, sep = "_")
  outwd <- file.path(outdir, "QuSAGE", label_gmt)
  if (!file.exists(outwd)) dir.create(outwd, recursive = TRUE)
  table_name  <<- file.path(outwd, paste0("table_",  label,  ".txt"))
  plotDC_name <<- file.path(outwd, paste0("plotDC_", label,  ".png"))
  plotCI_name <<- file.path(outwd, paste0("plotCI_", label,  ".png"))
}

main <- function () {
  print(date())
  date()
  # Load inputs via CLI
  GetInfoFromRDS(file_rds)
  SetPlotArg()
  for (f in list_gmt) {
    file_gmt <<- f
    print(file_gmt)
    GetGeneSets(file_gmt)
    for (i in unique(ident)) {
      print(paste("Cluster", i, date(), sep = " "))
      labels <- ident
      labels[labels!=i] = "B"
      labels[labels==i] = "C"
      contrast <- "C-B"
      qs.results <- qusage(eset, labels, contrast, MSIG.genesets)
      SetResultFileName(i)
      GetResultFiles(qs.results, table_name, plotDC_name, plotCI_name, i)
    }
    date()
  }
  cat(paste0("QuSAGE outputs saved under: ", file.path(outdir, "QuSAGE"), "\n"))
}

main()
