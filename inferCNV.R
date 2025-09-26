library(ggplot2)

# CLI: <observations.txt> <references.txt> <group.txt> <output_dir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript inferCNV.R <observations.txt> <references.txt> <group.txt> <output_dir>\n")
  quit(save = "no", status = 1)
}

observationsFile <- args[1]
referencesFile <- args[2]
groupFile <- args[3]
outdir <- args[4]
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

##### use CNV data to predict cancer cells and classify clusters ----->>>>>
cnv_result= read.csv( file = observationsFile, stringsAsFactors=FALSE ,sep=" ",check.name=F)
# this "cnv_result_2" is the CNV data including the cells chosen as reference (non-malignant cells)
cnv_result_2=read.csv(file=referencesFile,stringsAsFactors=FALSE ,sep=" ",check.name=F)

# the cnv result need to minus 1, considering it treats 1 as no normalization
cnv_result=cnv_result-1
cnv_result_2=cnv_result_2-1
cnv_result_combine=cbind(cnv_result,cnv_result_2)

cnv_meanSquare=apply(cnv_result_combine^2, 2, mean)

cancer_average=apply(cnv_result,1,mean)
cnv_corr=apply(cnv_result_combine,2,function(x) cor(x, cancer_average, method = "pearson"))
cnv_corr[is.na(cnv_corr)]=0

table(names(cnv_corr)==names(cnv_meanSquare))
# the name is not in the same order
# change the order of "cnv_corr" as that of names of "cnv_meanSquare"
cnv_corr<-cnv_corr[order(match(names(cnv_corr),names(cnv_meanSquare)))]
mergeOne=cbind(cnv_corr,cnv_meanSquare)
write.table(mergeOne, file.path(outdir, "CNV_cells.txt"), quote=F, sep="\t")
pdf(file.path(outdir, "correlation_against_meanSquare.pdf"), width=10, height=10, paper='special')
smoothScatter(cnv_corr, cnv_meanSquare, xlab="correlation", ylab="meanSquare", nrpoints=1000, cex=0.5)
abline(h=0.05,  lty=2,col="dodgerblue",bty="n")
abline(v=0.5,  lty=2, col="darkgreen",bty="n")
legend("left", legend=c("SNV signal 0.05","correlation with tumor cells 0.5"), col=c( "dodgerblue", "darkgreen"), lty=2, cex=1.5, bty="n")
dev.off()

#cell summary
group = read.delim(groupFile, stringsAsFactors=F,header=F)
celllist=group[,1]
cnv_meanSquare = cnv_meanSquare[order(match(names(cnv_meanSquare),celllist))]
cnv_result_combine = cnv_result_combine[order(match(colnames(cnv_result_combine),celllist))]

meandata = aggregate(cnv_meanSquare, list(group[,2]), FUN=mean)
cluster_mean = aggregate(t(as.matrix(cnv_result_combine)), list(group[,2]), FUN=mean)

cluster_corr = apply(cluster_mean[,-1],1,function(x) cor(x, cancer_average, method = "pearson"))
cnv_cluster_result=cbind(meandata,cluster_corr)
write.table(cnv_cluster_result, file.path(outdir, "CNV_clusters.txt"), quote=F, sep="\t")
gg=ggplot(cnv_cluster_result,aes(x=cnv_cluster_result[,2],y=cnv_cluster_result[,3]))+geom_point(shape=21,size=3)+geom_text(aes(y = cnv_cluster_result[,3] + .006, label = cnv_cluster_result[,1]))+geom_hline(yintercept = 0.5)+geom_vline(xintercept = 0.05)
ggplot2::ggsave(file.path(outdir, "CNVPlot.png"), gg, width=10, height=10)
cat(paste0("inferCNV summary saved to ", outdir, "\n"))
