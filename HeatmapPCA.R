## Written by Liana Engie, Dr. Scott E Fraser lab
## Last updated 5/11/2021

# Comparing nuclear noncoding RNA, bidirectional regions to ATAC-seq crossreferenced regions through heatmaps of both strands
# Uses sox10 nuclear RNA pulled from Trinh et al 2017 Biotagging paper: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89670
# Loads DeepTools computeMatrix tabular data and generates heat maps

# Input: computeMatrix from center of region, +/- 500bp on each side in 5bp bins. Each data set has a subset of bidirectional noncoding regions

library(data.table)
nofile <- fread("GalaxyComputeMatrix_sox10Reprod_Bidir_50Rm_500Flank_NoClusters_3kb_50bp.tabular")
c2file <- fread("GalaxyComputeMatrix_sox10Reprod_Bidir_50Rm_500Flank_Cluster2_3kb_50bp.tabular")

bins <- rep(c(paste(seq(-30,-1,by=1),"plus",sep=""),
            paste(seq(1,30,by=1),"plus",sep=""),
            paste(seq(-30,-1,by=1),"minus",sep=""),
            paste(seq(1,30,by=1),"minus",sep="")),
            2)
            
matrixtest <- rbind(c2file[,1:240],nofile[,1:240],use.names=FALSE)
col <- c(paste("cluster2",1:nrow(c2file),sep="_"),paste("nocluster",1:nrow(nofile),sep="_"))

row.names(matrixtest) <- col
colnames(matrixtest) <- bins

mds <- cmdscale(dist(matrixtest), k=3, eig=TRUE)
mds$eig
# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)
# plot the PCA
png(file="PCA_PropExplainedVariance.png")
barplot(eig_pc,
     las=1,
     xlab="Dimensions", 
     ylab="Proportion of explained variance (%)", y.axis=NULL,
     col="darkgrey")
dev.off()

mds <- cmdscale(dist(countmatrix2)) 
Warning message:
In dist(countmatrix2) : NAs introduced by coercion

png(file="PCA_Dim1vsDim2.png")
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
dev.off()
