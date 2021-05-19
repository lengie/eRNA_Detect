## Written by Liana Engie, Dr. Scott E Fraser lab
## Last updated 5/11/2021

# Comparing nuclear noncoding RNA, bidirectional regions to ATAC-seq crossreferenced regions
# Uses sox10 nuclear RNA pulled from Trinh et al 2017 Biotagging paper: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89670
# Loads DeepTools computeMatrix tabular data and generates heat maps

# Input: computeMatrix from center of region, +/- 500bp on each side in 5bp bins. Each data set has a subset of bidirectional noncoding regions


library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(fastcluster)
options(scipen=999)

one <- fread("GalaxyComputeMatrix_sox10RepBidirCluster1Overlaps.tabular")

oneplus <- cbind(one[,1:200],one[,401:600])
colnames(oneplus) <- c(rep("sox10nuc1plus",200),rep("sox10nuc2plus",200))
oneminus <- cbind(one[,201:400],one[,601:800])
colnames(oneminus) <- c(rep("sox10nuc1minus",200),rep("sox10nuc2minus",200))

oneall <- one[,1:800]
colnames(oneall) <- c(rep("sox10nuc1plus",200),rep("sox10nuc2plus",200),rep("sox10nuc1minus",200),rep("sox10nuc2minus",200))

two <- fread("GalaxyComputeMatrix_sox10RepBidirCluster2Overlaps.tabular")
three <- fread("GalaxyComputeMatrix_sox10RepBidirCluster3Overlaps.tabular")
four <- fread("GalaxyComputeMatrix_sox10RepBidirCluster4Overlaps.tabular")
five <- fread("GalaxyComputeMatrix_sox10RepBidirCluster5Overlaps.tabular")
six <- fread("GalaxyComputeMatrix_sox10RepBidirCluster6Overlaps.tabular")
seven <- fread("GalaxyComputeMatrix_sox10RepBidirCluster7Overlaps.tabular")
eight <- fread("GalaxyComputeMatrix_sox10RepBidirCluster8Overlaps.tabular")
nine <- fread("GalaxyComputeMatrix_sox10RepBidirCluster9Overlaps.tabular")
ten <- fread("GalaxyComputeMatrix_sox10RepBidirCluster10Overlaps.tabular")
all <- fread("GalaxyComputeMatrix_sox10RepBidirAllClustersOverlaps.tabular")

GenerateHeatmap <- function(data,label1,label2,label3,label4){
    fread(data)
    span <- length(data)-1
    plustable <- cbind(data[,1:span/4],one[,span/2+1:3*span/4]
    oneplus <- oneplus %>% replace(is.na(.), -1)
    colnames(plustable) <- c(rep(label1,span/4),rep(label2,span/4))
    minustable <- cbind(data[,span/4+1:span/2],one[,3*span/4+1:span]
    oneminus <- oneminus %>% replace(is.na(.), -1)
    colnames(minustable) <- c(rep(label3,span/4),rep(label4,span/4))
    
    mypaletteP <- brewer.pal(9,"BuGn") 
    mypaletteM <- brewer.pal(9,"RdPu")
    morecolsP <- colorRampPalette(mypaletteP)
    morecolsM <- colorRampPalette(mypaletteM)
    
    x11() #code to open it full screen?
                        
    hr <- flashClust::hclust(as.dist(1-cor(t(oneall), method="pearson")), method="ward.D2")
    heatmap.2(as.matrix(oneall),Rowv=as.dendrogram(hr),col=rev(morecolsM(50)),main="Sox10 Bidir +/-strands row scaled Ward",scale="row",trace="none",Colv=NA)
    heatmap.2(as.matrix(oneall),Rowv=as.dendrogram(hr),col=rev(morecolsP(50)),main="Sox10 Bidir +/- strands row scaled Ward",scale="row",trace="none",Colv=NA)

    heatmap.2(as.matrix(oneplus),Rowv=as.dendrogram(hr),col=rev(morecolsM(50)),main="Sox10 Bidir + row scaled Ward",scale="row",trace="none",Colv=NA)
    heatmap.2(as.matrix(oneminus),Rowv=as.dendrogram(hr),col=rev(morecolsM(50)),main="Sox10 Bidir - row scaled Ward",scale="row",trace="none",Colv=NA)

                        
}
