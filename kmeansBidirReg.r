### kmeansBidirReg.r
###
### Purpose: Take a bed file of regions of interest, as well as the createMatrix tabular data of transcription across bins spanning the ROIs, 
###          and output k-means clustering of the transcriptional profile 
###
###
### Written by Liana Engie
### Last updated: June 2021
###
### 
### Input:
### (Current) Output: 

library(GenomicFeatures)
library(GenomicAlignments) 
library(ggplot2)
library(dplyr) 
library(data.table)
library(pvclust)
library(parallel)
library(dendextend)
library(stringr)
options(scipen=999)

##regions under 1.5kb
# createMatrix data of transcription in bins spanning length of genomic location
under <- fread("Galaxy381-[100Rem500bpNearestSpanned_under1500bp_1.4kb_100bp__Heatmap_values].tabular",header=FALSE)
under15 <- t(under) %>% replace(is.na(.),0)

#1400bp centered on middle of region
binsun <- c(paste(seq(-7,-1,by=1),"00bp_plus",sep=""),
          paste(seq(1,7,by=1),"00bp_plus",sep=""),
          paste(seq(-7,-1,by=1),"00bp_minus",sep=""),
          paste(seq(1,7,by=1),"00bp_minus",sep="")) 
row.names(under15) <- binsun

# labeling by number to track
colnames(under15) <- paste("BidirReg",1:dim(under15)[2],sep="")
                      
# bed file of actual ROI genomic locations, minus one that Galaxy/deepTools had issue with
onefive <- fread("sox10_Zv9FlankedWiderNoncodingReprod100Rem500bpNearestSpanned<1.5kb.bed")
colnames(onefive) <- c("chr","start","end","ID","score","strand")
onefiveless <- onefive[-1356503,]
underGR <- GRanges(onefiveless)

# cluster list from sox10 paper
sox10_1 <- fread('ATACCrossRef/Sox10nuclear_cluster1.counts.txt')
sox10_2 <- fread('ATACCrossRef/Sox10nuclear_cluster2.counts.txt')
sox10_3 <- fread('ATACCrossRef/Sox10nuclear_cluster3.counts.txt')
sox10_4 <- fread('ATACCrossRef/Sox10nuclear_cluster4.counts.txt')
sox10_5 <- fread('ATACCrossRef/Sox10nuclear_cluster5.counts.txt')
sox10_6 <- fread('ATACCrossRef/Sox10nuclear_cluster6.counts.txt')
sox10_7 <- fread('ATACCrossRef/Sox10nuclear_cluster7.counts.txt')
sox10_8 <- fread('ATACCrossRef/Sox10nuclear_cluster8.counts.txt')
sox10_9 <- fread('ATACCrossRef/Sox10nuclear_cluster9.counts.txt')
sox10_10 <- fread('ATACCrossRef/Sox10nuclear_cluster10.counts.txt')
cluster1 <- GRanges(sox10_1)
cluster2 <- GRanges(sox10_2)
cluster3 <- GRanges(sox10_3)
cluster4 <- GRanges(sox10_4)
cluster5 <- GRanges(sox10_5)
cluster6 <- GRanges(sox10_6)
cluster7 <- GRanges(sox10_7)
cluster8 <- GRanges(sox10_8)
cluster9 <- GRanges(sox10_9)
cluster10 <- GRanges(sox10_10)

clist <- GRangesList(cluster10,cluster9,cluster8,cluster7,cluster6,cluster5,cluster4,cluster3,cluster2,cluster1)
cnames <- c("cluster10","cluster9","cluster8","cluster7","cluster6","cluster5","cluster4","cluster3","cluster2","cluster1") 
for(i in 1:10){
    # overlap between the bed file and each cluster
    overlaps <- findOverlaps(underGR,clist[[i]],ignore.strand=TRUE)
    if(length(overlaps)>0){
        # change the column name of the samples that overlap, so they're not 'bidirregion#' but "clusteri_#' and we can see how cluster overlaps are distributed in the new k-means 
        colnames(under15)[queryHits(overlaps)] <- paste(cnames[i],"_",1:length(overlaps),sep="")
        
        #save the list of regions that overlap with each cluster, as one regino may overlap with several clusters
        x <- paste(cnames[i],"list",sep="")
        assign(x,underGR[queryHits(overlaps),])
    }
}

# remove all columns with no reads at all
noZ <- under15[, colSums(under15 != 0) > 0]

# set colors for leaves of dendrogram
clusterno <- colnames(noZ)
clusterno <- replace(clusterno,grepl("Bidir",clusterno),11)
clusterno <- gsub('cluster([0-9]+)_.*','\\1',clusterno)
clustercol <- stringr::str_replace_all(clusterno, 
                                        setNames(c("gray","purple","black","gold","orange","red","magenta","green","green4","cyan","cornflowerblue"),
                                                 as.character(c(11,10,1:9))))

run.pvclust.iterations <- function(matrix,colorvector,iterations){
    for(i in 1:iterations){
        seed <- runif(1, min = 100, max = 999)
        set.seed(seed)
        print(seed) #maybe make a not-verbose option

        ind <- sample(dim(noZ)[2], 1000)
        samp <- noZ[,ind]
        colors_to_use <- clustercol[ind]
        clust <- pvclust(samp,method.hclust="centroid",method.dist="cor")
        plotd <- as.dendrogram(test)
        colors_to_use <- colors_to_use[order.dendrogram(plotd)]
        dendnew <- assign_values_to_leaves_edgePar(dend=plotd, value = colors_to_use, edgePar = "col") %>% set("labels_cex", 0)  
    
        pngfile <- paste("pvclust_ReprodUnd1500bpRandSamp",seed,".png",sep="")
            plot(dendnew,main="pvclust with leaves clustered by color") 
        dev.off()
    }
}


run.kmeans.iterations <- function(matrix,iterations){
    for(i in 1:iterations){
        seed <- runif(1, min = 100, max = 999)
        set.seed(seed)
        print(seed) 
        ind <- sample(dim(matrix)[2], 1000)
        samp <- noZ[,ind]
        clustercat <- clustercol[,ind]

        kClust <- kmeans(t(samp), centers=9, nstart = 100, iter.max = 10)
        pngfile <- paste("kmeans_ReprodUnd1500bpIter10N10k_",i,".png",sep="")
        png(pngfile,)
            fviz_cluster(kClusters, data=samp, geom="point", repel=TRUE, shape=clustercat) + ggtitle("k = 2") + scale_shape_manual('sox10 cluster designation', values=c(22,23,24))
        dev.off()
    }
}
