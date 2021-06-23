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

# sampling 10,000 lines 
samp <- under15[,sample(dim(under15)[2], 10000)]
test <- pvclust(as.dist(1-cor(t(samp))),method.hclust="centroid",method.dist="cor")
