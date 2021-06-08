### BidirNoncodRNA.r
###
### Purpose: Load RNA-seq data as BAM, remove coding regions via GTF annotations, merge strands into union-overlapped regions, explore these. 
###
###
### Written by Liana Engie
### Last updated: May 2021
###
### bidirncRNA(bamfile,gtffile)
### Input: string chromosome number, int input_start, int input_end, string strand (either "+" or "-")
### (Current) Output: bed6 file containing non-coding regions where RNA is read from both strands, consistently between two biological replicates

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments) 
library(ggplot2)
library(dplyr) 
library(data.table)
library(bedr)
options(scipen=999)

library(dplyr) 
library(GenomicFeatures)
library(GenomicRanges)
library(data.table)
library(ggplot2)
options(scipen=999)

## Useful functions
sumstat <- function(index){
    print(min(width(index)))
    print(median(width(index)))
    print(mean(width(index)))
    print(max(width(index)))
    print(length(index))}

# function to remove negative values
nonzero <- function(df){for(i in 1:nrow(df)){
    if(df$start[i]<1) df$start[i]=1}
return(df)}

##
# loading reproducible bidirectional regions and cluster list
repl <- fread("sox10_Zv9FlankedWiderNoncodingReprodOnly.bed")
colnames(repl) <- c("chr","start","end")
repl <- mutate(repl,size=end-start)
allten <- fread("ATACCrossRef/sox10_BiotaggingAllClusters.bed")
colnames(allten) <- c("chr","start","end","size","strand","counts","size")
allclusters <- GRanges(allten)

# removing regions under a certain length
rm50 <- dplyr::filter(repl,size>50)
rm100 <- dplyr::filter(repl,size>100)
rm500 <- dplyr::filter(repl,size>500)

repl50 <- GRanges(rm50)
repl100 <- GRanges(rm100)
repl500 <- GRanges(rm500)

grepl <- GRanges(repl)
# generate distance to nearest neighbor within reproducible bidirectional regions, thresholded
nearest <- distanceToNearest(grepl)
nearest50 <- distanceToNearest(repl50)
nearest100 <- distanceToNearest(repl100)
nearest500 <- distanceToNearest(repl500)

# finding out which have neighbors within 130bp
list <- which(mcols(nearest)$distance<130)
list50 <- which(mcols(nearest50)$distance<130)
list100 <- which(mcols(nearest100)$distance<130)
list500 <- which(mcols(nearest500)$distance<130)

# remove NA from nearest neighbor calculations
rm50NoNA <- rm50[queryHits(nearest50),]
rm100NoNA <- rm100[queryHits(nearest100),]
rm500NoNA <- rm500[queryHits(nearest500),]
replNoNA <- repl[queryHits(nearest),]

# make list of only indices that are close to each other
adj <- replNoNA[list,]
adj50 <- rm50NoNA[list50,]
adj100 <- rm100NoNA[list100,]
adj500 <- rm500NoNA[list500,]

# add flanks to those particular regions, but there will be some out of chrom bounds
adj$start <- adj$start - 150
adj$end <- adj$end + 150
adj50$start <- adj50$start -150
adj50$end <- adj50$end+150
adj100$start <- adj100$start -150
adj100$end <- adj100$end+150
adj500$start <- adj500$start -150
adj500$end <- adj500$end+150


# apply to all data frames
lapply(list(adj,adj50,adj100,adj500),function(x) x<-nonzero(x))

# combine together the old bidirectional regions and the new flanked regions
comb <- GenomicRanges::union(GRanges(adj),grepl)
comb50 <- GenomicRanges::union(GRanges(adj50),repl50)
comb100 <- GenomicRanges::union(GRanges(adj100),repl100)
comb500 <- GenomicRanges::union(GRanges(adj500),repl500)

#make sure no coding regions were brought back in
remcoding <- function(gr){
    hits <- GenomicRanges::setdiff(gr,coding,ignore.strand=TRUE)
    hits <- hits[!width(hits)==1]
    return(hits)
}
comb <- remcoding(comb)
comb50 <- remcoding(comb50)
comb100 <- remcoding(comb100)
comb500 <- remcoding(comb500)
       
# overlap comparison with all clusters and get unique "ground truth" regions that are covered
overlap <- findOverlaps(comb,allclusters)
overlap50 <- findOverlaps(comb50,allclusters)
overlap100 <- findOverlaps(comb100,allclusters)
overlap500 <- findOverlaps(comb500,allclusters)

unique <- unique(subjectHits(overlap))
length(unique)
length(unique)/length(allclusters)
#save indices 
write.table(allclusters[unique,],"sox10_Zv9FlankingWiderNoncodingReprodOverlaps150FlankCapturedClusters.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(allclusters[-unique,],"sox10_Zv9FlankingWiderNoncodingReprodOverlaps150FlankUncapturedClusters.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)

unique <- unique(subjectHits(overlap50))
length(unique)
length(unique)/length(allclusters)
write.table(allclusters[unique,],"sox10_Zv9FlankingWiderNoncodingReprodOverlaps50Rem150FlankCapturedClusters.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(allclusters[-unique,],"sox10_Zv9FlankingWiderNoncodingReprodOverlaps50Rem150FlankUncapturedClusters.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)

unique <- unique(subjectHits(overlap100))
length(unique)
length(unique)/length(allclusters)
write.table(allclusters[unique,],"sox10_Zv9FlankingWiderNoncodingReprodOverlaps100Rem150FlankCapturedClusters.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(allclusters[-unique,],"sox10_Zv9FlankingWiderNoncodingReprodOverlaps100Rem150FlankUncapturedClusters.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)

unique <- unique(subjectHits(overlap500))
length(unique)
length(unique)/length(allclusters)
write.table(allclusters[unique,],"sox10_Zv9FlankingWiderNoncodingReprodOverlaps500Rem150FlankCapturedClusters.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(allclusters[-unique,],"sox10_Zv9FlankingWiderNoncodingReprodOverlaps500Rem150FlankUncapturedClusters.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)

grlist <- GRangesList(comb,comb50,comb100,comb500)
int <- c("","50bpRem","100Rem","500Rem")
clist <- GRangesList(allclusters,cluster1,cluster2,cluster3,cluster4,cluster5,cluster6,cluster7,cluster8,cluster9,cluster10)
cnames <- c("allclusters","cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10") 
for(j in 1:length(grlist)){
    for(i in 1:length(clist)){
        overlaps <- findOverlaps(grlist[j],clist[i],ignore.strand=TRUE)
        save <- grlist[j][queryHits(overlaps),]
        print(sumstat(save))
        file <- paste("sox10_Zv9FlankedWiderNoncodingReprodOverlaps150Flank",int[j],cnames[i],".bed",sep="")
        write.table(save,file,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
    }
 }

# finding out which have neighbors within 500bp
list <- which(mcols(nearest)$distance<500)
list50 <- which(mcols(nearest50)$distance<500)
list100 <- which(mcols(nearest100)$distance<500)
list500 <- which(mcols(nearest500)$distance<500)

# make list of only indices that are within 500bp of each other
adj <- replNoNA[list,]
adj50 <- rm50NoNA[list50,]
adj100 <- rm100NoNA[list100,]
adj500 <- rm500NoNA[list500,]

# add flanks
adj$start <- adj$start - 500
adj$end <- adj$end + 500
adj50$start <- adj50$start - 500
adj50$end <- adj50$end + 500
adj100$start <- adj100$start - 500
adj100$end <- adj100$end + 500
adj500$start <- adj500$start - 500
adj500$end <- adj500$end + 500
adj <- nonzero(adj)
adj50 <- nonzero(adj50)
adj500<- nonzero(adj500)
adj100<- nonzero(adj100)

comb <- GenomicRanges::union(GRanges(adj),grepl)
comb50 <- GenomicRanges::union(GRanges(adj50),repl50)
comb100 <- GenomicRanges::union(GRanges(adj100),repl100)
comb500 <- GenomicRanges::union(GRanges(adj500),repl500)

# overlap comparison
overlap <- findOverlaps(comb,allclusters)
overlap50 <- findOverlaps(comb50,allclusters)
overlap100 <- findOverlaps(comb100,allclusters)
overlap500 <- findOverlaps(comb500,allclusters)

unique <- unique(subjectHits(overlap))
length(unique)
length(unique)/length(allclusters)

unique <- unique(subjectHits(overlap50))
length(unique)
length(unique)/length(allclusters)

unique <- unique(subjectHits(overlap100))
length(unique)
length(unique)/length(allclusters)

unique <- unique(subjectHits(overlap500))
length(unique)
length(unique)/length(allclusters)
