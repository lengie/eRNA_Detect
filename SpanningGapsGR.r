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

# overlap comparison
list50 <- findOverlaps(repl50,allclusters)
list100 <- findOverlaps(repl100,allclusters)
list500 <- findOverlaps(repl500,allclusters)

rm50NoNA <- rm50[queryHits(nearest50),]
rm100NoNA <- rm100[queryHits(nearest100),]
rm500NoNA <- rm500[queryHits(nearest500),]
replNoNA <- repl[queryHits(nearest),]

# make list of only indices that are close to each other
adj <- replNoNA[list,]
adj50 <- rm50NoNA[list50,]
adj100 <- rm100NoNA[list100,]
adj500 <- rm500NoNA[list500,]

# add flanks, but there will be some out of chrom bounds
adj$start <- adj$start - 150
adj$end <- adj$end + 150
adj50$start <- adj50$start -150
adj50$end <- adj50$end+150
adj100$start <- adj100$start -150
adj100$end <- adj100$end+150
adj500$start <- adj500$start -150
adj500$end <- adj500$end+150

# remove negative values
nonzero <- function(df){for(i in 1:nrow(df)){
    if(df$start[i]<1) df$start[i]=1}
return(df)}

# apply to all data frames
lapply(list(adj,adj50,adj100,adj500),function(x) x<-nonzero(x))

comb <- GenomicRanges::union(GRanges(adj),grepl)
comb50 <- GenomicRanges::union(GRanges(adj50),repl50)
comb100 <- GenomicRanges::union(GRanges(adj100),repl100)
comb500 <- GenomicRanges::union(GRanges(adj500),repl500)

# overlap comparison
overlap <- findOverlaps(comb,allclusters)
overlap50 <- findOverlaps(comb50,allclusters)
overlap100 <- findOverlaps(comb100,allclusters)
overlap500 <- findOverlaps(comb500,allclusters)

sumstat <- function(index){
    print(min(width(index)))
    print(median(width(index)))
    print(mean(width(index)))
    print(max(width(index)))
    print(length(index))}

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

grlist <- c(comb,comb50,comb100,comb500)
clist <- GRangesList(allclusters,cluster1,cluster2,cluster3,cluster4,cluster5,cluster6,cluster7,cluster8,cluster9,cluster10)
cnames <- c("allclusters","cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10")
for(j in 1:length(grlist)){
    for(i in 1:length(clist)){
        overlaps <- findOverlaps(grlist[j],clist[i],ignore.strand=TRUE)
        save <- grlist[j][queryHits(overlaps),]
        print(sumstat(save))
        file <- paste("sox10_Zv9FlankedWiderNoncodingReprodOverlaps150Flank",grlist[j],clist[i],".bed",sep="")
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
