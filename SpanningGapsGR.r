### SpanningGapsGR -- spanGaps()
###
### Purpose: Load RNA-seq data as BAM, remove coding regions via GTF annotations, merge strands into union-overlapped regions, explore these. 
###
###
### Written by Liana Engie
### Last updated: June 2021
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

#make sure no coding regions were brought back in
remcoding <- function(gr){
    hits <- GenomicRanges::setdiff(gr,coding,ignore.strand=TRUE)
    hits <- hits[!width(hits)==1]
    return(hits)
}

#make sure that adding flanking regions doesn't go beyond chromosome length
#(note: need to make this faster than nested for loops)
chrLimitCheck <- function(df,gen){
    temp <- c()
    limit <- getChromInfoFromUCSC(gen)
    for(i in 1:nrow(df)){
        for(j in 1:nrow(limit)){
            if(df$chr[i]==limit$chrom[j] & df$end[i]>limit$size[j]){
                df$end[i]=limit$size[j]  
               if(df$start[i]>limit$size[j]) temp=c(temp,i)
            }
        if(is.null(temp)==FALSE) df <- df[-temp,]
        }
    }
    return(df)
} 

#having some null issues on the above
chrLimitCheck <- function(df,gen){
    temp <- c()
    limit <- getChromInfoFromUCSC(gen)
    for(i in 1:nrow(df)){
        for(j in 1:nrow(limit)){
            if(df$chr[i]==limit$chrom[j] & df$end[i]>limit$size[j]){
                df$end[i]=limit$size[j]  
            }
        }
    }
    return(df)
} 


##
# loading reproducible bidirectional regions and cluster list
repl <- fread("sox10_Zv9FlankedWiderNoncodingReprodOnly.bed")
colnames(repl) <- c("chr","start","end")
repl <- mutate(repl,size=end-start)
allten <- fread("ATACCrossRef/sox10_BiotaggingAllClusters.bed")
colnames(allten) <- c("chr","start","end","size","strand","counts","size")
allclusters <- GRanges(allten)

spanGap <- function(repl,threshold,gap){
    # remove reads under threshold size
    rm <- dplyr::filter(repl,size>threshold)
    
    # distance to neighbors
    replgr <- GRanges(rm)
    nearest <- distanceToNearest(replgr)
    list <- which(mcols(nearest)$distance<gap)
    #remove NA
    rmNoNA <- rm[queryHits(nearest),]
    
    # add flanks to those close to each other
    adj <- rmNoNA[list,]
    adj$start <- adj$start - gap  #could have these be separate numbers if folks want to have sep parameters
    adj$end <- adj$end + gap
    #make sure we don't have negative indices
    adj <- nonzero(adj)
    adj <- chrLimitCheck(adj,"danRer7")
    
    #combine flanked regions into original and make sure there's no coding regions added back in
    comb <- GenomicRanges::union(GRanges(adj),grepl)
    comb <- remcoding(comb)
    return(comb)
}


## should make separate function for cluster overlap comparison       
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
