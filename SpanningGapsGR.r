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

library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments) 
library(ggplot2)
library(dplyr) 
library(data.table)
options(scipen=999)

## Useful functions
sumstat <- function(index){
    print(min(width(index)))
    print(median(width(index)))
    print(mean(width(index)))
    print(max(width(index)))
    print(length(index))}

# function to remove negative values
nonzeroGR <- function(gr){
    neg <- which(start(gr)<1)
    start(gr)[neg] <- 1
return(gr)}


#make sure no coding regions were brought back in
remcoding <- function(gr){
    hits <- GenomicRanges::setdiff(gr,coding,ignore.strand=TRUE)
    hits <- hits[!width(hits)==1]
    return(hits)
}

#make sure that adding flanking regions doesn't go beyond chromosome length
chrLimitCheckNoIntGR <- function(gr,limit){
    keep <- GRanges()
    for(i in 1:nrow(limit)){
        tmp <- gr[which(seqnames(gr)==as.character(limit$chrom[i])),]
        over <- which(end(tmp)>limit$size[i]))
        end(tmp)[over] <- limit$size[i]
        keep <- c(keep,tmp)
    }
    return(keep)
}
limit <- fread("danRer7.chrom.sizes")
colnames(limit) <- c("chrom","size")

##
# loading reproducible bidirectional regions and cluster list
repl <- fread("sox10_Zv9FlankedWiderNoncodingReprodOnly.bed")
colnames(repl) <- c("chr","start","end")
repl <- mutate(repl,size=end-start)
allten <- fread("ATACCrossRef/sox10_BiotaggingAllClusters.bed")
colnames(allten) <- c("chr","start","end","size","strand","counts","size")
allclusters <- GRanges(allten)


codes <- fread("Zv9ExonsUTRs500bpFlanking.bed")
colnames(codes) <- c("chr","start","end","size","strand","ID","score") 
coding <- GRanges(codes)

# this code adds a flanking buffer region to each read and combines the regions that now overlap. so all regions that do not merge with another now have a buffer
flankSpanGap <- function(repl,threshold,gap){
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

#this code does the same, but does not use the internet and starts with the GRanges object
flankSpanGapGR <- function(repl,gap){ 
    nearest <- distanceToNearest(repl)
    list <- which(mcols(nearest)$distance<gap)
    #remove NA
    rmNoNA <- repl[queryHits(nearest),]
    
    # add flanks to those close to each other
    adj <- rmNoNA[list,]
    start(adj) <- start(adj) - gap  #could have these be separate numbers if folks want to have sep parameters
    end(adj) <- end(adj) + gap
    #make sure we don't have negative indices
    adj <- nonzeroGR(adj)
    adj <- chrLimitCheckNoIntGR(adj,limit)
    
    #combine flanked regions into original and make sure there's no coding regions added back in
    comb <- GenomicRanges::union(adj,grepl)
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


spanNearestGap <- function(repl,threshold,gap,regionLimit){
    # remove reads under threshold size
    rm <- dplyr::filter(repl,size>threshold)
    
    # which neighbor is nearest, and how far?
    replgr <- GRanges(rm)
    ndist <- GenomicRanges::distanceToNearest(replgr)
    adj <- GRanges(rm[queryHits(ndist),]) #should have no NAs
    nearest <- GenomicRanges::nearest(adj) #should have no gaps from NAs  
    
    for(i in 1:length(nearest)){
        if(mcols(adj)[i,]<regionLimit && mcols(adj)[nearest[i],]<regionLimit && nearest[i]<gap){
            #preceding is nearest
            if(i>nearest[i]){ 
                start(adj)[i] = start(adj)[nearest(i)]
                }
            elseif(i<nearest[i]) #following is the nearest
                end(adj)[i] = end(adj)[nearest(i)] #will capture if several are close to each other in a row
        }
    }
    
    #make sure we don't have negative indices
    adj <- nonzero(adj)
    adj <- chrLimitCheck(adj,"danRer7")
    
    #combine flanked regions into original and make sure there's no coding regions added back in
    comb <- GenomicRanges::union(GRanges(adj),replgr)
    comb <- remcoding(comb)
    return(comb)
}
