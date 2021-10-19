### Read count, and read count by basepair, annd scatter plots
### 
### Last updated: October 2021

library(ggplot)
library(data.table)
library(GenomicRanges)



lowcount <- function(bed,bam){
    feat <- GRanges(bed)

    # the bam files are in Ensembl style not UCSC
    seqlevelsStyle(feat) <- "ensembl"

    strand(feat) <- '+'
    counts <- summarizeOverlaps(features=feat,reads=bam,ignore.strand=FALSE) 
    genecounts <- assay(counts)
    strand(feat) <- '-'
    counts <- summarizeOverlaps(features=feat,reads=bam,ignore.strand=FALSE)      
    genecounts <- data.frame(plus=genecounts,minus=assay(counts))

    print(paste("Number of regions: ",dim(genecounts),sep="")) #this is just the number of bidirectional regions
    print(paste("Regions with with 0 reads on plus strand: ",length(which(genecounts[,1]==0)),sep=""))
    print(paste("Regions with with 0 reads on minus strand: ",length(which(genecounts[,2]==0)),sep=""))
    widthplot <- data.frame(rpb=c(genecounts[,1]/width(feat) , genecounts[,2]/width(feat)),
                    strand=rep(c('+','-'),each=nrow(genecounts)))
    justreadsplot <- data.frame(rpb=c(genecounts[,1],genecounts[,2]),
                    strand=rep(c('+','-'),each=nrow(genecounts)))
    return(list(widthplot,justreadsplot))
}
