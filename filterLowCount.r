### Read count, and read count by basepair, annd scatter plots
### 
### Last updated: October 2021

library(ggplot)
library(data.table)
library(GenomicRanges)



lowcountStacked <- function(bed,bam){
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
    justreadsplot <- data.frame(counts=c(genecounts[,1],genecounts[,2]),
                    strand=rep(c('+','-'),each=nrow(genecounts)))
    return(list(widthplot,justreadsplot))
}


lowcountStrand <- function(bed,bam){
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
    widthplot <- data.frame(rpb_p = genecounts[,1]/width(feat) , 
                            rpb_m = genecounts[,2]/width(feat))
    justreadsplot <- data.frame(counts_p=genecounts[,1],
                                counts_m=genecounts[,2])
    return(list(widthplot,justreadsplot))
}


png("sox10_868_bidirStrandCount.png",width=1200,height=1080)
    ggplot(toplot68[[1]], aes(x=V1, y=V2)) +
      geom_point(size=0.25) + xlab("plus") + ylab("minus") +
      ggtitle("sox10 ribo-seq rep2 bidir regions read counts by length") +
      scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10")
dev.off() 

png("biorepplot_exonsbidirCodingRem.png",width=1200,height=1080)
    ggplot(toplot, aes(x=nuc1, y=nuc2,color=strand)) +
      geom_point(size=0.25) + xlab("Nuc1") + ylab("Nuc2") +
      ggtitle("sox10 biological replicates, exons' and bidir regions' read counts by length") +
      scale_color_manual(name="Strand",
                         labels=c("Plus","Minus"),
                         values=c("tan3","cadetblue4"))
dev.off() 
 
png("biorepplot_exonsbidirCodingRem.png",width=1000,height=800)
    ggplot(toplot, aes(x=nuc1, y=nuc2,color=strand)) +
      geom_point(size=0.25) + xlab("Nuc1") + ylab("Nuc2") +
      ggtitle("sox10 biological replicates, exons' and bidir regions' read counts by length") +
      scale_color_manual(name="Strand",
                         labels=c("Plus","Minus"),
                         values=c("tan3","cadetblue4")) +
      scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10")
dev.off() 
