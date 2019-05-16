### BidirNoncodRNA.r
###
### Purpose: Load RNA-seq data as BAM, remove coding regions via GTF annotations, merge strands into overlapped regions (get TPM of these regions?), explore these. 
###
###
### Written by Liana Engie
### Last updated: April 2019
###
### bidirncRNA(bamfile,gtffile)
### Input: string chromosome number, int input_start, int input_end, string strand (either "+" or "-")
### Output:

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments) 
library(ggplot2)
library(dplyr) 
library(data.table)
library(bedr)

bamfile <- "/auto/cmb-00/rr/engie/RNA/hets1.bam" 
gtffile <- "/auto/cmb-00/rr/engie/RNA/Danio_rerio.GRCz11.96.gtf" 
fragno <- 50.731011 

bamfile <- "/auto/cmb-00/rr/engie/RNA/hets2.bam" 

bidirncRNAwGTF{
    txdb <- makeTxDbFromGFF(gtffile,
                            format="gtf",
                            circ_seqs = character()
                           )
    seqlevelsStyle(txdb) <- "UCSC"
	coding <- cds(txdb)

	flag <- scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE)
    bamread <- readGAlignmentPairs(bamfile,
                                   param=ScanBamParam(flag=flag)
                                  )
    bamread <- granges(bamread)
	
	ncbam <- GenomicRanges::setdiff(bamread,coding,
                     ignore.strand=FALSE)
	#ignoring for now the * strand reads         
    #ncbam then creates a list of bamread that does not include any regions in coding
    
    ncbed <- data.frame(chr <- as.character(seqnames(ncbam)),
                        start <- as.integer(start(ncbam)-1),
                        end <- end(ncbam),
                        strand <- strand(ncbam)
                       )
    ncbed <- data.table(ncbed)
	colnames(ncbed) <- c("chr","start","end","strand")	#write.table(ncbed,file="noncodingHets1.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
	#merged <- bedr.merge.region(ncbed,verbose=FALSE)
	containschr <- grepl("chr",chr)
	justchr <- cbind(containschr,ncbed)
	justchr <- subset(justchr,containschr==TRUE)	
	remCol <- justchr[,containschr:=NULL]#write.table(remCol,file="noncodingHets1ONLYCHR.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
	#mergedchronly <- bedr.merge.region(remCol,verbose=FALSE)
	
    ## Finish merging in terminal
    
    outsidemerged <- fread("noncodingCHRONLY_mergedhets1.bed",data.table=TRUE,fill=TRUE)
	colnames(outsidemerged) <- c("chr","start","end")
    feat <- GRanges(seqnames=outsidemerged$chr,ranges=IRanges(start=outsidemerged$start,end=outsidemerged$end),ignore.strand=TRUE)
    bidir_read <- summarizeOverlaps(features=feat,reads=bamread,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE) 
    bidir_read_counts <- assay(bidir_read)
    outsidemerged <- cbind(outsidemerged,bidir_read_counts)
    outsidemerged <- mutate(outsidemerged, size = end-start)
    write.table(outsidemerged,file="ct711aHets1BidirRegions.csv",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

}

