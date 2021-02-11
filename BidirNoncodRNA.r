### BidirNoncodRNA.r
###
### Purpose: Load RNA-seq data as BAM, remove coding regions via GTF annotations, merge strands into overlapped regions (get TPM of these regions?), explore these. 
###
###
### Written by Liana Engie
### Last updated: February 2021
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

bam1file <- "ct711a_150804_hets_nuc1PrimaryReads.bam"
bamfile2 <- "/auto/cmb-00/rr/engie/RNA/hets2.bam" 
gtffile <- "/auto/cmb-00/rr/engie/RNA/Danio_rerio.GRCz11.96.gtf" 

bidirncRNAwGTF{
	txdb <- makeTxDbFromGFF(gtffile,
                               format="gtf",
                               circ_seqs = character()
                               )
	seqlevelsStyle(txdb) <- "UCSC"
	coding <- cds(txdb)

	flag <- scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE)
    	bamread <- readGAlignmentPairs(bamfile, param=ScanBamParam(flag=flag))
     	gbam <- GRanges(bamread)
	
   	overlaps <- findOverlaps(gbam1,coding,ignore.strand=FALSE) #could put the GRanges() in here to save storage space
	hits1 <- gbam1[-queryHits(overlaps),]
	underten1 <- hits1[width(hits1)<10000]
	ten1 <- hits1[width(hits1)>10000]
	write.table(ten1,file="/panfs/qcb-panasas/engie/GRCz11Star/ct711a_Hets1PrimaryNoncodingOver10k.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    	
	#needs to be a 6 column bed
	underten1df <- data.frame(chr = as.character(seqnames(underten1)),
                    		  start = start(underten1)-1,
                    	 	  end = end(underten1),
				  name = "n/a",
				  score = 0,
				  strand = strand(underten1)
				  )
	write.table(underten1df,file="/panfs/qcb-panasas/engie/GRCz11Star/ct711a_150804_hets_nuc1PrimaryNoncodingUnder10k.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

## Sort and merge outside of R
bedtools sort -i ct711a_150804_hets_nuc1PrimaryNoncodingUnder10k.bed > ct711a_150804_hets_nuc1PrimaryNoncodingUnder10kSorted.bed
bedtools merge -i ct711a_150804_hets_nuc1PrimaryNoncodingUnder10kSorted.bed -d 0 -s > ct711a_150804_hets_nuc1PrimaryNoncodingUnder10kMerged.bed
	
## Back in R
	outsidemerged1 <- fread("GRCz11Star/ct711a_150804_hets_nuc1PrimaryNoncodingUnder10kMerged.bed",data.table=TRUE,fill=TRUE)
	plus <- dplyr::filter(outsidemerged1,V4=="+")
	minus <- dplyr::filter(outsidemerged1,V4=="-")

	plus6 <- data.frame(plus$V1,plus$V2,plus$V3,ID="n/a",score=0,strand=plus$V4)
	minus6 <- data.frame(minus$V1,minus$V2,minus$V3,ID="n/a",score=0,strand=minus$V4)
	write.table(plus6,file="GRCz11Star/ct711a_150804_hets_nuc1PrimaryUnder10kMergedPlus.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
	write.table(minus6,file="GRCz11Star/ct711a_150804_hets_nuc1PrimaryUnder10kMergedMinus.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

## Intersect outside of R
intersectBed -sorted -S -a ct711a_150804_hets_nuc1PrimaryUnder10kMergedMinus.bed -b ct711a_150804_hets_nuc1PrimaryUnder10kMergedPlus.bed > ct711a_150804_hets_nuc1PrimaryUnder10kMergedIntersected.bed

 ## Back in R   
	intersected1 <- fread("GRCz11Star/ct711a_150804_hets_nuc1PrimaryUnder10kMergedIntersected.bed",data.table=TRUE,fill=TRUE)
	colnames(intersected1) <- c("chr","start","end","ID","score","strand")
	
	feat <- GRanges(seqnames=intersected1$chr,ranges=IRanges(start=intersected1$start,end=intersected1$end))  
	bidir_read1 <- summarizeOverlaps(features=feat,reads=bam1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
	bidir_read_counts1 <- assay(bidir_read1)
	intersected1 <- cbind(intersected1,bidir_read_counts1)
	intersected1 <- dplyr::mutate(intersected1, size = end-start)
	write.table(intersected1,file="ct711aHets1PrimaryReadsBidirRegions.csv",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t") 
}

#outside of R
bedtools sort -i ct711a_150804_hets_nuc1PrimaryUnpairedNoncoding.bed > ct711a_150804_hets_nuc1PrimaryUnpairedNoncodingSorted.bed 
bedtools merge -i ct711a_150804_hets_nuc1PrimaryUnpairedNoncodingSorted.bed -d 0 -s > ct711a_150804_hets_nuc1PrimaryUnpairedNoncodingMerged.bed

#saving after non-coding for graphical check
chr1 <- data.frame(chr=seqnames(gbam),
                    start <- as.integer(start(gbam)-1),
                    end <- end(gbam),
                    strand <- strand(gbam)
                    )
chr1$chr <- paste("chr",chr1$chr,sep="")

#saving after merging for graphical check
merged1 <- dplyr::filter(outsidemerged1,V1==1)
merged1$V1 <- paste("chr",merged1$V1,sep="")
merged1 <- data.frame(merged1$V1,merged1$V2,merged1$V3,ID="n/a",score=0,strand=merged1$V4)
write.table(merged1,file="GRCz11Star/ct711a_150804_hets_nuc1PrimaryUnder10kMergedChr1.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

#saving after intersecting for graphical check
one1 <- dplyr::filter(intersected1,chr==1)
one1$chr <- paste("chr",one1$chr,sep="")
write.table(one1,file="ct711a_150804_hets_nuc1PrimaryUnder10kMergedIntersectedChr1.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

intersectedchr1 <- dplyr::filter(intersected1,chr==1)
intersectedchr1 <- data.frame(chr=paste("chr",intersectedchr1$chr,sep=""),
				    start=intersectedchr1$start,
				    end=intersectedchr1$end,
				    score=intersectedchr1$reads,
				    ID=intersectedchr1$size,
				    strand=intersectedchr1$strand)
write.table(intersectedchr1,file="GRCz11Star/ct711aHets1PrimaryUnder10kBidirRegionsChr1.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

