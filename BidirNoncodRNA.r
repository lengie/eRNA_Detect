### BidirNoncodRNA.r
###
### Purpose: Load RNA-seq data as BAM, remove coding regions via GTF annotations, merge strands into union-overlapped regions, explore these. 
###
###
### Written by Liana Engie
### Last updated: August 2021
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
options(scipen=999)

# load bam files to process
bam1file <- "ct711a_150804_hets_nuc1PrimaryReads.bam"
bamfile2 <- "/auto/cmb-00/rr/engie/RNA/hets2.bam" 
gtffile <- "/auto/cmb-00/rr/engie/RNA/Danio_rerio.GRCz11.96.gtf" 


loadbg <- function(bg,strand){
    bg <- mutate(bg,strand)
    colnames(bg) <- c("chr","start","end","score","strand")
    return(bg)
}

codingGRange <- function(bamfile,gtffile,flank=500){
	#load coding regions to remove (coding exons and UTRs)
	txdb <- makeTxDbFromGFF(gtffile,
                               format="gtf",
                               circ_seqs = character()
                               )
	seqlevelsStyle(txdb) <- "UCSC"
	exons <- exons(txdb)
	utr5 <- fiveUTRsByTranscript(txdb)
	utr3 <- threeUTRsByTranscript(txdb)
 
	#adding a flank (default 500bp) to the UTRs
	utr5df <- as.data.frame(utr5)
	utr5minus <- dplyr::filter(utr5df,strand=="-")
	utr5plus <- dplyr::filter(utr5df,strand=="+")
	utr5minus <- data.frame(seqnames=utr5minus$seqnames,
				start=utr5minus$start,
				end=utr5minus$end+flank,
				strand=utr5minus$strand,
				exon_id=utr5minus$exon_id,
				exon_name=utr5minus$exon_name)
	utr5plus <- data.frame(seqnames=utr5plus$seqnames,
			       start=utr5plus$start-flank,
			       end=utr5plus$end,
			       strand=utr5plus$strand,
			       exon_id=utr5plus$exon_id,
			       exon_name=utr5plus$exon_name)

	utr3df <- as.data.frame(utr3)
	utr3minus <- dplyr::filter(utr3df,strand=="-")
	utr3plus <- dplyr::filter(utr53f,strand=="+")
	utr3minus <- data.frame(seqnames=utr3minus$seqnames,
				start=utr3minus$start-500,
				end=utr3minus$end,
				strand=utr3minus$strand,
				exon_id=utr3minus$exon_id,
				exon_name=utr3minus$exon_name)
	utr3plus <- data.frame(seqnames=utr5minus$seqnames,
			       start=utr5minus$start,
			       end=utr5minus$end+500,
			       strand=utr5minus$strand,
			       exon_id=utr5minus$exon_id,
			       exon_name=utr5minus$exon_name)
	
	#combing together all regions you want to remove from the bamfile
	coding <- c(exons,GRanges(utr5plus),GRanges(utr5minus),GRanges(utr3plus),GRanges(utr3minus))
	return(coding)
}

#read in bam file
flag <- scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE)
bamread <- readGAlignmentPairs(bamfile, param=ScanBamParam(flag=flag))
gbam <- GRanges(bamread)

ncOverlaps <- function(gbam,coding,filename){
	# find overlaps, remove them
   	overlaps <- findOverlaps(gbam,coding,ignore.strand=FALSE)
	hits <- gbam[-queryHits(overlaps),]
	# save the non-coding regions
	underten <- hits[width(hits)<10000]
	ten <- hits[width(hits)>10000]
	write.table(ten,file=paste(filename,"Over10k.bed",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    	
	#needs to be a 6 column bed
	undertendf <- data.frame(chr = as.character(seqnames(underten)),
                    		  start = start(underten)-1,
                    	 	  end = end(underten),
				  ID = 1:length(underten),
				  score = 0,
				  strand = strand(underten)
				  )
	#save the file
	write.table(undertendf,file=paste(filename,"Under10k.bed",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
	return(undertendf)
	}

## Sort and merge the individual strands outside of R, for each replicate
sortMerge <- function(filename,dist=0){
    system(sprintf("bedtools sort -i %s.bed > %sSorted.bed",filename,filename))
    system(paste("bedtools merge -i ",filename,"Sorted.bed -d ",dist," -S - > ",filename,"MergedMinus.bed",sep=""))
    system(paste("bedtools merge -i ",filename,"Sorted.bed -d ",dist," -S + > ",filename,"MergedPlus.bed",sep=""))
}

#merge the intersections for the wider
read_format <- function(file){
	strand <- fread(plusfile,data.table=TRUE,fill=TRUE)
	colnames(strand) <- c("seqnames","start","end","strand")
	return(GRanges(strand))
}
plus_file <- "SET_THESE_MergedPlus.bed"
minus_file <- "YOUR_DATASET_MergedMinus.bed"
plus <- read_format(plus_file)
minus <- read_format(minus_file)

overlap_format <- function(plus_strand,minus_strand,file){
	merge <- mergeByOverlaps(minus_strand,plus_strand,ignore.strand=TRUE)
	merged <- data.frame(chr=c(seqnames(merge$plus_strand),seqnames(merge$minus_strand)),
                        start=c(start(merge$plus_strand),start(merge$minus_strand)),
                        end=c(end(merge$plus_strand),end(merge$minus_strand)),
                        ID=1:nrow(merge),
                        score=0,
                        strand=c(strand(merge$plus_strand),strand(merge$minus_strand)))
        filename <- paste(file,"OverlapsToMerge.bed",sep="")
	write.table(merged,file=filename,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

    system(sprintf("bedtools sort -i %s > %sOverlapsMergedSorted.bed",filename,file))
    system(paste("bedtools merge -i ",file,"OverlapsMergedSorted.bed -d ",dist," > ",file,"OverlapsMerged.bed",sep=""))
}

replicateReprod <- function(){
	nuc1bidir <- fread("sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsMerged.bed")
	nuc2bidir <- fread("sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsMerged.bed")
	colnames(nuc1bidir) <- c("chr","start","end")
	colnames(nuc2bidir) <- c("chr","start","end")

	nuc1 <- fread("GSM2386488_sox10_nuc_combined.bed")
	nuc2 <- fread("GSM2386489_sox10_nuc_combined.bed")
	colnames(nuc1) <- c("chr","start","end","ID","score","strand")
	colnames(nuc2) <- c("chr","start","end","ID","score","strand")
	sox10_nuc1 <- GRanges(nuc1)
	sox10_nuc2 <- GRanges(nuc2)
	
	# make sure merge did not add an coding regions
	overlaps <- findOverlaps(X,coding,ignore.strand=FALSE)
	hits <- X[-queryHits(overlaps),]
	#REPEAT 

	feat <- GRanges(seqnames=nuc1bidir$chr,ranges=IRanges(start=nuc1bidir$start,end=nuc1bidir$end))  
	bidir_read1 <- summarizeOverlaps(features=feat,reads=sox10_nuc1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
	bidir_read_counts1 <- assay(bidir_read1)
	nuc1bidir <- cbind(nuc1bidir,bidir_read_counts1)
	nuc1bidir <- dplyr::mutate(nuc1bidir, size = end-start)
	write.table(IDasScore,"sox10_Zv9nuc1FlankedWiderBidirRegionsIDisScore.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)
	
	feat <- GRanges(seqnames=intersected1$chr,ranges=IRanges(start=intersected1$start,end=intersected1$end))  
	bidir_read1 <- summarizeOverlaps(features=feat,reads=bam1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
	bidir_read_counts1 <- assay(bidir_read1)
	intersected1 <- cbind(intersected1,bidir_read_counts1)
	intersected1 <- dplyr::mutate(intersected1, size = end-start)
	write.table(intersected1,file="ct711aHets1PrimaryReadsBidirRegions.csv",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t") 
}

intersected1 <- fread("ct711a_Homoz1PrimaryNoncodingUnder10kMergedIntersected.bed",data.table=TRUE,fill=TRUE)
colnames(intersected1) <- c("chr","start","end","ID","score","strand")
intersected2 <- fread("ct711a_Homoz2PrimaryNoncodingUnder10kMergedIntersected.bed",data.table=TRUE,fill=TRUE)
colnames(intersected2) <- c("chr","start","end","ID","score","strand")
hetsintersected1 <- fread("/panfs/qcb-panasas/engie/GRCz11Star/ct711a_150804_hets_nuc1PrimaryUnder10kMergedIntersected.bed",data.table=TRUE,fill=TRUE) #from jan 2021
colnames(hetsintersected1) <- c("chr","start","end","ID","score","strand")
hetsintersected2 <- fread("/panfs/qcb-panasas/engie/GRCz11Star/20_7_21 analyzed/ct711a_150804_hets_nuc2PrimaryUnder10kMergedIntersected.bed",data.table=TRUE,fill=TRUE) #from sept 2020
colnames(hetsintersected2) <- c("chr","start","end","ID","score","strand")

#unstranded
feathomo1 <- GRanges(seqnames=intersected1$chr,ranges=IRanges(start=intersected1$start,end=intersected1$end),strand=intersected1$strand)  #EDIT LATER ADDED STRAND
feathomo2 <- GRanges(seqnames=intersected2$chr,ranges=IRanges(start=intersected2$start,end=intersected2$end),strand=intersected2$strand)  
feathets1 <- GRanges(seqnames=hetsintersected1$chr,ranges=IRanges(start=hetsintersected1$start,end=hetsintersected1$end),strand=hetsintersected1$strand)  
feathets2 <- GRanges(seqnames=hetsintersected2$chr,ranges=IRanges(start=hetsintersected2$start,end=hetsintersected2$end),strand=hetsintersected2$strand)  

bidir_homo1 <- summarizeOverlaps(features=feathomo1,reads=ghomo1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
bidir_homo_counts1 <- assay(bidir_homo1)
bidir_homo2 <- summarizeOverlaps(features=feathomo2,reads=ghomo2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
bidir_homo_counts2 <- assay(bidir_homo2)
bidir_hets1 <- summarizeOverlaps(features=feathets1,reads=ghets1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
bidir_hets_counts1 <- assay(bidir_hets1)
bidir_hets2 <- summarizeOverlaps(features=feathets2,reads=ghets2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
bidir_hets_counts2 <- assay(bidir_hets2)

#plus strand
bidir_homo1plus <- summarizeOverlaps(features=feathomo1,reads=ghomo1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_homo_counts1plus <- assay(bidir_homo1plus)
bidir_homo2plus <- summarizeOverlaps(features=feathomo2,reads=ghomo2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_homo_counts2plus <- assay(bidir_homo2plus)
bidir_hets1plus <- summarizeOverlaps(features=feathets1,reads=ghets1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_hets_counts1plus <- assay(bidir_hets1plus)
bidir_hets2plus <- summarizeOverlaps(features=feathets2,reads=ghets2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_hets_counts2plus <- assay(bidir_hets2plus)

#minus strand
intersected1$strand <- "-"
intersected2$strand <- "-"
hetsintersected1$strand <- "-"
hetsintersected2$strand <- "-"
feathomo1 <- GRanges(seqnames=intersected1$chr,ranges=IRanges(start=intersected1$start,end=intersected1$end),strand=intersected1$strand)  
feathomo2 <- GRanges(seqnames=intersected2$chr,ranges=IRanges(start=intersected2$start,end=intersected2$end),strand=intersected2$strand)  
feathets1 <- GRanges(seqnames=hetsintersected1$chr,ranges=IRanges(start=hetsintersected1$start,end=hetsintersected1$end),strand=hetsintersected1$strand)  
feathets2 <- GRanges(seqnames=hetsintersected2$chr,ranges=IRanges(start=hetsintersected2$start,end=hetsintersected2$end),strand=hetsintersected2$strand)  

bidir_homo1minus <- summarizeOverlaps(features=feathomo1,reads=ghomo1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_homo_counts1minus <- assay(bidir_homo1minus)
bidir_homo2minus <- summarizeOverlaps(features=feathomo2,reads=ghomo2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_homo_counts2minus <- assay(bidir_homo2minus)
bidir_hets1minus <- summarizeOverlaps(features=feathets1,reads=ghets1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_hets_counts1minus <- assay(bidir_hets1minus)
bidir_hets2minus <- summarizeOverlaps(features=feathets2,reads=ghets2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_hets_counts2minus <- assay(bidir_hets2minus)

intersected1counts <- cbind(intersected1,bothstrands=bidir_homo_counts1, plus=bidir_homo_counts1plus, minus=bidir_homo_counts1minus)
intersected1counts <- dplyr::mutate(intersected1counts, size = end-start)
intersected2counts <- cbind(intersected2,bothstrands=bidir_homo_counts2, plus=bidir_homo_counts2plus, minus=bidir_homo_counts2minus)
intersected2counts <- dplyr::mutate(intersected2counts, size = end-start)

hetsintersected1counts <- cbind(hetsintersected1,bothstrands=bidir_hets_counts1,plus=bidir_hets_counts1plus, minus=bidir_hets_counts1minus)
hetsintersected1counts <- dplyr::mutate(hetsintersected1counts, size = end-start)
hetsintersected2counts <- cbind(hetsintersected2,bothstrands=bidir_hets_counts2, plus=bidir_hets_counts2plus, minus=bidir_hets_counts2minus)
hetsintersected2counts <- dplyr::mutate(hetsintersected2counts, size = end-start)

## Keep only regions that are reprodicuble over replicates NEED TO ADAPT
repmerge <- mergeByOverlaps(intersected1,intersected2,ignore.strand=TRUE)
repmergedf <- data.frame(chr=c(seqnames(repmerge$intersected1),seqnames(repmerge$intersected2)),
                        start=c(start(repmerge$intersected1),start(repmerge$intersected2)),
                        end=c(end(repmerge$intersected1),end(repmerge$intersected2)),
                        ID=c(paste("nuc1",1:length(repmerge),sep=""),paste("nuc2",1:length(repmerge),sep="")),
                        score=0,
                        strand=c(strand(repmerge$intersected1),strand(repmerge$intersected2)))
write.table(repmergedf,file="FlankedWiderNoncodingReprodPreUnion.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

bedtools sort -i FlankedWiderNoncodingReprodPreUnion.bed > FlankedWiderNoncodingReprodPreUnionSorted.bed
bedtools merge -i FlankedWiderNoncodingReprodPreUnionSorted.bed -d 0 > FlankedWiderNoncodingReprodOnly.bed

# TFBS alignment
library(BSgenome.Drerio.UCSC.danRer11)
library(motifmatchr)
library(TFBSTools) 

# downloaded single batch .txt for JASPAR code PFMs from JASPAR2020, raw JASPAR format
motifs <- readJASPARMatrix("JASPAR2020CoreTransfac/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt", matrixClass="PFM")
matchscorecalc <- matchMotifs(motifs,gHomo1,genome="danRer11",out="score") #make sure GRanges object is UCSC seqlevelsStyle
matches <- motifMatches(matchscorecalc)
scores <- assay(matchscorecalc)
