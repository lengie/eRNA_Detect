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

loadBAM <- function(file){
     flag <- scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE)
     reads <- readGAlignments(file,param=ScanBamParam(flag=flag))
     return(reads)
}


ncOverlaps <- function(gbam,coding,filename){
	seqlevelsStyle(gbam) <- "UCSC"
	seqlevelsStyle(coding) <- "UCSC"
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
	strand <- fread(file,data.table=TRUE,fill=TRUE)
	colnames(strand) <- c("seqnames","start","end","strand")
	return(GRanges(strand))
}
plus_file <- "SET_THESE_MergedPlus.bed"
minus_file <- "YOUR_DATASET_MergedMinus.bed"
plus <- read_format(plus_file)
minus <- read_format(minus_file)

overlap_format <- function(plus_strand,minus_strand,file,dist=0){
	merge <- mergeByOverlaps(minus_strand,plus_strand,ignore.strand=TRUE)
	merged <- data.frame(chr=c(seqnames(merge$plus_strand),seqnames(merge$minus_strand)),
                        start=c(start(merge$plus_strand),start(merge$minus_strand)),
                        end=c(end(merge$plus_strand),end(merge$minus_strand)),
                        ID=1:nrow(merge),
                        score=1:nrow(merge),
                        strand=c(strand(merge$plus_strand),strand(merge$minus_strand)))
        filename <- paste(file,"OverlapsToMerge.bed",sep="")
	write.table(merged,file=filename,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

    system(sprintf("bedtools sort -i %s > %sOverlapsMergedSorted.bed",filename,file))
    system(paste("bedtools merge -i ",file,"OverlapsMergedSorted.bed -d ",dist," > ",file,"OverlapsMerged.bed",sep=""))
}

overlap_format(plus,minus,"sox10_871_HiSatUnder10k")

# input is two GRanges objects and we are merging the overlapping regions in a strand specific manner
replicateReprod <- function(rep1,rep2,file,dist=0){
    repmergedf <- data.frame(chr=c(seqnames(rep1),seqnames(rep2)),
                        start=c(start(rep1),start(rep2)),
                        end=c(end(rep1),end(rep2)),
                        ID=1:(length(rep1)+length(rep2)),
                        score=1:(length(rep1)+length(rep2)),
                        strand='+') #just a stand in
    filename <- paste(file,"AllBidirReg.bed",sep="")
    write.table(repmergedf,file=filename,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
    
    system(sprintf("bedtools sort -i %s > %s_AllBidirRegSorted.bed",filename,file))
    system(paste("bedtools merge -i ",file,"_AllBidirRegSorted.bed -d ",dist," > ",file,"ReprodOnly.bed",sep=""))
}


readCounts <- function(df,gr){
	feat <- GRanges(seqnames=df$chr,ranges=IRanges(start=df$start,end=df$end),strand=dr$strand)  
	reads <- summarizeOverlaps(features=feat,reads=gr,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
	read_counts <- assay(reads)
        return(read_counts)	
}

twoRepCountsStranded <- function(df,rep1,rep2){
	df$strand='+'
	oneplus <- readCounts(df,rep1)
        twoplus <- readCounts(df,rep2)
	df$strand='-'
	oneminus <- readCounts(df,rep1)
        twominus <- readCounts(df,rep2)
        df <- dplyr::mutate(df,size=end-start)
        df <- cbind(df,oneplus,oneminus,twoplus,twominus)
        return(df)
}


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

intersected1 <- readCounts(intersected1,bam1)
write.table(intersected1,file="ct711aHets1PrimaryReadsBidirRegions.csv",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t") 
