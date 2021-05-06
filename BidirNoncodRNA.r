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

#will need to switch/reconcile
bedtools intersect -a sox10_Zv9nuc1FlankedNoncodingUnder10kMergedPlus6.bed -b sox10_Zv9nuc1FlankedNoncodingUnder10kMergedMinus6.bed -nonamecheck > sox10_Zv9nuc1FlankedNoncodingUnder10kMergedSeparatelyIntersected.bed
bedtools intersect -a sox10_Zv9nuc2FlankedNoncodingUnder10kMergedPlus6.bed -b sox10_Zv9nuc2FlankedNoncodingUnder10kMergedMinus6.bed -nonamecheck > sox10_Zv9nuc2FlankedNoncodingUnder10kMergedSeparatelyIntersected.bed
	
	
 ## Back in R   
	intersected1 <- fread("GRCz11Star/ct711a_150804_hets_nuc1PrimaryUnder10kMergedIntersected.bed",data.table=TRUE,fill=TRUE)
	colnames(intersected1) <- c("chr","start","end","ID","score","strand")
	
	#merge the intersections for the wider
	nuc1merge <- mergeByOverlaps(minus,plus,ignore.strand=TRUE)
	nuc2merge <- mergeByOverlaps(minus2,plus2,ignore.strand=TRUE)
	nuc1merged <- data.frame(chr=c(seqnames(nuc1merge$plus),seqnames(nuc1merge$minus)),
                        start=c(start(nuc1merge$plus),start(nuc1merge$minus)),
                        end=c(end(nuc1merge$plus),end(nuc1merge$minus)),
                        ID="n/a",
                        score=0,
                        strand=c(strand(nuc1merge$plus),strand(nuc1merge$minus)))
nuc2merged <- data.frame(chr=c(seqnames(nuc2merge$plus2),seqnames(nuc2merge$minus2)),
                        start=c(start(nuc2merge$plus2),start(nuc2merge$minus2)),
                        end=c(end(nuc2merge$plus2),end(nuc2merge$minus2)),
                        ID="n/a",
                        score=0,
                        strand=c(strand(nuc2merge$plus2),strand(nuc2merge$minus2)))
write.table(nuc1merged,file="sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsToMerge.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
write.table(nuc2merged,file="sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsToMerge.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
 
# outside of R
bedtools sort -i sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsToMerge.bed > sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsToMergeSorted.bed
bedtools merge -i sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsToMergeSorted.bed -d 0 > sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsMerged.bed

bedtools sort -i sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsToMerge.bed > sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsToMergeSorted.bed
bedtools merge -i sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsToMergeSorted.bed -d 0 > sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsMerged.bed

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


# TFBS alignment
library(BSgenome.Drerio.UCSC.danRer11)
library(motifmatchr)
library(TFBSTools) 

# downloaded single batch .txt for JASPAR code PFMs from JASPAR2020, raw JASPAR format
motifs <- readJASPARMatrix("JASPAR2020CoreTransfac/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt", matrixClass="PFM")
matchscorecalc <- matchMotifs(motifs,gHomo1,genome="danRer11",out="score") #make sure GRanges object is UCSC seqlevelsStyle
matches <- motifMatches(matchscorecalc)
scores <- assay(matchscorecalc)
