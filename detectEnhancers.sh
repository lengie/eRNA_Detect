## detectEnhancers.sh written by Liana Engie last edited 2020/7/21
## run in the folder where your data is and where you want the output files to be created

#!/bin/bash
#SBATCH --ntasks=6
#SBATCH --time=03:00:00

R
library(GenomicAlignments)
library(GenomicFeatures)
library(dplyr)
library(data.table)

bam2 <- readGAlignmentPairs("GRCz11Star/ct711a_150804_hets_nuc2PrimaryReads.bam",param=ScanBamParam(flag=flag))
gbam2 <- GRanges(bam2)
overlaps2 <- findOverlaps(gbam2,coding,ignore.strand=FALSE)
hits2 <- gbam2[-queryHits(overlaps2),]
underten2 <- hits2[width(hits2)<10000]
ten2 <- hits2[width(hits2)>10000)]
write.table(ten2,file="/panfs/qcb-panasas/engie/GRCz11Star/ct711a_Hets2PrimaryNoncodingOver10k.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

underten2df <- data.frame(chr = as.character(seqnames(underten2)),
                    start = start(underten2)-1,
                    end = end(underten2),
                    name = "n/a",
                    score = 0,
                    strand = strand(underten2)
                    )
write.table(underten2df,file="/panfs/qcb-panasas/engie/GRCz11Star/ct711a_150804_hets_nuc2PrimaryNoncodingUnder10k.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
q()
N

bedtools sort -i ct711a_150804_hets_nuc2PrimaryNoncodingUnder10k.bed > ct711a_150804_hets_nuc2PrimaryNoncodingUnder10kSorted.bed
bedtools merge -i ct711a_150804_hets_nuc2PrimaryNoncodingUnder10kSorted.bed -d 0 -s > ct711a_150804_hets_nuc2PrimaryNoncodingUnder10kMerged.bed

R
library(data.table)
library(dplyr)

outsidemerged2 <- fread("GRCz11Star/ct711a_150804_hets_nuc2PrimaryNoncodingUnder10kMerged.bed",data.table=TRUE,fill=TRUE)
plus2 <- dplyr::filter(outsidemerged2,V4=="+")
minus2 <- dplyr::filter(outsidemerged2,V4=="-")

plus62 <- data.frame(plus2$V1,plus2$V2,plus2$V3,ID="n/a",score=0,strand=plus2$V4)
minus62 <- data.frame(minus2$V1,minus2$V2,minus2$V3,ID="n/a",score=0,strand=minus2$V4)
write.table(plus62,file="GRCz11Star/ct711a_150804_hets_nuc2PrimaryUnder10kMergedPlus.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(minus62,file="GRCz11Star/ct711a_150804_hets_nuc2PrimaryUnder10kMergedMinus.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
q()
N

intersectBed -sorted -S -a ct711a_150804_hets_nuc2PrimaryUnder10kMergedMinus.bed -b ct711a_150804_hets_nuc2PrimaryUnder10kMergedPlus.bed > ct711a_150804_hets_nuc2PrimaryUnder10kMergedIntersected.bed

R
library(data.table)
library(dplyr)

intersected2 <- fread("GRCz11Star/ct711a_150804_hets_nuc2PrimaryUnder10kMergedIntersected.bed",data.table=TRUE,fill=TRUE)
colnames(intersected2) <- c("chr","start","end","ID","score","strand")
one2 <- dplyr::filter(intersected2,chr==1)
one2$chr <- paste("chr",one2$chr,sep="")
write.table(one2,file="ct711a_150804_hets_nuc2PrimaryUnder10kMergedIntersectedChr1.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
feat <- GRanges(seqnames=intersected2$chr,ranges=IRanges(start=intersected2$start,end=intersected2$end))  
bidir_read2 <- summarizeOverlaps(features=feat,reads=bam2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
bidir_read_counts2 <- assay(bidir_read2)
intersected2 <- cbind(intersected2,bidir_read_counts2)
intersected2 <- dplyr::mutate(intersected2, size = end-start)
write.table(intersected2,file="GRCz11Star/ct711aHets2PrimaryReadsBidirRegions.csv",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")