### Kim et al 2010 paper uses AB SOLiD Systems for RNA-sequencing, and all data sets are reported with single 5' end basepair
###
### Code for extending the bedgraph of bigwig files

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments) 
library(ggplot2)
library(dplyr) 
library(data.table)
library(rtracklayer)
options(scipen=999)
 

## load the annotated coding regions
codingfile <- fread("MMusculus37_ExonsUTRs500bpFlanking.bed")
colnames(codingfile) <- c("chr","start","end","ID","score","strand")
coding <- GRanges(codingfile)

## load the files
loadbg <- function(bg,strand){
    bg <- mutate(bg,strand)
    colnames(bg) <- c("chr","start","end","score","strand")
    return(bg)
}

plusfile <- "GSM530210_hr0c+.bedgraph"
minusfile <- "GSM530210_hr0c-.bedgraph"
plus14file <- "GSM530214_h0_b2+.bedgraph"
minus14file <- "GSM530214_h0_b2-.bedgraph"

plus <- fread(plusfile)
minus <- fread(minusfile)
plus <- loadbg(plus,"+")
minus <- loadbg(minus,"-")

plus14 <- fread(plus14file)
minus14 <- fread(minus14file)
plus14 <- loadbg(plus14,"+")
minus14 <- loadbg(minus14,"-")

plus12 <- fread("GSM530212_hr6c+.bedgraph")
minus12 <- fread("GSM530212_hr6c-.bedgraph")
plus13 <- fread("GSM530213_hr6_b3+.bedgraph",header=FALSE)
minus13 <- fread("GSM530213_hr6_b3-.bedgraph",header=FALSE)
plus16 <- fread("GSM530216_h6_b2+.bedgraph")
minus16 <- fread("GSM530216_h6_b2-.bedgraph")
plus12 <- loadbg(plus12,"+")
minus12 <- loadbg(minus12,"-")
plus13 <- loadbg(plus13,"+")
minus13 <- loadbg(minus13,"-")
plus16 <- loadbg(plus16,"+")
minus16 <- loadbg(minus16,"-")
# probably should have written the code to just include fread ah well

## expand the reads in the 3' direction to be full fragments
expand <- function(bg,length,strand){
    if(strand=="plus"){
        bg$end <- bg$end + length}
    else if(strand=="minus"){
        bg$start <- bg$start - length}
    else{return("Strand must be 'plus' or 'minus'!")}
    return(bg)
}

plus <- expand(plus,84,"plus")
minus <- expand(minus,84,"minus")
plus12 <- expand(plus12,84,"plus")
...
plus13 <- expand(plus13,84,"plus")
minus13 <- expand(minus13,84,"minus")

# wrote new chrlimit check that uses a preloaded file    
## save the individual expanded strand files
checkbed6 <- function(bg,filename){
    bg <- chrLimitCheckNoInt(bg,limit)
    bg <- mutate(bg,ID=1:nrow(bg))
    bg <- bg[,c("chr","start","end","ID","score","strand")]
    write.table(bg,paste(filename,"_readsextended85.bed",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}
checkbed6(plus,"GSM530210+")
checkbed6(minus,"GSM530210-")
checkbed6(plus12,"GSM530212+")
checkbed6(minus12,"GSM530212-")
checkbed6(plus13,"GSM530213+")
checkbed6(minus13,"GSM530213-")
checkbed6(plus14,"GSM530214+")
checkbed6(minus14,"GSM530214-")
checkbed6(plus16,"GSM530216+")
checkbed6(minus16,"GSM530216-")

## make into GRanges for to remove coding regions
gr10 <- GRanges(rbind(plus,minus))
gr12 <- GRanges(rbind(plus12,minus12))
gr13 <- GRanges(rbind(plus13,minus13))
gr14 <- GRanges(rbind(plus14,minus14))
gr16 <- GRanges(rbind(plus16,minus16))

## Remove coding regions
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
				  score = mcols(underten)$score,
				  strand = strand(underten)
				  )
	#save the file
	write.table(undertendf,file=paste(filename,"Under10k.bed",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
	return(undertendf)
	}

underten10 <- ncOverlaps(gr10,coding,"GSM530210_FlankedNoncoding")
underten12 <- ncOverlaps(gr12,coding,"GSM530212_FlankedNoncoding")
underten13 <- ncOverlaps(gr13,coding,"GSM530213_FlankedNoncoding")
underten14 <- ncOverlaps(gr14,coding,"GSM530214_FlankedNoncoding")
underten16 <- ncOverlaps(gr16,coding,"GSM530216_FlankedNoncoding")
 

## outside of R -- merge the individual strands and sum scores of the merged regions
bedtools sort -i GSM530210_FlankedNoncodingUnder10k.bed > GSM530210_FlankedNoncodingSorted.bed 
bedtools merge -i GSM530210_FlankedNoncodingSorted.bed -d 0 -S - -c 5 -o sum > GSM530210_FlankedNoncodingScorePlusOnly.bed 
bedtools merge -i GSM530210_FlankedNoncodingSorted.bed -d 0 -S + -c 5 -o sum > GSM530210_FlankedNoncodingScoreMinusOnly.bed 

bedtools sort -i GSM530214_FlankedNoncodingUnder10k.bed > GSM530214_FlankedNoncodingSorted.bed 
bedtools merge -i GSM530214_FlankedNoncodingSorted.bed -d 0 -S - -c 5 -o sum > GSM530214_FlankedNoncodingScorePlusOnly.bed 
bedtools merge -i GSM530214_FlankedNoncodingSorted.bed -d 0 -S + -c 5 -o sum > GSM530214_FlankedNoncodingScoreMinusOnly.bed 

bedtools sort -i GSM530212_FlankedNoncodingUnder10k.bed > GSM530212_FlankedNoncodingSorted.bed 
bedtools merge -i GSM530212_FlankedNoncodingSorted.bed -d 0 -S - -c 5 -o sum > GSM530212_FlankedNoncodingScorePlusOnly.bed 
bedtools merge -i GSM530212_FlankedNoncodingSorted.bed -d 0 -S + -c 5 -o sum > GSM530212_FlankedNoncodingScoreMinusOnly.bed 

bedtools sort -i GSM530213_FlankedNoncodingUnder10k.bed > GSM530213_FlankedNoncodingSorted.bed 
bedtools merge -i GSM530213_FlankedNoncodingSorted.bed -d 0 -S - -c 5 -o sum > GSM530213_FlankedNoncodingScorePlusOnly.bed 
bedtools merge -i GSM530213_FlankedNoncodingSorted.bed -d 0 -S + -c 5 -o sum > GSM530213_FlankedNoncodingScoreMinusOnly.bed 
 
bedtools sort -i GSM530216_FlankedNoncodingUnder10k.bed > GSM530216_FlankedNoncodingSorted.bed 
bedtools merge -i GSM530216_FlankedNoncodingSorted.bed -d 0 -S - -c 5 -o sum > GSM530216_FlankedNoncodingScorePlusOnly.bed 
bedtools merge -i GSM530216_FlankedNoncodingSorted.bed -d 0 -S + -c 5 -o sum > GSM530216_FlankedNoncodingScoreMinusOnly.bed 

## move back into R to reformat so the strands can be further intersected
read_format_score <- function(file,str){
	strand <- fread(file,data.table=TRUE,fill=TRUE)
	colnames(strand) <- c("seqnames","start","end","score")
        strand <- mutate(strand,strand=str)
	return(GRanges(strand))
}
plus <- read_format_score("GSM530210_FlankedNoncodingScorePlusOnly.bed","+")
minus <- read_format_score("GSM530210_FlankedNoncodingScoreMinusOnly.bed","-")
plus14 <- read_format_score("GSM530214_FlankedNoncodingScorePlusOnly.bed","+")
minus14 <- read_format_score("GSM530214_FlankedNoncodingScoreMinusOnly.bed","-")
plus12 <- read_format_score("GSM530212_FlankedNoncodingScorePlusOnly.bed","+")
minus12 <- read_format_score("GSM530212_FlankedNoncodingScoreMinusOnly.bed","-")
plus13 <- read_format_score("GSM530213_FlankedNoncodingScorePlusOnly.bed","+")
minus13 <- read_format_score("GSM530213_FlankedNoncodingScoreMinusOnly.bed","-")
plus16 <- read_format_score("GSM530216_FlankedNoncodingScorePlusOnly.bed","+")
minus16 <- read_format_score("GSM530216_FlankedNoncodingScoreMinusOnly.bed","-")


overlap_format <- function(plus_strand,minus_strand){
	merge <- mergeByOverlaps(minus_strand,plus_strand,ignore.strand=TRUE)
	merged <- data.frame(chr=c(seqnames(merge$plus_strand),seqnames(merge$minus_strand)),
                        start=c(start(merge$plus_strand),start(merge$minus_strand)),
                        end=c(end(merge$plus_strand),end(merge$minus_strand)),
                        ID=1:(length(merge$plus_strand) + length(merge$minus_strand)),
                        score=c(-merge[,2],merge[,4]),
                        strand=c(strand(merge$plus_strand),strand(merge$minus_strand)))
	return(merged)
}
rep1merged <- overlap_format(plus,minus)
rep12merged <- overlap_format(plus12,minus12)
rep13merged <- overlap_format(plus13,minus13)
rep14merged <- overlap_format(plus14,minus14)
rep16merged <- overlap_format(plus16,minus16)
write.table(rep1merged,file="GSM530210_FlankedNoncodingScoreOverlapsToMerge.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
write.table(rep14merged,file="GSM530214_FlankedNoncodingScoreOverlapsToMerge.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
write.table(rep12merged,file="GSM530212_FlankedNoncodingScoreOverlapsToMerge.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
write.table(rep13merged,file="GSM530213_FlankedNoncodingScoreOverlapsToMerge.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
write.table(rep16merged,file="GSM530216_FlankedNoncodingScoreOverlapsToMerge.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

# outside of R do another merge, between the overlapping regions of each strand
bedtools sort -i GSM530210_FlankedNoncodingScoreOverlapsToMerge.bed > GSM530210_FlankedNoncodingScoreOverlapsToMergeSorted.bed
bedtools merge -i GSM530210_FlankedNoncodingScoreOverlapsToMergeSorted.bed -d 0 -c 5 -o sum > GSM530210_FlankedNoncodingScoreOverlapsMerged.bed

bedtools sort -i GSM530212_FlankedNoncodingScoreOverlapsToMerge.bed > GSM530212_FlankedNoncodingScoreOverlapsToMergeSorted.bed
bedtools merge -i GSM530212_FlankedNoncodingScoreOverlapsToMergeSorted.bed -d 0 -c 5 -o sum > GSM530212_FlankedNoncodingScoreOverlapsMerged.bed

bedtools sort -i GSM530213_FlankedNoncodingScoreOverlapsToMerge.bed > GSM530213_FlankedNoncodingScoreOverlapsToMergeSorted.bed
bedtools merge -i GSM530213_FlankedNoncodingScoreOverlapsToMergeSorted.bed -d 0 -c 5 -o sum > GSM530213_FlankedNoncodingScoreOverlapsMerged.bed

bedtools sort -i GSM530214_FlankedNoncodingScoreOverlapsToMerge.bed > GSM530214_FlankedNoncodingScoreOverlapsToMergeSorted.bed
bedtools merge -i GSM530214_FlankedNoncodingScoreOverlapsToMergeSorted.bed -d 0 -c 5 -o sum > GSM530214_FlankedNoncodingScoreOverlapsMerged.bed

bedtools sort -i GSM530216_FlankedNoncodingScoreOverlapsToMerge.bed > GSM530216_FlankedNoncodingScoreOverlapsToMergeSorted.bed
bedtools merge -i GSM530216_FlankedNoncodingScoreOverlapsToMergeSorted.bed -d 0 -c 5 -o sum > GSM530216_FlankedNoncodingScoreOverlapsMerged.bed

 
 
## extending the bigwig files directly 
extendBigwig <- function(bw,step,strand){
    dt <- data.table()
    if(strand=='-'){
        for(i in 1:length(bw)){
            dt <- rbind(dt,list(chr=as.character(seqnames(bw)[i]),start=start(bw)[i],end=end(bw)[i],score=mcols(bw)$score[i]))
            for(j in 1:step-1){
                dt <- rbind(dt,list(chr=as.character(seqnames(bw)[i]),start=start(bw)[i]+j,end=start(bw)[i]+j,score=mcols(bw)$score[i]))
            }
        }
    }else if(strand=='+'){
        for(i in 1:length(bw)){
            dt <- rbind(dt,data.table(chr=as.character(seqnames(bw)[i]),start=start(bw)[i],end=end(bw)[i],score=mcols(bw)$score[i]))
            for(j in 1:step-1){
                dt <- rbind(dt,list(chr=as.character(seqnames(bw)[i]),start=start(bw)[i]-j,end=start(bw)[i]-j,score=mcols(bw)$score[i]))
            }
        }
    }else{return("strand must either be '-' or '+'")}
    dt <- dt[, lapply(score, sum), by=list(chr, start, end)]
    gr <- GRanges(as.data.frame(dt))
}
