### detectEnhancers.R
###
### Purpose: Take RNA-seq data and output locations and FPKMs of regions with bidirectional reads within a given region 
###
###
### Written by Liana Engie
### Last updated: April 2017
###
### detectEnhancers(chromosome, input_start, input_end,strand)
### Input: string chromosome number, 
### Output:
library(Rsamtools)
library(biomaRt)
library(GenomicAlignments)
library(GenomicFeatures)
library(ggplot2)
library(data.table) 

clusters <- findReadRegions(chromosome,input_start,input_end,strand,bed){
	if(strand!="-"||strand!="+"){
		print('Strand should be "-" or "+"')
		break 
	}
	
	interval <- GRanges(seqnames=chromosome,range=IRanges(input_start,input_end),strand=strand)
	list <- subsetByOverlaps(bed,interval) #GRanges obj
	if(length(ranges(list))==0){
		print('No clusters of reads in range.')
	}#clusters <- ranges(list) #IRanges obj
	clusters <- list
	
	#now we have the clusters of reads. Get the read counts
	#Remove read counts that are 0
	#Find read counts of the opposite strand
	#Calculate TPM
	#Report TPM, keep if it's over 2
}

range <- reverseStrand(range){
	anti <- strand(range)
	for(i=1:length(runValue(anti))){
		if(runValues(anti)[i]=="-"){
			runValues(anti)[i] <- "+"
		}else if (runValues(anti)[i]=="+"){runValues(anti)[i] <- "-"} #don't want to mess with the *
	}
	strand(range) <- anti
}

TPM <- TPMCalc(clusters, reads){
	counts <- summarizeOverlaps(features=clusters,reads=reads,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	sense_counts <- assay(counts)
	#just read strand from GRanges obj
	opp_strand <- reverseStrand(clusters)
	
	counts <- summarizeOverlaps(features=opp_strand,reads=reads,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	anti_counts <- assay(counts)
	
	#could we use an lapply instead?
	#maybe if we appended the read counts as metadata in the GRanges obj...
	rpk <- c()
	for(i=1:length(ranges(cluster))){
		rpk[i] <- sense_counts[i]/width(list)[i] #or do I have to use append 
		 
	}
	scaling <- sum(rpk)/1000000 #except this should be for the whole sample
	TPM <- rpk/scaling 
	
}

FPKM <- FPKMCalc(clusters,reads,scaling_factor){
	counts <- summarizeOverlaps(features=clusters,reads=reads,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	sense_counts <- assay(counts)
	#just read strand from GRanges obj
	opp_strand <- reverseStrand(clusters)
	
	counts <- summarizeOverlaps(features=opp_strand,reads=reads,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	anti_counts <- assay(counts)
	
	
	FPKM <- data.frame()
	colnames(FPKM) <- c("start","width","sense","antisense")
}

hets1 <- "/auto/cmb-00/rr/engie/RNA/Aligned.sortedByCoord.out.bam" 
hets2 <- "/auto/cmb-00/rr/engie/RNA/hets2.bam"
hets1frag <- 50.731011 
hets2frag <- 56.315336 
flag <- scanBamFlag(isNotPrimaryRead=FALSE, isDuplicate=FALSE)
hetsread1 <- readGAlignmentPairs(hets1,param=ScanBamParam(flag=flag))
hetsread2 <- readGAlignmentPairs(hets2,param=ScanBamParam(flag=flag))
file <- '/auto/cmb-00/rr/engie/RNA/merged1.bed'
bed1 <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
bed1R <- GRanges(seqnames=bed1$V1,ranges=IRanges(start=bed1$V2, end=bed1$V3),strand=bed1$V4)

clusters <- findReadRegions("chr1",39693870,39800883,"-",bed1R)
MAML_TPM <- TPMCalc(clusters,hetsread1)
MAML_TPM2 <- TPMCalc(clusters,hetsread2)

bach2 <- findReadRegions("chr20",24073750,24223950, "+",bed1R)
bach2FPKM <- FPKMCalc(bach2,hetsread1,hets1frag)
bach2FPKMcheck <- FPKMCalc(bach2,hetsread3,hets2frag)

